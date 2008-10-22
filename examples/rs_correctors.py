# c: 05.05.2008, r: 05.05.2008
import sys
sys.path.append( '.' )

from sfepy.fem.periodic import *

# c: 05.05.2008, r: 05.05.2008
def define_regions( filename ):
    """Define various subdomain for a given mesh file. This function is called
    below."""
    regions = {}
    is3d = False
    
    regions['Y'] = ('all', {})

    eog = 'elements of group %d'
    if filename.find( 'osteonT1' ) >= 0:
        mat_ids = [11, 39, 6, 8, 27, 28, 9, 2, 4, 14, 12, 17, 45, 28, 15]
        regions['Ym'] = (' +e '.join( (eog % im) for im in  mat_ids ), {})
        wx = 0.865
        wy = 0.499

    regions['Yc'] = ('r.Y -e r.Ym', {})

    # Sides.
    regions['Left'] = ('nodes in (x < -%.3f)' % wx, {})
    regions['Right'] = ('nodes in (x > %.3f)' % wx, {})
    regions['Bottom'] = ('nodes in (y < -%.3f)' % wy, {})
    regions['Top'] = ('nodes in (y > %.3f)' % wy, {})
    regions['Corners'] = ("""nodes in
                            ((x < -%.3f) & (y < -%.3f))
                          | ((x >  %.3f) & (y < -%.3f))
                          | ((x >  %.3f) & (y >  %.3f))
                          | ((x < -%.3f) & (y >  %.3f))
                          """ % ((wx, wy) * 4), {})
    return is3d, regions, mat_ids

##
# c: 05.05.2008, r: 05.05.2008
def get_pars( ts, coor, region, ig, mat_ids = [] ):
    """Define material parameters:
         $D_ijkl$ (elasticity),
       in a given region."""
    dim = coor.shape[1]
    sym = (dim + 1) * dim / 2

    m2i = region.domain.mat_ids_to_i_gs
    matrix_igs = [m2i[im] for im in mat_ids]

    out = {}

    # in 1e+10 [Pa]
    lam = 1.7
    mu = 0.3
    o = nm.array( [1.] * dim + [0.] * (sym - dim), dtype = nm.float64 )
    oot = nm.outer( o, o )
    out['D'] = lam * oot + mu * nm.diag( o + 1.0 )

    if ig not in matrix_igs: # channels
        out['D'] *= 1e-1

    return out
    
##
# Mesh file.
filename_mesh = 'examples/osteonT1_11.mesh'

##
# Define regions (subdomains, boundaries) - $Y$, $Y_i$, ...
# depending on a mesh used.
is3d, regions, mat_ids = define_regions( filename_mesh )

if is3d:
    dim, geom = 3, '3_4'
else:
    dim, geom = 2, '2_3'

##
# Define fields: 'displacement' in $Y$,
# 'pressure_m' in $Y_m$.
field_1 = {
    'name' : 'displacement',
    'dim' : (dim,1),
    'domain' : 'Y',
    'bases' : {'Y' : '%s_P1' % geom}
}

field_2 = {
    'name' : 'pressure_m',
    'dim' : (1,1),
    'domain' : 'Ym',
    'bases' : {'Ym' : '%s_P1' % geom}
}

##
# Define corrector variables: unknown displaements: uc, test: vc
# displacement-like variables: Pi, Pi1, Pi2
variables = {
    'uc'       : ('unknown field',   'displacement', 0),
    'vc'       : ('test field',      'displacement', 'uc'),
    'Pi'       : ('parameter field', 'displacement', 'uc'),
    'Pi1'      : ('parameter field', 'displacement', 'uc'),
    'Pi2'      : ('parameter field', 'displacement', 'uc'),
}

##
# Periodic boundary conditions.
if dim == 3:
    epbc_10 = {
        'name' : 'periodic_x',
        'region' : ['Left', 'Right'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'match_x_plane',
    }
    epbc_11 = {
        'name' : 'periodic_y',
        'region' : ['Near', 'Far'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'match_y_plane',
    }
    epbc_12 = {
        'name' : 'periodic_z',
        'region' : ['Top', 'Bottom'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'match_z_plane',
    }
else:
    epbc_10 = {
        'name' : 'periodic_x',
        'region' : ['Left', 'Right'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'match_y_line',
    }
    epbc_11 = {
        'name' : 'periodic_y',
        'region' : ['Top', 'Bottom'],
        'dofs' : {'uc.all' : 'uc.all'},
        'match' : 'match_x_line',
    }
    
##
# Dirichlet boundary conditions.
ebcs = {
    'fixed_u' : ('Corners', {'uc.all' : 0.0}),
}

##
# Material defining constitutive parameters of the microproblem.
material_1 = {
    'name' : 'm',
    'mode' : 'function',
    'region' : 'Y',
    'function' : 'get_pars',
    'extra_args' : {'mat_ids' : mat_ids},
}

##
# Numerical quadratures for volume (i3 - order 3) integral terms.
integral_1 = {
    'name' : 'i3',
    'kind' : 'v',
    'quadrature' : 'gauss_o3_d%d' % dim,
}

##
# Steady state correctors $\bar{\omega}^{rs}$.
equations = {
    'eq_1' : 
    """dw_lin_elastic.i3.Y( m.D, vc, uc )
       = - dw_lin_elastic.i3.Y( m.D, vc, Pi )""",
}

##
# FE assembling options.
fe = {
    'chunk_size' : 100000,
    'cache_override' : True,
}

##
# Solvers.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack', # Direct solver.
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 2,
    'eps_a'      : 1e-8,
    'eps_r'      : 1e-2,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 0.99999,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

############################################
# Mini-application below, computing the homogenized elastic coefficients.

##
# c: 28.02.2007, r: 13.02.2008
def build_op_pi( var_name, problem, ir, ic ):
    """\Pi^{rs}_i = \delta_{ri} y_s. """
    var = problem.variables[var_name]
    coor = var.field.get_coor()[:,:-1]

    pi = nm.zeros_like( coor )
    pi[:,ir] = coor[:,ic]
    pi.shape = (pi.shape[0] * pi.shape[1],)

    return pi

##
# c: 05.05.2008, r: 05.05.2008
def create_pis( problem, variables, var_name ):
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    pis = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = build_op_pi( var_name, problem, ir, ic )
            pis[ir,ic] = pi
    return pis

##
# c: 05.05.2008, r: 05.05.2008
def  solve_steady_correctors_rs( problem, equations, variables, pis,
                               ofn_trunk, post_process_hook = None,
                               file_per_var = False ):
    """Compute the steady state correctors $\bar{\omega}^{rs}$"""
    from sfepy.base.base import Struct
    
    dim = problem.domain.mesh.dim

    problem.set_variables( variables )
    problem.set_equations( equations )

    problem.time_update()

    states_rs = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = pis[ir,ic]
            # Non-state variables must be assigned manually.
            problem.variables['Pi'].data_from_data( pi )

            state = problem.create_state_vector()
            problem.apply_ebc( state )
            state = problem.solve()
            assert problem.variables.has_ebc( state )
            states_rs[ir,ic] = state

            problem.save_state( ofn_trunk + '_steady_rs_%d%d.vtk' % (ir, ic),
                               state, post_process_hook = post_process_hook,
                               file_per_var = file_per_var )
    return Struct( name = 'Steady RS correctors',
                   states_rs = states_rs,
                   di = problem.variables.di )

##
# c: 05.03.2008, r: 05.03.2008
def iter_sym( dim ):
    for ii in xrange( dim ):
        yield ii, ii
    for ir in xrange( 0, dim ):
        for ic in xrange( ir + 1, dim ):
            yield ir, ic

def coef_e( problem, corrs_rs, pis ):
    """Homogenized elastic coefficient $E_{ijkl}$."""
    from sfepy.fem import eval_term_op

    coef_term = 'dw_lin_elastic.i3.Y( m.D, Pi1, Pi2 )'

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2
    coef = nm.zeros( (sym, sym), dtype = nm.float64 )

    indx = corrs_rs.di.indx['uc']
    for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
        omega1 = corrs_rs.states_rs[irr,icr][indx]
        pi1 = pis[irr,icr] + omega1
        # Non-state variables must be assigned manually.
        problem.variables['Pi1'].data_from_data( pi1 )
            
        for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
            omega2 = corrs_rs.states_rs[irc,icc][indx]
            pi2 = pis[irc,icc] + omega2
            # Non-state variables must be assigned manually.
            problem.variables['Pi2'].data_from_data( pi2 )

            # Variables have their data, so evaluate the term.
            val = eval_term_op( None, coef_term, problem, call_mode = 'd_eval' )

            coef[ir,ic] = val
    return coef

##
# c: 05.05.2008, r: 05.05.2008
def main():
    from sfepy.base.base import spause
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.fem import ProblemDefinition
    from sfepy.base.ioutils import get_trunk

    nm.set_printoptions( precision = 3 )

    spause( r""">>>
First, this file will be read in place of an input
(problem description) file.
Press 'q' to quit the example, press any other key to continue...""" )
    required, other = get_standard_keywords()
    # Use this file as the input file.
    conf = ProblemConf.from_file( __file__, required, other )
    print conf
    spause( r""">>>
...the read input.
['q'/other key to quit/continue...]""" )

    spause( r""">>>
Now the input will be used to create a ProblemDefinition instance.
['q'/other key to quit/continue...]""" )
    problem = ProblemDefinition.from_conf( conf,
                                          init_variables = False,
                                          init_equations = False )
    print problem
    spause( r""">>>
...the ProblemDefinition instance.
['q'/other key to quit/continue...]""" )


    spause( r""">>>
The homogenized elastic coefficient $E_{ijkl}$ is expressed
using $\Pi$ operators, computed now. In fact, those operators are permuted
coordinates of the mesh nodes.
['q'/other key to quit/continue...]""" )
    pis = create_pis( problem, conf.variables, 'Pi' )
    print pis
    spause( r""">>>
...the $\Pi$ operators.
['q'/other key to quit/continue...]""" )

    ofn_trunk = get_trunk( conf.filename_mesh ) + '_out'
    spause( r""">>>
Next, $E_{ijkl}$ needs so called steady state correctors $\bar{\omega}^{rs}$,
computed now. The results will be saved in: %s_*.vtk
['q'/other key to quit/continue...]""" % ofn_trunk )

    corrs_rs = solve_steady_correctors_rs( problem, conf.equations,
                                        conf.variables, pis, ofn_trunk )
    print corrs_rs
    spause( r""">>>
...the $\bar{\omega}^{rs}$ correctors.
['q'/other key to quit/continue...]""" )


    spause( r""">>>
Finally, $E_{ijkl}$ can be computed.
['q'/other key to quit/continue...]""" )
    c_e = coef_e( problem, corrs_rs, pis )
    print r""">>>
The homogenized elastic coefficient $E_{ijkl}$, symmetric storage
with rows, columns in 11, 22, 12 ordering:"""
    print c_e
    
if __name__ == '__main__':
    main()
