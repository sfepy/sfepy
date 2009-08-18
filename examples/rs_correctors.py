# c: 05.05.2008, r: 05.05.2008
import sys
sys.path.append( '.' )

from sfepy.fem.periodic import *
from sfepy.homogenization.utils import define_box_regions

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

    # Sides and corners.
    regions.update( define_box_regions( 2, (wx, wy) ) )
    return is3d, regions, mat_ids

##
# c: 05.05.2008, r: 05.05.2008
def get_pars(ts, coor, mode=None, region=None, ig=None, mat_ids = []):
    """Define material parameters:
         $D_ijkl$ (elasticity),
       in a given region."""
    if mode == 'qp':
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

        for key, val in out.iteritems():
            out[key] = nm.tile(val, (coor.shape[0], 1, 1))

        return out
    
##
# Mesh file.
filename_mesh = 'osteonT1_11.mesh'

##
# Define regions (subdomains, boundaries) - $Y$, $Y_i$, ...
# depending on a mesh used.
is3d, regions, mat_ids = define_regions( filename_mesh )

functions = {
    'get_pars' : (lambda ts, coors, mode=None, region=None, ig=None:
                  get_pars(ts, coors, mode, region, ig, mat_ids=mat_ids),),
    'match_x_plane' : (match_x_plane,),
    'match_y_plane' : (match_y_plane,),
    'match_z_plane' : (match_z_plane,),
    'match_x_line' : (match_x_line,),
    'match_y_line' : (match_y_line,),
}

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
    'region' : 'Y',
    'function' : 'get_pars',
}

##
# Numerical quadratures for volume (i3 - order 3) integral terms.
integral_1 = {
    'name' : 'i3',
    'kind' : 'v',
    'quadrature' : 'gauss_o3_d%d' % dim,
}

##
# Homogenized coefficients to compute.
coefs = {
    'E' : {'requires' : ['pis', 'corrs_rs'],
           'variables' : ['Pi1', 'Pi2'],
           'expression' : 'dw_lin_elastic.i3.Y( m.D, Pi1, Pi2 )'},
}

##
# Data required to compute the homogenized coefficients.
all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:dim] ]
requirements = {
    'pis' : {
        'variables' : ['uc'],
    },
    ##
    # Steady state correctors $\bar{\omega}^{rs}$.
    'corrs_rs' : {
         'requires' : ['pis'],
         'variables' : ['uc', 'vc', 'Pi'],
         'ebcs' : ['fixed_u'],
         'epbcs' : all_periodic,
         'equations' : {'eq' : """dw_lin_elastic.i3.Y( m.D, vc, uc )
                                = - dw_lin_elastic.i3.Y( m.D, vc, Pi )"""},
         'save_name' : 'corrs_elastic',
    },
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
    'kind' : 'ls.scipy_direct', # Direct solver.
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
# c: 05.05.2008, r: 28.11.2008
def main():
    import os
    from sfepy.base.base import spause
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.fem import eval_term_op, ProblemDefinition
    import sfepy.homogenization.coefs_elastic as ce

    nm.set_printoptions( precision = 3 )

    spause( r""">>>
First, this file will be read in place of an input
(problem description) file.
Press 'q' to quit the example, press any other key to continue...""" )
    required, other = get_standard_keywords()
    required.remove( 'equations' )
    # Use this file as the input file.
    conf = ProblemConf.from_file( __file__, required, other )
    print conf.to_dict().keys()
    spause( r""">>>
...the read input as a dict (keys only for brevity).
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
    req = conf.requirements['pis']
    mini_app = ce.ShapeDimDim( 'pis', problem, req )
    pis = mini_app()
    print pis
    spause( r""">>>
...the $\Pi$ operators.
['q'/other key to quit/continue...]""" )

    spause( r""">>>
Next, $E_{ijkl}$ needs so called steady state correctors $\bar{\omega}^{rs}$,
computed now.
['q'/other key to quit/continue...]""" )
    req = conf.requirements['corrs_rs']

    save_name = req.get( 'save_name', '' )
    cwd = os.path.split( os.path.join( os.getcwd(), __file__ ) )[0]
    name = os.path.join( cwd, save_name )

    mini_app = ce.CorrectorsElasticRS( 'steady rs correctors', problem, req )
    save_hook = mini_app.make_save_hook( name, problem.output_format )
    corrs_rs = mini_app( data = {'pis': pis}, save_hook = save_hook )
    print corrs_rs
    spause( r""">>>
...the $\bar{\omega}^{rs}$ correctors.
The results are saved in: %s*.%s
['q'/other key to quit/continue...]""" % (name, problem.output_format) )

    spause( r""">>>
Then the volume of the domain is needed.
['q'/other key to quit/continue...]""" )
    volume = eval_term_op( None, 'd_volume.i3.Y( uc )', problem )
    print volume

    spause( r""">>>
...the volume.
['q'/other key to quit/continue...]""" )

    spause( r""">>>
Finally, $E_{ijkl}$ can be computed.
['q'/other key to quit/continue...]""" )
    mini_app = ce.ElasticCoef( 'homogenized elastic tensor',
                               problem, conf.coefs['E'] )
    c_e = mini_app( volume, data = {'pis': pis, 'corrs_rs' : corrs_rs} )
    print r""">>>
The homogenized elastic coefficient $E_{ijkl}$, symmetric storage
with rows, columns in 11, 22, 12 ordering:"""
    print c_e
    
if __name__ == '__main__':
    main()
