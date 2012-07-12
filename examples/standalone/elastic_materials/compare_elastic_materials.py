"""
Compare various elastic materials w.r.t. uniaxial tension/compression test.

Requires Matplotlib.
"""
from optparse import OptionParser
import sys
sys.path.append( '.' )

import numpy as nm

def define():
    """Define the problem to solve."""

    filename_mesh = 'el3.mesh'

    options = {
        'nls' : 'newton',
        'ls' : 'ls',
        'ts' : 'ts',
        'save_steps' : -1,
    }

    functions = {
        'linear_tension' : (linear_tension,),
        'linear_compression' : (linear_compression,),
        'empty' : (lambda ts, coor, mode, region, ig: None,),
    }

    field_1 = {
        'name' : 'displacement',
        'dtype' : nm.float64,
        'shape' : (3,),
        'region' : 'Omega',
        'approx_order' : 1,
    }

    # Coefficients are chosen so that the tangent stiffness is the same for all
    # material for zero strains.
    # Young modulus = 10 kPa, Poisson's ratio = 0.3
    material_1 = {
        'name' : 'solid',

        'values' : {
            'K'  : 8.333, # bulk modulus
            'mu_nh' : 3.846, # shear modulus of neoHookean term
            'mu_mr' : 1.923, # shear modulus of Mooney-Rivlin term
            'kappa' : 1.923, # second modulus of Mooney-Rivlin term
            'lam' : 5.769, # Lame coefficients for LE term
            'mu_le' : 3.846,
        }
    }

    material_2 = {
        'name' : 'load',
        'function' : 'empty'
    }

    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    regions = {
        'Omega' : ('all', {}),
        'Bottom' : ('nodes in (z < 0.1)', {}),
        'Top' : ('nodes in (z > 2.9)', {}),
    }

    ebcs = {
        'fixb' : ('Bottom', {'u.all' : 0.0}),
        'fixt' : ('Top', {'u.[0,1]' : 0.0}),
    }

    ##
    # Balance of forces.
    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'order' : 1,
    }
    integral_3 = {
        'name' : 'isurf',
        'kind' : 's',
        'order' : 2,
    }
    equations = {
        'linear' : """dw_lin_elastic_iso.i1.Omega( solid.lam, solid.mu_le, v, u )
                      = dw_surface_ltr.isurf.Top( load.val, v )""",
        'neoHookean' : """dw_tl_he_neohook.i1.Omega( solid.mu_nh, v, u )
                        + dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )
                        = dw_surface_ltr.isurf.Top( load.val, v )""",
        'Mooney-Rivlin' : """dw_tl_he_neohook.i1.Omega( solid.mu_mr, v, u )
                           + dw_tl_he_mooney_rivlin.i1.Omega( solid.kappa, v, u )
                           + dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )
                           = dw_surface_ltr.isurf.Top( load.val, v )""",
    }

    ##
    # Solvers etc.
    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 5,
        'eps_a'      : 1e-10,
        'eps_r'      : 1.0,
        'macheps'    : 1e-16,
        'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp': 0.001,
        'ls_on'      : 1.1,
        'ls_min'     : 1e-5,
        'check'      : 0,
        'delta'      : 1e-6,
        'is_plot'    : False,
        'problem'    : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
    }

    solver_2 = {
        'name' : 'ts',
        'kind' : 'ts.simple',

        't0'    : 0,
        't1'    : 1,
        'dt'    : None,
        'n_step' : 101, # has precedence over dt!
    }

    return locals()

##
# Pressure tractions.
def linear_tension(ts, coor, mode=None, **kwargs):
    if mode == 'qp':
        val = nm.tile(0.1 * ts.step, (coor.shape[0], 1, 1))
        return {'val' : val}

def linear_compression(ts, coor, mode=None, **kwargs):
    if mode == 'qp':
        val = nm.tile(-0.1 * ts.step, (coor.shape[0], 1, 1))
        return {'val' : val}


def store_top_u( displacements ):
    """Function _store() will be called at the end of each loading step. Top
    displacements will be stored into `displacements`."""
    def _store( problem, ts, state ):

        top = problem.domain.regions['Top']
        top_u = problem.get_variables()['u'].get_state_in_region( top )
        displacements.append( nm.mean( top_u[:,-1] ) )

    return _store

def solve_branch(problem, branch_function):
    from sfepy.applications import solve_evolutionary

    displacements = {}
    for key, eq in problem.conf.equations.iteritems():
        problem.set_equations( {key : eq} )

        load = problem.get_materials()['load']
        load.set_function(branch_function)

        out = []
        solve_evolutionary(problem, save_results=False,
                           step_hook=store_top_u(out))
        displacements[key] = nm.array( out, dtype = nm.float64 )
    return displacements

usage = """%prog [options]"""
helps = {
    'no_plot' : 'do not show plot window',
}

def main():
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.fem import ProblemDefinition
    from sfepy.base.plotutils import plt

    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-n', '--no-plot',
                      action="store_true", dest='no_plot',
                      default=False, help=helps['no_plot'])
    options, args = parser.parse_args()

    required, other = get_standard_keywords()
    # Use this file as the input file.
    conf = ProblemConf.from_file( __file__, required, other )

    # Create problem instance, but do not set equations.
    problem = ProblemDefinition.from_conf( conf,
                                           init_equations = False )

    # Solve the problem. Output is ignored, results stored by using the
    # step_hook.
    u_t = solve_branch(problem, linear_tension)
    u_c = solve_branch(problem, linear_compression)

    # Get pressure load by calling linear_*() for each time step.
    ts = problem.get_timestepper()
    load_t = nm.array([linear_tension(ts, nm.array([[0.0]]), 'qp')['val']
                       for aux in ts.iter_from( 0 )],
                      dtype=nm.float64).squeeze()
    load_c = nm.array([linear_compression(ts, nm.array([[0.0]]), 'qp')['val']
                       for aux in ts.iter_from( 0 )],
                      dtype=nm.float64).squeeze()

    # Join the branches.
    displacements = {}
    for key in u_t.keys():
        displacements[key] = nm.r_[u_c[key][::-1], u_t[key]]
    load = nm.r_[load_c[::-1], load_t]


    if plt is None:
        print 'matplotlib cannot be imported, printing raw data!'
        print displacements
        print load
    else:
        legend = []
        for key, val in displacements.iteritems():
            plt.plot( load, val )
            legend.append( key )

        plt.legend( legend, loc = 2 )
        plt.xlabel( 'tension [kPa]' )
        plt.ylabel( 'displacement [mm]' )
        plt.grid( True )

        plt.gcf().savefig( 'pressure_displacement.png' )

        if not options.no_plot:
            plt.show()

if __name__ == '__main__':
    main()
