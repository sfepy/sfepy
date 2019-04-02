import numpy as nm

from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

from sfepy.terms import register_term
from sfepy.solvers import  register_solver
from sfepy.discrete.dg.examples.run_dg_examples import get_cfl_setup


# import various ICs
from sfepy.discrete.dg.my_utils.inits_consts import ghump, gsmooth, left_par_q, left_cos, superic, three_step_u, \
    sawtooth_q, const_q, quadr_cub

from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm


# import TSSs
from sfepy.discrete.dg.dg_tssolver import TVDRK3StepSolver, RK4StepSolver, EulerStepSolver
from sfepy.discrete.dg.dg_limiters import IdentityLimiter, Moment1DLimiter

register_term(AdvectDGFluxTerm)
register_solver(TVDRK3StepSolver)
register_solver(RK4StepSolver)
register_solver(EulerStepSolver)


def get_common(approx_order, CFL, t0, t1, limiter, get_ic):

    regions = {
        'Omega' : 'all',
    }

    fields = {
        'density' : ('real', 'scalar', 'Omega', '1d', 'DG', 'legendre') #
    }

    variables = {
        'u' : ('unknown field', 'density', 0, 1),
        'v' : ('test field',    'density', 'u'),
    }

    functions = {
        'get_ic' : (get_ic,)
    }

    ics = {
        'ic' : ('Omega', {'u.0' : 'get_ic'}),
    }

    equations = {
        'Advection' : """
                       dw_volume_dot.i.Omega(v, u)
                       + dw_s_dot_mgrad_s.i.Omega(a.val, u[-1], v)
                       - dw_dg_advect_laxfrie_flux.i.Omega(a.val, v, u[-1]) = 0
                      """
    }

    solvers = {
        "tss" : ('ts.tvd_runge_kutta_3',
                             {"t0": t0,
                              "t1": t1,
                              'limiter' : IdentityLimiter}),
        'nls' : ('nls.newton',{} ),
        'ls'  : ('ls.scipy_direct', {})
    }

    options = {
        'ts' : 'tss',
        'nls' : 'newton',
        'ls' : 'ls',
        'save_times' : 100,
        'pre_process_hook' : get_cfl_setup(CFL)
    }
    return locals()
