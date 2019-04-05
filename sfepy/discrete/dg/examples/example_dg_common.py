import numpy as nm

from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)

from sfepy.terms import register_term
from sfepy.solvers import  register_solver

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

def get_cfl_setup(CFL):

    def setup_cfl_condition(problem):
        """
        Sets up CFL condition for problem ts_conf in problem
        :param problem: discrete.problem.Problem
        :return:
        """
        ts_conf = problem.ts_conf
        mesh = problem.domain.mesh
        dim = mesh.dim
        first_field = list(problem.fields.values())[0]
        first_field_name = list(problem.fields.keys())[0]
        approx_order = first_field.approx_order
        mats = problem.create_materials('a')
        velo = problem.conf_materials['material_a__0'].values["val"]
        max_velo = nm.max(nm.linalg.norm(velo))
        dx = nm.min(problem.domain.mesh.cmesh.get_volumes(dim))
        order_corr = 1. / (2 * approx_order + 1)
        dt = dx / max_velo * CFL * order_corr
        # time_steps_N = int((tf - t0) / dt) * 2
        tn = int(nm.ceil((ts_conf.t1 - ts_conf.t0) / dt))
        dtdx = dt / dx

        ts_conf += Struct(dt=dt, n_step=tn)
        output("Preprocessing hook setup_cfl_condition:")
        output("Approximation order of field {}({}) is {}".format(first_field_name, first_field.family_name, approx_order))
        output("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
        output("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
        output("CFL coefficient was {0} and order correction 1/{1} = {2}".format(CFL,  (2 * approx_order + 1), order_corr))
        output("Courant number c = max(norm(a)) * dt/dx = {0}".format(max_velo * dtdx))
        output("------------------------------------------")
        output("Time stepping solver is {}".format(ts_conf.kind))


    return setup_cfl_condition

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
