import numpy as nm

from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)

from sfepy.terms import register_term
from sfepy.solvers import register_solver

# import various ICs
from sfepy.discrete.dg.my_utils.inits_consts import ghump, gsmooth, left_par_q, left_cos, superic, three_step_u, \
    sawtooth_q, const_q, quadr_cub


from sfepy.discrete.dg.dg_tssolver import TVDRK3StepSolver, RK4StepSolver, EulerStepSolver

from sfepy.discrete.dg.dg_limiters import IdentityLimiter, MomentLimiter1D

from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm, NonlinScalarDotGradTerm, NonlinearHyperDGFluxTerm
from sfepy.discrete.dg.dg_terms import DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

register_term(AdvectDGFluxTerm)
register_term(NonlinScalarDotGradTerm)
register_term(NonlinearHyperDGFluxTerm)
register_term(DiffusionDGFluxTerm)
register_term(DiffusionInteriorPenaltyTerm)

register_solver(TVDRK3StepSolver)
register_solver(RK4StepSolver)
register_solver(EulerStepSolver)

functions = {}

diffusion_schemes_implicit = {"symmetric" :  "- dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                                    "- dw_dg_diffusion_flux.i.Omega(D.val, v, u)",
                     "non-symmetric": "- dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                                      "+ dw_dg_diffusion_flux.i.Omega(D.val, v, u)",
                     "incomplete": " dw_dg_diffusion_flux.i.Omega(D.val, u, v)"}

diffusion_schemes_explicit = {"symmetric" :  "- dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)" +
                                    "- dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])",
                     "non-symmetric": "- dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)" +
                                      "+ dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])",
                     "incomplete": " dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)"}



def local_register_function(fun):
    try:
        functions.update({fun.__name__: (fun,)})

    except AttributeError:  # Already a sfepy Function.
        fun = fun.function
        functions.update({fun.__name__: (fun,)})

    return fun


def get_cfl_setup(CFL=None, dt=None):

    if CFL is None and dt is None:
        raise ValueError("Specifiy either CFL or dt in CFL setup")

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
        if dt is None:
            _dt = dx / max_velo * CFL * order_corr
        else:
            _dt = dt
        # time_steps_N = int((tf - t0) / dt) * 2
        tn = int(nm.ceil((ts_conf.t1 - ts_conf.t0) / _dt))
        dtdx = _dt / dx

        ts_conf.dt = _dt
        ts_conf.n_step = tn
        output("Preprocessing hook setup_cfl_condition:")
        output("Approximation order of field {}({}) is {}".format(first_field_name, first_field.family_name, approx_order))
        output("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
        output("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, _dt))
        if dt is None:
            output("CFL coefficient was {0} and order correction 1/{1} = {2}".format(CFL,  (2 * approx_order + 1), order_corr))
        else:
            output("CFL coefficient {0} was ignored, dt specified directly".format(CFL))
        output("Courant number c = max(norm(a)) * dt/dx = {0}".format(max_velo * dtdx))
        output("------------------------------------------")
        output("Time stepping solver is {}".format(ts_conf.kind))

    return setup_cfl_condition


def build_transient_diffusion_advection_2D(mesh_hook, approx_order, CFL, t0, t1,
                                           dt=None,
                                           velo=(1, 0),
                                           diffusion_coef=0.02,
                                           regions={},
                                           diffusion_scheme_name="symmetric",
                                           ic_fun=lambda x, y: nm.zeros(shape=x.shape),
                                           source_fun=lambda x, y, t: nm.zeros(shape=x.shape),
                                           bc_funs={},
                                           ebcs=None,
                                           epbcs=None,
                                           **kwargs):

    functions = {}
    def local_register_function(fun):
        functions.update({fun.__name__: (fun,)})
        return fun

    for key, val in kwargs.items():
        if callable(val):
            local_register_function(val)

    @local_register_function
    def get_ic(x, ic=None):
        return get_ic(x[..., 0:1], x[..., 1:])

    @local_register_function
    def get_source(ts, coors, mode="qp", **kwargs):
        if mode == "qp":
            t = ts.dt * ts.step
            x_1 = coors[..., 0]
            x_2 = coors[..., 1]
            res = source_fun(x_1, x_2, t)
            return {"val": res[..., None, None]}

    @local_register_function
    def get_bc(ts, coors, bc, problem):
        return bc_funs(ts, coors, bc, problem)

    materials = {
        'a': ({'val': [velo], '.flux': 0.0},),
        'D': ({'val': [diffusion_coef], '.Cw': 1.},),
        'g': 'get_source'
    }

    filename_mesh = mesh_hook

    _regions = regions
    regions = {
        'Omega' : 'all',
    }
    regions.update(_regions)

    fields = {
        'density' : ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u' : ('unknown field', 'density', 0, 1),
        'v' : ('test field',    'density', 'u'),
    }

    ics = {
        'ic' : ('Omega', {'u.0' : 'get_ic'}),
    }

    dgebcs = ebcs if ebcs is not None else {}

    dgepbcs = epbcs if epbcs is not None else {}

    equations = {
        'Advection' : "  dw_volume_dot.i.Omega(v, u)" +

                      "  + dw_s_dot_mgrad_s.i.Omega(a.val, u[-1], v)" +
                      "  - dw_dg_advect_laxfrie_flux.i.Omega(a.val, v, u[-1])" +

                      "  - dw_laplace.i.Omega(D.val, v, u[-1])" +
                      " + " + diffusion_schemes[diffusion_scheme_name] +
                      " - " + str(diffusion_coef) + " * dw_dg_interior_penal.i.Omega(D.Cw, v, u[-1])"+

                      " + dw_volume_lvf.i.Omega(g.val, v)" +
                      " = 0 "
    }

    solvers = {
        "tss" : ('ts.tvd_runge_kutta_3',
                             {"t0": t0,
                              "t1": t1,
                              'limiter' : IdentityLimiter}),
        'nls' : ('nls.newton', {} ),
        'ls'  : ('ls.scipy_direct', {})
    }

    options = {
        'ts' : 'tss',
        'nls' : 'newton',
        'ls' : 'ls',
        'save_times' : 100,
        'pre_process_hook' : get_cfl_setup(CFL=CFL) if CFL is not None else get_cfl_setup(dt=dt)
    }
    return locals()


def get_1Dmesh_hook(XS, XE, n_nod):
    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            coors = nm.linspace(XS, XE, n_nod).reshape((n_nod, 1))
            conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
            mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
            descs = ['1_2']

            mesh = Mesh.from_data('laplace_1d', coors, None,
                                  [conn], [mat_ids], descs)
            return mesh

        elif mode == 'write':
            pass

    return mesh_hook


def get_gen_block_mesh_hook(dims, shape, centre, mat_id=0, name='block',
                   coors=None, verbose=True):
    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            mesh = gen_block_mesh(dims, shape, centre, mat_id=mat_id, name=name,
                   coors=coors, verbose=verbose)
            return mesh

        elif mode == 'write':
            pass

    return mesh_hook