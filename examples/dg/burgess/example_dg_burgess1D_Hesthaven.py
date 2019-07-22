"""
Based on

Jan S Hesthaven. Discontinuous Galerkin methods Lecture 8.
Brown University, Jan.Hesthaven@Brown.edu


"""

from examples.dg.example_dg_common import *


def define(filename_mesh=None, approx_order=1, Cw=1,
           diffusion_coef=0.002, flux=0.0, CFL = .1, diff_scheme_name="symmetric"):


    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    example_name = "burgess_hesthaven"
    dim = 1

    if filename_mesh is None:
        filename_mesh = get_1Dmesh_hook(-1, 1, 20)

    t0 = 0.
    t1 = 1.
    velo = 1.

    regions = {
        'Omega': 'all',
        'left' : ('vertices in x == -1', 'vertex'),
        'right': ('vertices in x == 1', 'vertex'),
    }

    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')  #
    }

    variables = {
        'u': ('unknown field', 'f', 0, 1),
        'v': ('test field', 'f', 'u'),
    }

    integrals = {
        'i': 2 * approx_order,
    }

    def analytic_sol(coors, t):
        x = coors[..., 0]
        eps = diffusion_coef
        tanh = nm.tanh
        res = -tanh(-1/4*(2*t - 2*x - 1)/eps) + 1
        return res

    @local_register_function
    def sol_fun(ts, coors, mode="qp", **kwargs):
        t = ts.time
        if mode == "qp":
            return {"u": analytic_sol(coors, t)[..., None, None]}

    @local_register_function
    def bc_fun(ts, coors, bc, problem):
        x = coors[..., 0]
        res = nm.zeros(nm.shape(x))
        pi = nm.pi
        t = ts.time

        eps = diffusion_coef
        tanh = nm.tanh

        if bc.diff == 1:
            if "left" in bc.name:
                res[:] = 1/2*(tanh(-1/4*(2*t + 1)/eps)**2 - 1)/eps
            elif "right" in bc.name:
                res[:] = 1/2*(tanh(-1/4*(2*t - 3)/eps)**2 - 1)/eps
        else:
            if "left" in bc.name:
                res[:] = -tanh(-1/4*(2*t + 1)/eps) + 1
            elif "right" in bc.name:
                res[:] = -tanh(-1/4*(2*t - 3)/eps) + 1

        return res

    @local_register_function
    def get_ic(x, ic=None):
        tanh = nm.tanh
        eps = diffusion_coef
        return -tanh(1/4*(2*x + 1)/eps) + 1

    @local_register_function
    def adv_fun(u):
        vu = velo * u[..., None]
        return vu

    @local_register_function
    def adv_fun_d(u):
        v1 = velo * nm.ones(u.shape + (1,))
        return v1

    @local_register_function
    def burg_fun(u):
        vu = 1/2 * u[..., None] ** 2
        return vu

    @local_register_function
    def burg_fun_d(u):
        v1 = u[..., None]
        return v1


    materials = {
        'a'     : ({'val': [velo], '.flux': flux},),
        'nonlin': ({'.fun': adv_fun, '.dfun': adv_fun_d},),
        'burg'  : ({'.fun': burg_fun, '.dfun': burg_fun_d},),
        'D'     : ({'val': [diffusion_coef], '.Cw': Cw},),
    }

    ics = {
        'ic': ('Omega', {'u.0': 'get_ic'}),
    }

    dgebcs = {
        'u_left' : ('left', {'u.all': 'bc_fun', 'grad.u.all': 'bc_fun'}),
        'u_right' : ('right', {'u.all': 'bc_fun', 'grad.u.all': 'bc_fun'}),

    }

    equations = {
                     # temporal der
        'balance':   "dw_volume_dot.i.Omega(v, u)" +
                     #  non-linear "advection"
                     " + dw_ns_dot_grad_s.i.Omega(burg.fun, burg.dfun, u[-1], v)" +
                     " - dw_dg_nonlinear_laxfrie_flux.i.Omega(a.flux, burg.fun, burg.dfun, v, u[-1])" +
                     #  diffusion
                     " - dw_laplace.i.Omega(D.val, v, u[-1])"
                     " + dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)" +
                     " + dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])" +
                     " - " + str(diffusion_coef) + "*dw_dg_interior_penal.i.Omega(D.Cw, v, u[-1])"
                     " = 0"
    }

    solvers = {
        "tss": ('ts.euler',
                {"t0"     : t0,
                 "t1"     : t1,
                 # 'limiter': IdentityLimiter,
                 'verbose': True}),
        'nls': ('nls.newton', {}),
        'ls' : ('ls.scipy_direct', {})
    }

    options = {
        'ts'              : 'tss',
        'nls'             : 'newton',
        'ls'              : 'ls',
        'save_times'      : 100,
        'output_format'   : 'msh',
        'pre_process_hook': get_cfl_setup(CFL)
    }

    return locals()


globals().update(define())
