"""

Simple example for second order ODE

    du^2/dx^2 = g

    u(0) = u(1) = 0
    du/dx(0) = du/dx(0) = 2*pi

"""
from examples.dg.example_dg_common import *


def define(filename_mesh=None, approx_order=1, Cw=100, diffusion_coef=1, CFL=0.4,
           use_symbolic=False, transient=False):
    t0 = 0
    t1 = 1

    example_name =  "trns_per_diff1D" if transient else "per_diff1D"
    dim = 1

    if filename_mesh is None:
        filename_mesh = get_1Dmesh_hook(0, 1, 2)


    materials = {
        'D': ({'val': [diffusion_coef], '.Cw': Cw},),
        'g': 'source_fun'
    }

    regions = {
        'Omega' : 'all',
        'Gamma' : ('vertices of surface', 'facet'),
        'left': ('vertices in x == 0', 'vertex'),
        'right': ('vertices in x == 1', 'vertex')
    }

    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u' : ('unknown field', 'f', 0),
        'v' : ('test field',    'f', 'u'),
    }

    # dgebcs = {
    #     'u_left': ('left', {'u.all': "bc_fun", 'grad.u.all': "bc_fun"}),
    #     'u_right': ('right', {'u.all': "bc_fun", 'grad.u.all': "bc_fun"}),
    # }

    dgepbc_1 = {
        'name'  : 'u_rl',
        'region': ['right', 'left'],
        'dofs': {'u.all': 'u.all'},
        'match': 'match_y_line',
    }

    integrals = {
        'i' : 2 * approx_order,
    }

    if transient:
        equations = {
            'Advection': " dw_volume_dot.i.Omega(v, u) " +

                         " - dw_laplace.i.Omega(D.val, v, u[-1]) " +
                         " + dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                         " + dw_dg_diffusion_flux.i.Omega(D.val, v, u)" +
                         " - " + str(diffusion_coef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u[-1])" +

                         " + dw_volume_lvf.i.Omega(g.val, v)"
                         " = 0"
        }

        solvers = {
            "tss": ('ts.tvd_runge_kutta_3',
                    {"t0": t0,
                     "t1": t1,
                     'limiter': IdentityLimiter,
                     'verbose': True}),
            'nls': ('nls.newton', {}),
            'ls': ('ls.scipy_direct', {})
        }

        options = {
            'ts': 'tss',
            'nls': 'newton',
            'ls': 'ls',
            'save_times': 100,
            'active_only': False,
            'output_format': 'vtk',
            'pre_process_hook': get_cfl_setup(CFL)
        }
    else:
        equations = {
            'Temperature': " - dw_laplace.i.Omega(D.val, v, u) " +
                           " + dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                           " + dw_dg_diffusion_flux.i.Omega(D.val, v, u)" +
                           " - " + str(diffusion_coef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u)" +
                           " + dw_volume_lvf.i.Omega(g.val, v) = 0"
        }
        solvers = {
            'ls': ('ls.auto_direct', {}),
            'newton': ('nls.newton',
                       {'i_max': 1,
                        'eps_a': 1e-10,
                        }),
        }

        options = {
            'nls': 'newton',
            'ls': 'ls',
            'output_format': 'vtk',

        }



    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    @local_register_function
    def bc_fun(ts, coors, bc, problem):
        x = coors[..., 0]
        res = nm.zeros(nm.shape(x))
        pi = nm.pi

        if bc.diff == 1:
            if "left" in bc.name:
                res[:] = 2*pi
            elif "right" in bc.name:
                res[:] = 2*pi
        else:
            # bc are zero
            pass
        return res

    @local_register_function
    def source_fun(ts, coors, mode="qp", **kwargs):
        # t = ts.dt * ts.step
        eps = diffusion_coef
        sin = nm.sin
        cos = nm.cos
        pi = nm.pi
        if mode == "qp":
            x = coors[..., 0]
            res = 4*pi**2*eps*sin(2*pi*x)
            return {"val": res[..., None, None]}

    def analytic_sol(coors, t=0):
        x = coors[..., 0]
        sin = nm.sin
        pi = nm.pi
        res = sin(2*pi*x)
        return res

    @local_register_function
    def sol_fun(ts, coors, mode="qp", **kwargs):
        t = ts.time
        if mode == "qp":
            return {"u": analytic_sol(coors, t)[..., None, None]}

    return locals()


globals().update(define())
