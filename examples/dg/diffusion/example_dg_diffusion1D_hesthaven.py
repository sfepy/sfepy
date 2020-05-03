"""

Simple example for second order ODE

    du^2/dx^2 = g

    u(0) = u(1) = 0
    du/dx(0) = du/dx(0) = 2*pi

"""
from examples.dg.example_dg_common import *


def define(filename_mesh=None,
           approx_order=2,

          adflux=None,
           limit=False,

           cw=200,
           diffcoef=1,
           diffscheme="symmetric",

           cfl=0.001,
           dt=None,
           ):

    transient = True
    t0 = 0
    t1 = 0.1

    example_name = "hest_diff1"
    dim = 1

    mstart = 0
    mend = 2 * nm.pi

    if filename_mesh is None:
        filename_mesh = get_1Dmesh_hook(mstart, mend, 80)


    materials = {
        'D': ({'val': [diffcoef], '.Cw': cw},),
        'a': ({'val': 1.0, '.flux': 0.0},),
    }

    regions = {
        'Omega' : 'all',
        'Gamma' : ('vertices of surface', 'facet'),
        'left': ('vertices in x == 0', 'vertex'),
        'right': ('vertices in x > {}'.format(mend - 10e-5), 'vertex')
    }

    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u' : ('unknown field', 'f', 0, 1),
        'v' : ('test field',    'f', 'u'),
    }

    # dgebcs = {
    #     'u_left': ('left', {'u.all': 0, 'grad.u.all': 0}),
    #     'u_right': ('right', {'u.all': 0, 'grad.u.all': 0}),
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

    equations = {
        'diffusion': " dw_volume_dot.i.Omega(v, u) " +

                     " - dw_laplace.i.Omega(D.val, v, u[-1]) " +
                     " + dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)" +
                     " + dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])" +
                     " - dw_dg_interior_penalty.i.Omega(D.val, D.Cw, v, u[-1])" +

                     # " - dw_dg_diffusion_fluxHest1.i.Omega(D.val, v, u[-1]) " +

                     " = 0"
    }

    solvers = {
        "tss": ('ts.runge_kutta_4',
                {"t0": t0,
                 "t1": t1,
                 'limiter': IdentityLimiter,
                 'verbose': False}),
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
        'pre_process_hook': get_cfl_setup(cfl)
    }


    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    def analytic_sol(coors, t=0):
        x = coors[..., 0]
        try:
            res = nm.exp(-t) * nm.sin(x)
        except ValueError:
            res = nm.exp(-t)[None, :] * nm.sin(x)[:, None]

        return res

    @local_register_function
    def sol_fun(ts, coors, mode="qp", **kwargs):
        t = ts.time
        if mode == "qp":
            return {"u": analytic_sol(coors, t)[..., None, None]}

    @local_register_function
    def get_ic(x, ic=None):
        return nm.sin(x)

    ics = {
        'ic': ('Omega', {'u.0': 'get_ic'}),
    }

    return locals()


globals().update(define())
