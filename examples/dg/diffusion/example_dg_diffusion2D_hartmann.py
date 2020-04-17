"""
Based on

Ralf Hartmann. Numerical Analysis of Higher Order Discontinuous Galerkin Finite Element methods.
Institute of Aerodynamic and Flow Technology DLR (German Aerospace Center). 13. Oct. 2008
"""
from examples.dg.example_dg_common import *


def define(filename_mesh=None,
           approx_order=2,

          adflux=None,
           limit=False,

           cw=100,
           diffcoef=1,
           diffscheme="symmetric",

           cfl=0.4,
           dt=None,
           ):

    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    example_name = "diffusion_only_hartmann"
    dim = 2

    if filename_mesh is None:
        filename_mesh = "mesh/mesh_tens_2D_01_20.vtk"

    regions = {
        'Omega'     : 'all',
        'left' : ('vertices in x == 0', 'edge'),
        'right': ('vertices in x == 1', 'edge'),
        'top' : ('vertices in y == 1', 'edge'),
        'bottom': ('vertices in y == 0', 'edge')
    }
    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')  #
    }

    variables = {
        'u': ('unknown field', 'f', 0, 1),
        'v': ('test field', 'f', 'u'),
    }

    def analytic_sol(coors, t):
        x_1 = coors[..., 0]
        x_2 = coors[..., 1]
        sin = nm.sin
        pi = nm.pi
        exp = nm.exp
        res = sin(1/2*pi*x_1)*sin(1/2*pi*x_2)
        return res


    @local_register_function
    def sol_fun(ts, coors, mode="qp", **kwargs):
        t = ts.time
        if mode == "qp":
            return {"u": analytic_sol(coors, t)[..., None, None]}

    @local_register_function
    def bc_funs(ts, coors, bc, problem):
        # return 2*coors[..., 1]
        t = ts.dt*ts.step
        x_1 = coors[..., 0]
        x_2 = coors[..., 1]
        res = nm.zeros(x_1.shape)
        sin = nm.sin
        cos = nm.cos
        exp = nm.exp
        pi = nm.pi
        if bc.diff == 0:
            if "left" in bc.name:
                res[:] = 0
            elif "bottom" in bc.name:
                res[:] = 0
            elif "right" in bc.name:
                res[:] = sin(1/2*pi*x_2)
            elif "top" in bc.name:
                res[:] = sin(1/2*pi*x_1)

        elif bc.diff == 1:
            if "left" in bc.name:
                res = nm.stack((1/2*pi*sin(1/2*pi*x_2), res), axis=-2)
            elif "bottom" in bc.name:
                res = nm.stack((res, 1/2*pi*sin(1/2*pi*x_1)), axis=-2)
            elif "right" in bc.name:
                res = nm.stack((res, 1/2*pi*cos(1/2*pi*x_2)), axis=-2)
            elif "top" in bc.name:
                res = nm.stack((1/2*pi*cos(1/2*pi*x_1), res), axis=-2)

        return res

    @local_register_function
    def source_fun(ts, coors, mode="qp", **kwargs):
        if mode == "qp":
            t = ts.dt * ts.step
            x_1 = coors[..., 0]
            x_2 = coors[..., 1]
            sin = nm.sin
            cos = nm.cos
            exp = nm.exp
            pi = nm.pi
            eps = diffcoef
            res = 1/2*pi**2*eps*sin(1/2*pi*x_1)*sin(1/2*pi*x_2)
            return {"val": res[..., None, None]}

    materials = {
        'D'     : ({'val': [diffcoef], '.Cw': cw},),
        'g'     : 'source_fun'
    }

    ics = {
        'ic': ('Omega', {'u.0': 0}),
    }

    dgebcs = {
        'u_left' : ('left', {'u.all': 'bc_funs', 'grad.u.all': "bc_funs"}),
        'u_right' : ('right', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
        'u_bottom' : ('bottom', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
        'u_top' : ('top', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),

    }

    integrals = {
        'i': 2 * approx_order,
    }

    equations = {
        'Advection':   " - dw_laplace.i.Omega(D.val, v, u) " +
                       " + dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                       " + dw_dg_diffusion_flux.i.Omega(D.val, v, u)" +
                       " - " + str(diffcoef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u)" +
                       " + dw_volume_lvf.i.Omega(g.val, v) = 0"
    }

    solver_0 = {
        'name': 'ls',
        'kind': 'ls.scipy_direct',
        # 'kind': 'ls.mumps',

    }

    solver_1 = {
        'name': 'newton',
        'kind': 'nls.newton',

        'i_max': 5,
        'eps_a': 1e-8,
        'eps_r': 1.0,
        'macheps': 1e-16,
        'lin_red': 1e-2,  # Linear system error < (eps_a * lin_red).
        'ls_red': 0.1,
        'ls_red_warp': 0.001,
        'ls_on': 0.99999,
        'ls_min': 1e-5,
        'check': 0,
        'delta': 1e-6,
    }

    options = {
        'nls': 'newton',
        'ls': 'ls',
        'output_format': 'msh',
        # 'pre_process_hook': get_cfl_setup(cfl)
    }
    return locals()


globals().update(define())
pass
