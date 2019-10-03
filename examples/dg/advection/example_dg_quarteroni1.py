"""
Based on

Antonietti, P., & Quarteroni, A. (2013). Numerical performance of discontinuous and stabilized
    continuous Galerkin methods for convection–diffusion problems.
Numerical Methods for Hyperbolic Equations, 75–85. https://doi.org/10.1201/b14172-9
"""

from examples.dg.example_dg_common import *


def define(filename_mesh=None, approx_order=1, flux=0, CFL=None, dt=None,
           Cw=100, diffusion_coef=1, diff_scheme_name="symmetric"):

    CFL=None
    dt=None

    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    example_name = "quart1"
    dim = 2

    if filename_mesh is None:
        filename_mesh = get_gen_block_mesh_hook((1., 1.), (20, 20), (.5, .5))

    velo = [1., 1.]

    angle = 0.0  # - nm.pi / 5
    rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                     [nm.sin(angle), nm.cos(angle)]])
    velo = nm.sum(rotm.T * nm.array(velo), axis=-1)[:, None]

    regions = {
        'Omega'     : 'all',
        'left' : ('vertices in x == 0', 'edge'),
        'right': ('vertices in x == 1', 'edge'),
        'top' : ('vertices in y == 1', 'edge'),
        'bottom': ('vertices in y == 0', 'edge')
    }

    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u': ('unknown field', 'f', 0),
        'v': ('test field', 'f', 'u'),
    }

    integrals = {
        'i': 2 * approx_order,
    }


    @local_register_function
    def bc_funs(ts, coors, bc, problem):
        # return 2*coors[..., 1]
        t = ts.dt*ts.step
        x_1 = coors[..., 0]
        x_2 = coors[..., 1]
        res = nm.zeros(nm.shape(x_1))

        sin = nm.sin
        cos = nm.cos
        exp = nm.exp
        pi = nm.pi

        if bc.diff == 0:
            if "left" in bc.name:
                res[:] = 0
            elif "right" in bc.name:
                res[:] = 0
            elif "bottom" in bc.name:
                res[:] = 0 #-2*sin(2*pi*x_1)
            elif "top" in bc.name:
                res[:] = 0

        elif bc.diff == 1:
            if "left" in bc.name:
                res = nm.stack((-2*pi*(x_2**2 - x_2),
                                    res),
                               axis=-2)
            elif "right" in bc.name:
                res = nm.stack((-2*pi*(x_2**2 - x_2), res,),
                               axis=-2)
            elif "bot" in bc.name:
                res = nm.stack((res,
                                sin(2*pi*x_1)),
                               axis=-2)
            elif "top" in bc.name:
                res = nm.stack((res,
                                -sin(2*pi*x_1)),
                               axis=-2)

        return res


    @local_register_function
    def source_fun(ts, coors, mode="qp", **kwargs):
        # t = ts.dt * ts.step
        eps = diffusion_coef
        sin = nm.sin
        cos = nm.cos
        exp = nm.exp
        sqrt = nm.sqrt
        pi = nm.pi
        if mode == "qp":
            x_1 = coors[..., 0]
            x_2 = coors[..., 1]
            res = -2*pi*(x_2**2 - x_2)*cos(2*pi*x_1)\
                  - 2*(2*pi**2*(x_2**2 - x_2)*sin(2*pi*x_1) - sin(2*pi*x_1))*eps\
                  - (2*x_2 - 1)*sin(2*pi*x_1)
            return {"val": res[..., None, None]}


    def analytic_sol(coors, t):
        x_1 = coors[..., 0]
        x_2 = coors[..., 1]
        sin = nm.sin
        pi = nm.pi
        res = -(x_2 ** 2 - x_2) * sin(2 * pi * x_1)
        return res


    @local_register_function
    def sol_fun(ts, coors, mode="qp", **kwargs):
        t = ts.time
        if mode == "qp":
            return {"u": analytic_sol(coors, t)[..., None, None]}


    dgebcs = {
        'u_left' : ('left', {'u.all': "bc_funs", 'grad.u.all' : "bc_funs"}),
        'u_top'  : ('top', {'u.all': "bc_funs", 'grad.u.all' :  "bc_funs"}),
        'u_bot'  : ('bottom', {'u.all': "bc_funs", 'grad.u.all' :  "bc_funs"}),
        'u_right': ('right', {'u.all': "bc_funs", 'grad.u.all' :  "bc_funs"}),
    }

    materials = {
        'a'     : ({'val': [velo], '.flux': flux},),
        'D'     : ({'val': [diffusion_coef], '.Cw': Cw},),
        'g'     : 'source_fun'
    }

    equations = {
        'balance': """
                    + dw_s_dot_mgrad_s.i.Omega(a.val, u, v)
                    - dw_dg_advect_laxfrie_flux.i.Omega(a.flux, a.val, v, u)
                   """
                   +
                   " - dw_laplace.i.Omega(D.val, v, u) " +
                   " + dw_dg_diffusion_flux.i.Omega(D.val, u, v)" +
                   " + dw_dg_diffusion_flux.i.Omega(D.val, v, u)" +
                   " - " + str(diffusion_coef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u)" +
                   " + dw_volume_lvf.i.Omega(g.val, v) = 0"

    }

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 5,
        'eps_a'      : 1e-8,
        'eps_r'      : 1.0,
        'macheps'   : 1e-16,
        'lin_red'    : 1e-2,  # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp' : 0.001,
        'ls_on'      : 0.99999,
        'ls_min'     : 1e-5,
        'check'     : 0,
        'delta'     : 1e-6,
    }

    options = {
        'nls'             : 'newton',
        'ls'              : 'ls',
        'output_format'   : 'msh',
        # 'pre_process_hook': get_cfl_setup(CFL)
    }
    return locals()


globals().update(define())
pass

