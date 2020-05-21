"""
Based on

Kučera, V. (n.d.). Higher order methods for the solution of compressible flows.
Charles University. p. 21 eq. (1.39)
"""


from examples.dg.example_dg_common import *
from sfepy import data_dir



mesh_center = (0, 0)
mesh_size = (2, 2)

def define(filename_mesh=None,
           approx_order=2,

           adflux=0,
           limit=False,

           cw=10,
           diffcoef=0.002,
           diffscheme="symmetric",

           cfl=None,
           dt=1e-5,
           ):


    functions = {}
    def local_register_function(fun):
        try:
            functions.update({fun.__name__: (fun,)})

        except AttributeError:  # Already a sfepy Function.
            fun = fun.function
            functions.update({fun.__name__: (fun,)})

        return fun

    example_name = "kucera_dt_simp"
    dim = 2

    if filename_mesh is None:
        filename_mesh = data_dir + "/meshes/2d/square_tri2.mesh"

    t0 = 0.
    t1 = .1
    if dt is None and cfl is None:
        dt = 1e-5

    velo = [1., 1.]

    angle = 0  # - nm.pi / 5
    rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                     [nm.sin(angle), nm.cos(angle)]])
    velo = nm.sum(rotm.T * nm.array(velo), axis=-1)[:, None]
    burg_velo = velo.T / nm.linalg.norm(velo)

    regions = {
        'Omega': 'all',
        'left' : ('vertices in x == -1', 'edge'),
        'right': ('vertices in x == 1', 'edge'),
        'top'  : ('vertices in y == 1', 'edge'),
        'bottom': ('vertices in y == -1', 'edge')
    }

    fields = {
        'f': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u': ('unknown field', 'f', 0, 1),
        'v': ('test field', 'f', 'u'),
    }

    integrals = {
        'i': 5,
    }

    def analytic_sol(coors, t):
        x_1 = coors[..., 0]
        x_2 = coors[..., 1]
        sin = nm.sin
        pi = nm.pi
        exp = nm.exp
        res = -(exp(-t) - 1)*(sin(5*x_1*x_2) + sin(-4*x_1*x_2 + 4*x_1 + 4*x_2))
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
        sin = nm.sin
        cos = nm.cos
        exp = nm.exp
        if bc.diff == 0:
            if "left" in bc.name:
                res = -(exp(-t) - 1)*(sin(-5*x_2) + sin(8*x_2 - 4))
            elif "bottom" in bc.name:
                res = -(exp(-t) - 1) * (sin(-5 * x_1) + sin(8 * x_1 - 4))
            elif "right" in bc.name:
                res = -(exp(-t) - 1)*(sin(4) + sin(5*x_2))
            elif "top" in bc.name:
                res = -(exp(-t) - 1)*(sin(4) + sin(5*x_1))

        elif bc.diff == 1:
            if "left" in bc.name:
                res = nm.stack(((4*(x_2 - 1)*cos(4) - 5*x_2*cos(5*x_2))*(exp(-t) - 1),
                                    -5*(exp(-t) - 1)*cos(5*x_2)),
                               axis=-2)
            elif "bottom" in bc.name:
                res = nm.stack(((5*cos(-5*x_1) - 8*cos(8*x_1 - 4))*(exp(-t) - 1),
                                -(5*x_1*cos(-5*x_1) - 4*(x_1 - 1)*cos(8*x_1 - 4))*(exp(-t) - 1)),
                               axis=-2)

            elif "right" in bc.name:
                res = nm.stack(((4*(x_2 - 1)*cos(4) - 5*x_2*cos(5*x_2))*(exp(-t) - 1),
                                -5*(exp(-t) - 1)*cos(5*x_2)),
                               axis=-2)
            elif "top" in bc.name:
                res = nm.stack((-5*(exp(-t) - 1)*cos(5*x_1),
                                (4*(x_1 - 1)*cos(4) - 5*x_1*cos(5*x_1))*(exp(-t) - 1)),
                               axis=-2)

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
            res = (
                    + (5 * x_1 * cos(5 * x_1 * x_2) - 4 * (x_1 - 1) * cos(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1) ** 2 * (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2))
                    + (5 * x_2 * cos(5 * x_1 * x_2) - 4 * (x_2 - 1) * cos(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1) ** 2 * (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2))
                    - ((25 * x_1 ** 2 * sin(5 * x_1 * x_2) - 16 * (x_1 - 1) ** 2 * sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1)
                     + (25 * x_2 ** 2 * sin(5 * x_1 * x_2) - 16 * (x_2 - 1) ** 2 * sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1)) * diffcoef
                    + (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * exp(-t)
            )
            return {"val": res[..., None, None]}

    def adv_fun(u):
        vu = velo.T * u[..., None]
        return vu

    def adv_fun_d(u):
        v1 = velo.T * nm.ones(u.shape + (1,))
        return v1

    def burg_fun(u):
        vu = .5*burg_velo * u[..., None] ** 2
        return vu

    def burg_fun_d(u):
        v1 = burg_velo * u[..., None]
        return v1


    materials = {
        'a'     : ({'val': [velo], '.flux':adflux},),
        'D'     : ({'val': [diffcoef], '.Cw': cw},),
        'g'     : 'source_fun'
    }

    ics = {
        'ic': ('Omega', {'u.0': 0}),
    }

    dgebcs = {
        'u_left' : ('left', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
        'u_right' : ('right', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
        'u_bottom' : ('bottom', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
        'u_top' : ('top', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),

    }

    equations = {
      'balance':
         "dw_volume_dot.i.Omega(v, u)" +
         #  non-linear "advection"
         " - dw_ns_dot_grad_s.i.Omega(burg_fun, burg_fun_d, u[-1], v)" +
         " + dw_dg_nonlinear_laxfrie_flux.i.Omega(a.flux, burg_fun, burg_fun_d, v, u[-1])" +
         #  diffusion
         " + dw_laplace.i.Omega(D.val, v, u[-1])" +
         diffusion_schemes_explicit[diffscheme] +
         " + dw_dg_interior_penalty.i.Omega(D.val, D.Cw, v, u[-1])"
         # source
         + " - dw_volume_lvf.i.Omega(g.val, v)"
         " = 0"
    }

    solvers = {
        "tss.tvd_runge_kutta_3": ('ts.tvd_runge_kutta_3',
                                  {"t0": t0,
                                   "t1": t1,
                                   'limiters': {
                                       "f": MomentLimiter2D} if limit else {},
                                   'verbose': False}),
        "tss.euler": ('ts.euler',
                {"t0"     : t0,
                 "t1"     : t1,
                 'limiters': {"f": MomentLimiter2D} if limit else {},
                 'verbose': False}),
        'nls': ('nls.newton', {}),
        'ls' : ('ls.scipy_direct', {})
    }

    options = {
        'ts'              : 'tss.euler',
        'nls'             : 'nls.newton',
        'ls'              : 'ls.mumps',
        'save_times'      : 100,
        'output_format'   : 'msh',
		'file_format'     : 'gmsh-dg',
        'pre_process_hook': get_cfl_setup(CFL=cfl, dt=dt)
    }

    return locals()

globals().update(define())
