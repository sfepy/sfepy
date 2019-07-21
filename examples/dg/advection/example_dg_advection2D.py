from examples.dg.example_dg_common import *

mesh_center = (0.5, 0.25)
mesh_size = (1.0, 0.5)

def define(filename_mesh=None, approx_order=1, Cw=100, CFL=0.4, dt=None, angle=0,
           diffusion_coef=1, diff_scheme_name="symmetric", flux=0):

    example_name = "adv_2Dt0"
    dim = 2

    if filename_mesh is None:
        filename_mesh = "../mesh/mesh_tens_2D_01_20.vtk"

    t0 = 0.
    t1 = 1.

    # get_common(approx_order, CFL, t0, t1, None, get_ic)
    rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                     [nm.sin(angle), nm.cos(angle)]])
    velo = nm.sum(rotm.T * nm.array([0., 0.]), axis=-1)[:, None]
    materials = {
        'a': ({'val': [velo], '.flux': 0.0},),
    }

    regions = {
        'Omega'     : 'all',
        'left': ('vertices in x == 0', 'edge'),
        'right': ('vertices in x == 1', 'edge'),
        'top': ('vertices in y == 1', 'edge'),
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
        res = gsmooth(x_1) * gsmooth(x_2)
        return res

    @local_register_function
    def sol_fun(ts, coors, mode="qp", **kwargs):
        t = ts.time
        if mode == "qp":
            return {"u": analytic_sol(coors, t)[..., None, None]}

    def get_ic(x, ic=None):
        return gsmooth(x[..., 0:1]) * gsmooth(x[..., 1:])


    functions = {
        'get_ic': (get_ic,)
    }

    ics = {
        'ic': ('Omega', {'u.0': 'get_ic'}),
    }

    dgepbc_1 = {
        'name': 'u_rl',
        'region': ['right', 'left'],
        'dofs': {'u.all': 'u.all'},
        'match': 'match_y_line',
    }


    integrals = {
        'i': 2 * approx_order,
    }



    equations = {
        'Advection': """
                       dw_volume_dot.i.Omega(v, u)
                       + dw_s_dot_mgrad_s.i.Omega(a.val, u[-1], v)
                       - dw_dg_advect_laxfrie_flux.i.Omega(a.flux, a.val, v, u[-1]) = 0
                      """
    }

    solvers = {
        "tss": ('ts.tvd_runge_kutta_3',
                {"t0"     : t0,
                 "t1"     : t1,
                 'limiter': IdentityLimiter,
                 'verbose': True}),
        'nls': ('nls.newton', {}),
        'ls' : ('ls.scipy_direct', {})
    }

    options = {
        'ts'              : 'tss',
        'nls'             : 'newton',
        'ls'              : 'ls',
        'save_times'      : 100,
        'active_only'     : False,
        'output_format'   : 'msh',
        'pre_process_hook': get_cfl_setup(CFL) if dt is None else get_cfl_setup(dt=dt)
    }

    if dt is None:
        del dt

    return locals()


globals().update(define())
