from examples.dg.example_dg_common import *
from sfepy.discrete.dg.dg_terms import DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

example_name = "tdiffusion_only1"
dim = int(example_name[example_name.index("D") - 1])

filename_mesh = "../mesh/mehs_tens_2D_01_20.vtk"

approx_order = 2
t0 = 0.
t1 = .1
CFL = .4
diffusion_coef = 2*1e-2
Cw = 100

regions = {
    'Omega'     : 'all',
    'left' : ('vertices in x == 0', 'edge'),
    'right': ('vertices in x == 1', 'edge'),
    'top' : ('vertices in y == 1', 'edge'),
    'bottom': ('vertices in y == 0', 'edge')
}
fields = {
    'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')  #
}

variables = {
    'u': ('unknown field', 'density', 0, 1),
    'v': ('test field', 'density', 'u'),
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
                 + (25 * x_2 ** 2 * sin(5 * x_1 * x_2) - 16 * (x_2 - 1) ** 2 * sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1)) * diffusion_coef
                + (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * exp(-t)
        )
        return {"val": res[..., None, None]}

materials = {
    'D'     : ({'val': [diffusion_coef], '.Cw': Cw},),
    'g'     : 'source_fun'
}

ics = {
    'ic': ('Omega', {'u.0': 'get_ic'}),
}

dgebcs = {
    'u_left' : ('left', {'u.all': 'bc_funs', 'grad.u.all': (-a, 0)}),
    'u_right' : ('right', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
    'u_bottom' : ('bottom', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),
    'u_top' : ('top', {'u.all': 'bc_funs', 'grad.u.all': 'bc_funs'}),

}

integrals = {
    'i': 2 * approx_order,
}

diff_scheme_name = "symmetric"

equations = {
    'Advection': " dw_volume_dot.i.Omega(v, u) " +

                 " - dw_laplace.i.Omega(D.val, v, u[-1]) " +
                 " + " + diffusion_schemes_explicit[diff_scheme_name] +
                 " - " + str(diffusion_coef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u[-1])" +

                 " + dw_volume_lvf.i.Omega(g.val, v)"
                 " = 0"
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
    'pre_process_hook': get_cfl_setup(CFL)
}
