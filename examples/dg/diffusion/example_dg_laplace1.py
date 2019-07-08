"""

Problem specs for

            div(grad u) = 0

on rectangle
                    u = 0
                    u_x = 0
    [0,b]┌---------------------------┐[a, b]
         |                           |
         |                           |
u_x = -a |         u(x,y)            | u_x = 0
u = 0    |                           | u = 0
         |                           |
    [0,0]└---------------------------┘[a, 0]
                    u_y = b
                    u = 0

solution to this is 1/2*x**2 - 1/2*y**2 - a*x + b*y
"""

from examples.dg.example_dg_common import *
from sfepy.discrete.dg.dg_terms import DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

example_name = "laplace1"
dim = 2

filename_mesh = get_gen_block_mesh_hook((1., 1.), (16, 16), (.5, .5))

approx_order = 1
diffusion_coef = 1
Cw = 10000
a = 1
b = 1
c = 0

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
    res = 1/2*x_1**2 - 1/2*x_2**2 - a*x_1 + b*x_2 + c
    return res


@local_register_function
def sol_fun(ts, coors, mode="qp", **kwargs):
    t = ts.time
    if mode == "qp":
        return {"u": analytic_sol(coors, t)[..., None, None]}

materials = {
    'D'     : ({'val': [diffusion_coef], '.Cw': Cw},),
}

dgebcs = {
    'u_left' : ('left', {'u.all': 1, 'grad.u.all': (0, 0)}),
    'u_right' : ('right', {'u.all': 0, 'grad.u.all': (0, 0)}),
    'u_bottom' : ('bottom', {'u.all': 0, 'grad.u.all': (0, 0)}),
    'u_top' : ('top', {'u.all': 0, 'grad.u.all': (0, 0)}),

}

integrals = {
    'i': 2 * approx_order,
}

diff_scheme_name = "symmetric"

equations = {
    'Advection': "dw_laplace.i.Omega(D.val, v, u) " +
                 diffusion_schemes_implicit[diff_scheme_name] +
                 " - " + str(diffusion_coef) + " * dw_dg_interior_penal.i.Omega(D.Cw, v, u)" +
                 " = 0"
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
