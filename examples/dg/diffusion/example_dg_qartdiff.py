# from __future__ import absolute_import
import numpy as nm
import sympy as sm

from examples.dg.example_dg_common import *
from sfepy.base.base import Struct
from sfepy.discrete import Function
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from examples.dg.diffusion.fem_conv_test import eval_expr, define_functions, SimpleExpression

example_name = "quartdiff1"
dim = 2

filename_mesh = get_gen_block_mesh_hook((1., 1.), (20, 20), (.5, .5))


approx_order = 1
diffusion_coef = 1
Cw = 10000
use_symbolic = True

materials = {
    'D': ({'val': [diffusion_coef], '.Cw': Cw},),
    'g': 'source_fun'
}

regions = {
    'Omega': 'all',
    'Gamma': ('vertices of surface', 'facet'),
}

fields = {
    'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
}

variables = {
    'u': ('unknown field', 'density', 0),
    'v': ('test field', 'density', 'u'),
}

if use_symbolic:
    estr = '-(y ** 2 - y) * sin(2 * pi * x)'
    expression = SimpleExpression(estr)
    expr = expression.define_expression()
    bc_fun, source_fun, sol_fun = expression.define_functions(expr)

    bc_fun = local_register_function(bc_fun)
    source_fun = local_register_function(source_fun)
    sol_fun = local_register_function(sol_fun)

else:
    @local_register_function
    def bc_fun(ts, coors, bc, problem):
        return analytic_sol(coors, ts.time)


    @local_register_function
    def source_fun(ts, coors, mode="qp", **kwargs):
        # t = ts.dt * ts.step
        eps = diffusion_coef
        sin = nm.sin
        cos = nm.cos
        pi = nm.pi
        if mode == "qp":
            x_1 = coors[..., 0]
            x_2 = coors[..., 1]
            res = -2 * pi * (x_2 ** 2 - x_2) * cos(2 * pi * x_1) \
                  - 2 * (2 * pi ** 2 * (x_2 ** 2 - x_2) * sin(2 * pi * x_1) - sin(2 * pi * x_1)) * eps \
                  - (2 * x_2 - 1) * sin(2 * pi * x_1)
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
    'fix': ('Gamma', {'u.all': 'bc_fun', 'grad.u.all' : "bc_fun"}),
}

integrals = {
    'i': 2 * approx_order,
}

diff_scheme_name = "symmetric"

equations = {
    'Temperature':  " - dw_laplace.i.Omega(D.val, v, u) " +
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
}