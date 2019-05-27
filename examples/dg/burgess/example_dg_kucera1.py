from discrete.functions import Functionize
from examples.dg.example_dg_common import *
from toolz import reduce
from operator import mul
import numpy as nm

from sfepy.discrete.dg.dg_terms import NonlinScalarDotGradTerm, NonlinearHyperDGFluxTerm
from sfepy.discrete.dg.dg_terms import DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

register_term(NonlinScalarDotGradTerm)
register_term(NonlinearHyperDGFluxTerm)
register_term(DiffusionDGFluxTerm)
register_term(DiffusionInteriorPenaltyTerm)

example_name = "kucera1"
dim = 2  # int(example_name[example_name.index("D") - 1])

# filename_mesh = "mesh/messedquad2_diamond.vtk"
filename_mesh = "../mesh/square_tri2.mesh"

approx_order = 2
t0 = 0.
t1 = .5
CFL = .4

n_el_nod = int(reduce(mul, map(lambda i: approx_order + i + 1, range(dim))) /
               reduce(mul, range(1, dim + 1)))  # number of DOFs per element

regions = {
    'Omega': 'all',
    'left' : ('vertices in x == -1', 'edge'),
    'right': ('vertices in x == 1', 'edge'),
    'top' : ('vertices in y == 1', 'edge'),
    'bottom': ('vertices in y == -1', 'edge')
}

fields = {
    'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')  #
}

variables = {
    'u': ('unknown field', 'density', 0, 1),
    'v': ('test field', 'density', 'u'),
}

angle = - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
# velo = nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
velo = nm.array([[1., 1.]]).T


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


def get_ic(x, ic=None):
    return gsmooth(x[..., 0:1] - .4) * gsmooth(x[..., 1:] - .4)


def adv_fun(u):
    vu = velo.T * u[..., None]
    return vu


def adv_fun_d(u):
    v1 = velo.T * nm.ones(u.shape + (1,))
    return v1


burg_velo = velo.T / nm.linalg.norm(velo)


def burg_fun(u):
    vu = .5*burg_velo * u[..., None] ** 2
    return vu


def burg_fun_d(u):
    v1 = burg_velo * u[..., None]
    return v1


functions = {
    'get_ic'    : (get_ic,),
    'burg_fun'  : (burg_fun,),
    'burg_fun_d': (burg_fun_d,),
    'bc_funs' : (bc_funs,),
    'source_fun': (source_fun,)
}

diffusion_coef = 0.002
materials = {
    'a'     : ({'val': [velo], '.flux': 0.0},),
    'nonlin': ({'.fun': adv_fun, '.dfun': adv_fun_d},),
    'burg'  : ({'.fun': burg_fun, '.dfun': burg_fun_d},),
    'D'     : ({'val': [diffusion_coef], '.Cw': 1.},),
    # 'g'     : ({'function': source_fun},)
}


rhs = {
    'name' : 'g',
    'function' : 'source_fun',
}

ics = {
    'ic': ('Omega', {'u.0': 'get_ic'}),
}

dgebcs = {
    'u_left' : ('left', {'u.all': 'bc_funs', 'gradu.all': 'bc_funs'}),
    'u_right' : ('right', {'u.all': 'bc_funs', 'gradu.all': 'bc_funs'}),
    'u_bottom' : ('bottom', {'u.all': 'bc_funs', 'gradu.all': 'bc_funs'}),
    'u_top' : ('top', {'u.all': 'bc_funs', 'gradu.all': 'bc_funs'}),

}

integrals = {
    'i': 2 * approx_order,
}

equations = {
    'Advection': "dw_volume_dot.i.Omega(v, u)" +
                 # non-linear advection
                 " + dw_ns_dot_grad_s.i.Omega(burg.fun, burg.dfun, u[-1], v)" +
                 " - dw_dg_nonlinear_laxfrie_flux.i.Omega(a.flux, burg.fun, burg.dfun, v, u[-1])" +
                 #  diffusion
                 " - dw_laplace.i.Omega(D.val, v, u[-1]) + dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])"
                 " - "
                 + str(diffusion_coef) + "*"
                 + "dw_dg_interior_penal.i.Omega(D.Cw, v, u[-1])"
                 + " + dw_volume_lvf.i.Omega(g.val, v)"
                 " = 0"
}

solvers = {
    "tss": ('ts.euler',
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
