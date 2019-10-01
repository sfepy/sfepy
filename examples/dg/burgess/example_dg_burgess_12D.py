from examples.dg.example_dg_common import *
from examples.dg.burgess.example_dg_burgess_diff_2D import *

example_name = "burgess_12D"
dim = int(example_name[example_name.index("D") - 1])

filename_mesh = "../mesh/mesh_tens_12D_20.vtk"

approx_order = 1
t0 = 0.
t1 = 1
CFL = .4

regions = {
    'Omega'     : 'all',
    'Gamma_Left': ('vertices in (x < 0.055)', 'cell'),
}

fields = {
    'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')  #
}

variables = {
    'u': ('unknown field', 'density', 0, 1),
    'v': ('test field', 'density', 'u'),
}


def get_ic(x, ic=None):
    return gsmooth(x[..., 0:1])  # * gsmooth(x[..., 1:])


angle = - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
velo = nm.array([[1., 0.]]).T


def adv_fun(u):
    vu = velo.T * u[..., None]
    return vu


def adv_fun_d(u):
    v1 = velo.T * nm.ones(u.shape + (1,))
    return v1


burg_velo = velo.T / nm.linalg.norm(velo)


def burg_fun(u):
    vu = burg_velo * u[..., None] ** 2
    return vu


def burg_fun_d(u):
    v1 = 2 * burg_velo * u[..., None]
    return v1


functions = {
    'get_ic'    : (get_ic,),
    'burg_fun'  : (burg_fun,),
    'burg_fun_d': (burg_fun_d,)
}

materials = {
    'a'     : ({'val': [velo], '.flux': 0.0},),
}
ics = {
    'ic': ('Omega', {'u.0': 'get_ic'}),
}

# ebcs = {
#     'u_left' : ('Gamma_Left', {'u.all' : .5}),
#     # 'u_righ' : ('Gamma_Right', {'u.all' : -0.3}),
# }

integrals = {
    'i': 2 * approx_order,
}

equations = {
    'Advection': "dw_volume_dot.i.Omega(v, u)" +
                 "+ dw_ns_dot_grad_s.i.Omega(burg_fun, burg_fun_d, u[-1], v)" +
                 "- dw_dg_nonlinear_laxfrie_flux.i.Omega(burg_fun, burg_fun_d, v, u[-1]) = 0"
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
