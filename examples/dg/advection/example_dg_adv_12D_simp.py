from examples.dg.example_dg_common import *

example_name = "adv_12D_simp"
dim = int(example_name[example_name.index("D") - 1])

filename_mesh = "..\mesh\mesh_simp_12D_01_20double.vtk"

approx_order = 1
t0 = 0.
t1 = .2
cfl = .4

# get_common(approx_order, cfl, t0, t1, None, get_ic)
angle = - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
velo = nm.array([[1., 0.]]).T
# velor = nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
materials = {
    'a': ({'val': [velo], '.flux': 0.0},),
}

regions = {
    'Omega'      : 'all',
    'Gamma_Left' : ('vertices in (x < 0.005)', 'facet'),
    'Gamma_Right': ('vertices in (x > 0.955)', 'facet'),

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


from sfepy.discrete.fem.periodic import match_y_line

functions = {
    'get_ic'      : (get_ic,),
    'match_y_line': (match_y_line,)
}

ics = {
    'ic': ('Omega', {'u.0': 'get_ic'}),
}
#
dgebcs = {
    'u_left' : ('Gamma_Left', {'u.all' : .5}),
    'u_right' : ('Gamma_Right', {'u.all' : -0.3}),
}

# dgepbc_1 = {
#     'name'  : 'u_rl',
#     'region': ['Gamma_Right', 'Gamma_Left'],
#     'dofs'  : {'u.all': 'u.all'},
#     'match' : 'match_y_line',
# }

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
    'pre_process_hook': get_cfl_setup(cfl)
}
