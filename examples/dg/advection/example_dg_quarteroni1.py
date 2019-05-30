from examples.dg.example_dg_common import *

example_name = "adv_2D_tens"
dim = int(example_name[example_name.index("D") - 1])

filename_mesh = get_gen_block_mesh_hook((1., 1.), (3, 3), (.5, .5))

approx_order = 0
t0 = 0.
t1 = .2
CFL = .4

# get_common(approx_order, CFL, t0, t1, None, get_ic)
angle = 0.0  # - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
velo = nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
materials = {
    'a': ({'val': [velo], '.flux': 0.0},),
}

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
    'u': ('unknown field', 'density', 0),
    'v': ('test field', 'density', 'u'),
}

@local_register_function
def get_ic(x, ic=None):
    return gsmooth(x[..., 0:1]) * gsmooth(x[..., 1:])

ics = {
    'ic': ('Omega', {'u.0': 'get_ic'}),
}

integrals = {
    'i': 2 * approx_order,
}

dgebcs = {
    'u_left' : ('left', {'u.all': 1}),
    'u_top'  : ('top', {'u.all': 1}),
    'u_bot'  : ('bottom', {'u.all': 1}),
    'u_right': ('right', {'u.all': 1}),

}

equations = {
    'balance': """
                   dw_s_dot_mgrad_s.i.Omega(a.val, u, v) 
                   - dw_dg_advect_laxfrie_flux.i.Omega(a.flux, a.val, v, u) 
                   = 0
                   
                   
                  """
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
