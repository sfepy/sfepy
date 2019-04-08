from sfepy.discrete.dg.examples.example_dg_common import *

example_name = "adv_12D_simp"
dim = int(example_name[example_name.index("D") - 1])

filename_mesh = "mesh/simp_12D_mesh.vtk"

approx_order = 1
t0 = 0.
t1 = .2
CFL = .4


# get_common(approx_order, CFL, t0, t1, None, get_ic)
angle = - nm.pi/5
rotm = nm.array([[nm.cos(angle),  -nm.sin(angle)],
                 [nm.sin(angle),  nm.cos(angle)]])
velo = nm.array([[1., 0.]]).T
# velor = nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
materials = {
    'a' : ({'val': [velo], '.flux': 0.0},),
}

regions = {
    'Omega' : 'all',
}

fields = {
    'density' : ('real', 'scalar', 'Omega', str(approx_order)+'d', 'DG', 'legendre') #
}

variables = {
    'u' : ('unknown field', 'density', 0, 1),
    'v' : ('test field',    'density', 'u'),
}

def get_ic(x, ic=None):
    return gsmooth(x[..., 0:1])# * gsmooth(x[..., 1:])

functions = {
    'get_ic' : (get_ic,)
}

ics = {
    'ic' : ('Omega', {'u.0' : 'get_ic'}),
}

integrals = {
    'i' : 2 * approx_order,
}

equations = {
    'Advection' : """
                   dw_volume_dot.i.Omega(v, u)
                   + dw_s_dot_mgrad_s.i.Omega(a.val, u[-1], v)
                   - dw_dg_advect_laxfrie_flux.i.Omega(a.val, v, u[-1]) = 0
                  """
}

solvers = {
    "tss" : ('ts.tvd_runge_kutta_3',
                         {"t0": t0,
                          "t1": t1,
                          'limiter' : IdentityLimiter,
                          'verbose' : True}),
    'nls' : ('nls.newton',{} ),
    'ls'  : ('ls.scipy_direct', {})
}

options = {
    'ts' : 'tss',
    'nls' : 'newton',
    'ls' : 'ls',
    'save_times' : 100,
    'output_format' : 'msh',
    'pre_process_hook' : get_cfl_setup(CFL)
}