from examples.dg.example_dg_common import *

from sfepy.discrete.dg.dg_basis import get_n_el_nod

example_name = "adv_2D_mquad"
dim = 2  # int(example_name[example_name.index("D") - 1])

filename_mesh = "../mesh/messedquad2_diamond.vtk"

approx_order = 1
t0 = 0.
t1 = 1
CFL = .4

n_el_nod = get_n_el_nod(approx_order, dim, extended="quad" in filename_mesh)

# get_common(approx_order, CFL, t0, t1, None, get_ic)
angle = - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
# velo = nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
velo = nm.array([[1., 1.]]).T

materials = {
    'a': ({'val': [velo], '.flux': 0.0},),
}

regions = {
    'Omega': 'all',
    # 'Gamma_Left': ('vertices in (x < 0.055)', 'cell'),
}

fields = {
    'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')  #
}

variables = {
    'u': ('unknown field', 'density', 0, 1),
    'v': ('test field', 'density', 'u'),
}


def get_ic(x, ic=None):
    return gsmooth(x[..., 0:1]) * gsmooth(x[..., 1:])


# select IC based on mesh used
if "regularquad_diamond" in filename_mesh:
    icarray = nm.zeros((16 * n_el_nod,))
    icarray[2 * n_el_nod] = 1
elif "regulartri_diamond" in filename_mesh:
    icarray = nm.zeros((32 * n_el_nod,))
    icarray[15 * n_el_nod] = 1
    icarray[13 * n_el_nod] = 1
elif "messedtri2_diamond" in filename_mesh:
    icarray = nm.zeros((144 * n_el_nod,))
    icarray[19 * n_el_nod] = 1
    icarray[17 * n_el_nod] = 1
elif "messedquad2_diamond" in filename_mesh:
    icarray = nm.zeros((64 * n_el_nod,))
    icarray[44 * n_el_nod] = 1

functions = {
    'get_ic': (get_ic,)
}

ics = {
    'ic': ('Omega', {'u.0': icarray}),
}

integrals = {
    'i': 2 * approx_order,
}

equations = {
    'Advection': """
                   dw_volume_dot.i.Omega(v, u)
                   + dw_s_dot_mgrad_s.i.Omega(a.val, u[-1], v)
                   - dw_dg_advect_laxfrie_flux.i.Omega(a.val, v, u[-1]) = 0
                  """
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
