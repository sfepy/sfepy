from examples.dg.example_dg_common import *
from toolz import reduce
from operator import mul

from sfepy.discrete.dg.dg_terms import NonlinScalarDotGradTerm, NonlinearHyperDGFluxTerm

register_term(NonlinScalarDotGradTerm)
register_term(NonlinearHyperDGFluxTerm)

example_name = "burgess_2D_mquad"
dim = 2 #int(example_name[example_name.index("D") - 1])

filename_mesh = "mesh/messedquad2_diamond.vtk"

approx_order = 1
t0 = 0.
t1 = 1
CFL = .4

n_el_nod =  int(reduce(mul, map(lambda i: approx_order + i + 1, range(dim))) /
                         reduce(mul, range(1, dim + 1)))  # number of DOFs per element

# get_common(approx_order, CFL, t0, t1, None, get_ic)
angle = - nm.pi/5
rotm = nm.array([[nm.cos(angle),  -nm.sin(angle)],
                 [nm.sin(angle),  nm.cos(angle)]])
# velo = nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
velo = nm.array([[1., 1.]]).T



regions = {
    'Omega' : 'all',
    # 'Gamma_Left': ('vertices in (x < 0.055)', 'cell'),
}

fields = {
    'density' : ('real', 'scalar', 'Omega', str(approx_order)+'d', 'DG', 'legendre') #
}

variables = {
    'u' : ('unknown field', 'density', 0, 1),
    'v' : ('test field',    'density', 'u'),
}



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

def get_ic(x, ic=None):
    return gsmooth(x[..., 0:1]) * gsmooth(x[..., 1:])

def adv_fun(u):
    vu = velo.T * u[..., None]
    return vu

def adv_fun_d(u):
    v1 = velo.T * nm.ones(u.shape + (1,))
    return v1

burg_velo = velo.T / nm.linalg.norm(velo)

def burg_fun(u):
    vu = burg_velo * u[..., None]**2
    return vu

def burg_fun_d(u):
    v1 = 2 * burg_velo * u[..., None]
    return v1

functions = {
    'get_ic' : (get_ic,),
    'burg_fun': (burg_fun,),
    'burg_fun_d': (burg_fun_d,)
}

materials = {
    'a' : ({'val': [velo], '.flux': 0.0},),
    'nonlin' :({'.fun' : adv_fun, '.dfun' : adv_fun_d},),
    'burg': ({'.fun': burg_fun, '.dfun': burg_fun_d},),

}

ics = {
    'ic' : ('Omega', {'u.0' : icarray}),
}

integrals = {
    'i' : 2 * approx_order,
}

equations = {
    'Advection' : "dw_volume_dot.i.Omega(v, u)" +
                  "+ dw_ns_dot_grad_s.i.Omega(burg.fun, burg.dfun, u[-1], v)" +
                  "- dw_dg_nonlinear_laxfrie_flux.i.Omega(a.flux, burg.fun, burg.dfun, v, u[-1]) = 0"
}


solvers = {
    "tss" : ('ts.euler',
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
    'active_only' : False,
    'output_format' : 'msh',
    'pre_process_hook' : get_cfl_setup(CFL)
}