# c: 21.09.2008
import os
import numpy as nm

from sfepy import data_dir
from sfepy.fem import MeshIO

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'
## filename_mesh = data_dir + '/meshes/2d/special/circle_in_square_small.mesh'
## filename_mesh = data_dir + '/meshes/3d/special/cube_sphere.mesh'
## filename_mesh = data_dir + '/meshes/2d/special/cube_cylinder.mesh'

omega = 1
omega_squared = omega**2

conf_dir = os.path.dirname(__file__)
io = MeshIO.any_from_filename(filename_mesh, prefix_dir=conf_dir)
bbox, dim = io.read_bounding_box( ret_dim = True )

geom = {3 : '3_4', 2 : '2_3'}[dim]

x_left, x_right = bbox[:,0]

regions = {
    'Y' : ('all', {}),
    'Y1' : ('elements of group 1', {}),
    'Y2' : ('elements of group 2', {}),
    'Y2_Surface': ('r.Y1 *n r.Y2', {'can_cells' : False}),
    'Left' : ('nodes in (x < %f)' % (x_left + 1e-3), {}),
    'Right' : ('nodes in (x > %f)' % (x_right - 1e-3), {}),
}

material_2 = {
    'name' : 'inclusion',

    # epoxy
    'function' : 'get_inclusion_pars',
}

def get_inclusion_pars(ts, coor, mode=None, region=None, ig=None):
    """TODO: implement proper 3D -> 2D transformation of constitutive
    matrices."""
    if mode == 'qp':
        n_nod, dim = coor.shape
        sym = (dim + 1) * dim / 2

        dielectric = nm.eye( dim, dtype = nm.float64 )
        # !!!
        coupling = nm.ones( (dim, sym), dtype = nm.float64 )
        #    coupling[0,1] = 0.2

        out = {
            # Lame coefficients in 1e+10 Pa.
            'lam' : 0.1798,
            'mu' : 0.148,
            # dielectric tensor
            'dielectric' : dielectric,
            # piezoelectric coupling
            'coupling' : coupling,
            'density' : 0.1142, # in 1e4 kg/m3
        }

        for key, val in out.iteritems():
            out[key] = nm.tile(val, (coor.shape[0], 1, 1))
        return out

functions = {
    'get_inclusion_pars' : (get_inclusion_pars,),
}

field_0 = {
    'name' : 'displacement',
    'dtype' : nm.float64,
    'shape' : dim,
    'region' : 'Y',
    'approx_order' : 1,
}

field_2 = {
    'name' : 'potential',
    'dtype' : nm.float64,
    'shape' : (1,),
    'region' : 'Y',
    'approx_order' : 1,
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
    'phi' : ('unknown field', 'potential', 1),
    'psi' : ('test field', 'potential', 'phi'),
}

ebcs = {
    'u1' : ('Left', {'u.all' : 0.0}),
    'u2' : ('Right', {'u.0' : 0.1}),
    'phi' : ('Y2_Surface', {'phi.all' : 0.0}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d%d' % dim,
}

equations = {
    '1' : """- %f * dw_mass_vector.i1.Y( inclusion.density, v, u )
             + dw_lin_elastic_iso.i1.Y( inclusion.lam, inclusion.mu, v, u )
             - dw_piezo_coupling.i1.Y2( inclusion.coupling, v, phi )
           = 0""" % omega_squared,
    '2' : """dw_diffusion.i1.Y( inclusion.dielectric, psi, phi )
           + dw_piezo_coupling.i1.Y2( inclusion.coupling, u, psi )
           = 0""",
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 100000
}

##
# Solvers etc.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'    : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp': 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'      : 0,
    'delta'      : 1e-6,
    'is_plot'    : False,
    'problem'    : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}
