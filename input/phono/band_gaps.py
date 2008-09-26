# c: 25.09.2007, r: 10.09.2008
import os
import numpy as nm
from sfepy.fem.meshio import MeshIO

#fileName_mesh = 'database/phono/cube_sphere.mesh'
#fileName_mesh = 'database/phono/cube_cylinder.mesh'
filename_mesh = 'database/phono/mesh_circ21.mesh'
#fileName_mesh = 'database/phono/mesh_circ21_small.mesh'

cwd = os.path.split( os.path.join( os.getcwd(), __file__ ) )[0]

options = {
    'save_eig_vectors' : (10, 0),
    'eig_range' : (0, 30), # -> freq_range = eigs[slice(*eig_range)][[0, -1]]
    'freq_margins' : (10, 10), # % of freq_range

    'feps' : 1e-10, # frequency
    'zeps' : 1e-12, # zero finding
    'teps' : 1e-1, # eigenmomentum threshold
    'teps_rel' : True, # eigenmomentum threshold is relative w.r.t. largest one
    'freq_step' : 0.01, # % of freq_range
#    'eig_vector_transform' : ('select_in_plane', 'z', 1e-1),
#    'plot_tranform' : ('clip', (-20, 20)),
    'plot_tranform' : ('normalize', (-1, 1)),
    'squared' : False,

    'output_dir' : os.path.join( cwd, 'output/' ),

    'fig_name' : 'band_gaps.pdf',
    
#    'method' : 'eig.sgscipy', # 'eig.sgscipy' (default) or 'eig.symeig'

    'eigenmomentum' : {'var' : 'up',
                       'regions' : ['Y2'],
                       'term' : '%s * di_volume_integrate.i1.%s( %s )'},
    # Used to compute average density.
    'region_to_material' : {'Y1' : 'matrix',
                            'Y2' : 'inclusion',},
    'volume' : 'd_volume.i1.%s( uy )',
    'eig_problem' : 'simple',
}

regions = {
    'Y' : ('all', {}),
    'Y1' : ('elements of group 1', {}),
    'Y2' : ('elements of group 2', {}),
    'Y2_Surface': ('r.Y1 *n r.Y2', {'can_cells' : False}),
}

material_1 = {
    'name' : 'matrix',
    'mode' : 'here',
    'region' : 'Y1',

    # aluminium
    'lame' : {'lambda' : 5.898, 'mu' : 2.681}, # in 1e+10 Pa
    'density' : 0.2799, # in 1e4 kg/m3
}

material_2 = {
    'name' : 'inclusion',
    'mode' : 'here',
    'region' : 'Y2',

    # epoxy
    'lame' : {'lambda' : 0.1798, 'mu' : 0.148}, # in 1e+10 Pa
    'density' : 0.1142, # in 1e4 kg/m3
}

dim = MeshIO.any_from_filename( filename_mesh ).read_dimension()
geom = {3 : '3_4', 2 : '2_3'}[dim]

field_0 = {
    'name' : 'displacement_Y',
    'dim' : (dim,1),
    'domain' : 'Y',
    'bases' : {'Y' : '%s_P1' % geom}
}

field_1 = {
    'name' : 'displacement_Y2',
    'dim' : (dim,1),
    'domain' : 'Y2',
    'bases' : {'Y2' : '%s_P1' % geom}
}

variables = {
    'u' : ('unknown field', 'displacement_Y2', 0),
    'v' : ('test field', 'displacement_Y2', 'u'),
    'uy' : ('parameter field', 'displacement_Y', None),
    'up' : ('parameter field', 'displacement_Y2', 'u'),
}

ebc_1 = {
    'name' : 'ZeroSurface',
    'region' : 'Y2_Surface',
    'dofs' : {'u.all' : 0.0},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d%d' % dim,
}

equations = {
    'lhs' : """dw_lin_elastic_iso.i1.Y2( inclusion.lame, v, u )""",
    'rhs' : """dw_mass_vector.i1.Y2( inclusion.density, v, u )""",
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 100000
}

def clip( data, plot_range ):
    return nm.clip( data, *plot_range )

def normalize( data, plot_range ):
    aux = nm.arctan( data )
    return clip( aux, plot_range )

##
# 02.10.2007, c
def select_in_plane( vec, shape, normal_direction, eps ):
    n_nod, dim = shape
    dir_vecs = {2 : {'x': 0, 'y' : 1, 'z' : 1},
               3 : {'x': 0, 'y' : 1, 'z' : 2}}
    ident = nm.eye( dim, dtype = nm.float64 )
    dir_vec = ident[:,dir_vecs[dim][normal_direction]]

    proj = nm.dot( nm.reshape( vec, (n_nod, dim) ), dir_vec )
    if nm.any( nm.abs( proj ) > eps ):
        return nm.zeros_like( vec ), True
    else:
        return vec, False
