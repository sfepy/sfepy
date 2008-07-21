# 25.09.2007, c
# last revision: 08.04.2008
"""
u1 is a dummy variable used unly for volume computation.
"""

#fileName_mesh = 'database/phono/cube_sphere.mesh'
#fileName_mesh = 'database/phono/cube_cylinder.mesh'
file_name_mesh = 'database/phono/mesh_circ21.mesh'
#fileName_mesh = 'database/phono/mesh_circ21_small.mesh'

options = {
    'save_eig_vectors' : (10, 10),
    'eig_range' : (0, 30), # -> freq_range = eigs[slice(*eig_range)][[0, -1]]
    'freq_margins' : (10, 10), # % of freq_range
    'feps' : 1e-10, # frequency
    'zeps' : 1e-12, # zero finding
    'teps' : 1e-3, # eigenmomentum threshold
    'freq_step' : 0.01, # % of freq_range
#    'eig_vector_transform' : ('select_in_plane', 'z', 1e-1),
#    'plot_tranform' : ('clip', (-20, 20)),
    'plot_tranform' : ('normalize', (-1, 1)),
    'squared' : False,
    
#    'method' : 'eig.sgscipy', # 'eig.sgscipy' (default) or 'eig.symeig'
}

# Whole domain $Y$.
region_1000 = {
    'name' : 'Y',
    'select' : 'all',
}

# Domain $Y_1$.
region_1 = {
    'name' : 'Y1',
    'select' : 'elements of group 1',
}

# Domain $Y_2$.
region_2 = {
    'name' : 'Y2',
    'select' : 'elements of group 2',
}

# Surface of $Y_2$.
region_100 = {
    'name' : 'Y2_Surface',
    'select' : 'r.Y1 *n r.Y2',
    'can_cells' : False,
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

if file_name_mesh.find( 'cube_' ) >= 0:
    dim, geom = 3, '3_4'
else:
    dim, geom = 2, '2_3'

field_0 = {
    'name' : 'displacement_Y1',
    'dim' : (dim,1),
    'flags' : (),
    'domain' : 'Y1',
    'bases' : {'Y1' : '%s_P1' % geom}
}

field_1 = {
    'name' : 'displacement_Y2',
    'dim' : (dim,1),
    'flags' : (),
    'domain' : 'Y2',
    'bases' : {'Y2' : '%s_P1' % geom}
}

field_2 = {
    'name' : 'eigen_direction',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Y2',
    'bases' : {'Y2' : '%s_P1' % geom}
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : 'displacement_Y2',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : 'displacement_Y2',
    'dual' : 'u',
}
variable_3 = {
    'name' : 'u1',
    'kind' : 'parameter field',
    'field' : 'displacement_Y1',
    'like' : 'u',
}
variable_4 = {
    'name' : 'uc',
    'kind' : 'parameter field',
    'field' : 'eigen_direction',
    'like' : None,
}
variable_5 = {
    'name' : 'd',
    'kind' : 'parameter field',
    'field' : 'eigen_direction',
    'like' : None,
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
    'chunk_size' : 1000
}

import numpy as nm
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
