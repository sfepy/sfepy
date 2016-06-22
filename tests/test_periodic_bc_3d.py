# 04.06.2007, c
# last revision: 25.02.2008
from __future__ import absolute_import
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/various_formats/small3d.mesh'

material_1 = {
    'name' : 'coef',
    'values' : {'coef' : 1.0},
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}
region_1 = {
    'name' : 'Left',
    'select' : 'vertices in (x < -0.499)',
    'kind' : 'facet',
}
region_2 = {
    'name' : 'Right',
    'select' : 'vertices in (x > 0.499)',
    'kind' : 'facet',
}
region_3 = {
    'name' : 'Near',
    'select' : 'vertices in (y < -0.499)',
    'kind' : 'facet',
}
region_4 = {
    'name' : 'Far',
    'select' : 'vertices in (y > 0.499)',
    'kind' : 'facet',
}
region_5 = {
    'name' : 'Bottom',
    'select' : 'vertices in (z < -0.499)',
    'kind' : 'facet',
}
region_6 = {
    'name' : 'Top',
    'select' : 'vertices in (z > 0.499)',
    'kind' : 'facet',
}

field_1 = {
    'name' : '3_displacement',
    'dtype' : 'real',
    'shape' : (3,),
    'region' : 'Omega',
    'approx_order' : 2,
}

field_2 = {
    'name' : 'pressure',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_displacement',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_displacement',
    'dual' : 'u',
}
variable_3 = {
    'name' : 'p',
    'kind' : 'unknown field',
    'field' : 'pressure',
    'order' : 1,
}
variable_4 = {
    'name' : 'q',
    'kind' : 'test field',
    'field' : 'pressure',
    'dual' : 'p',
}

ebcs = {}
epbc_10 = {
    'name' : 'rl',
    'region' : ['Left', 'Right'],
    'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
    'match' : 'match_x_plane',
}
epbc_12 = {
    'name' : 'tb',
    'region' : ['Top', 'Bottom'],
    'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
    'match' : 'match_z_plane',
}
epbc_13 = {
    'name' : 'nf',
    'region' : ['Near', 'Far'],
    'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
    'match' : 'match_y_plane',
}

from sfepy.discrete.fem.periodic import match_x_plane, match_y_plane, match_z_plane

functions = {
    'match_x_plane' : (match_x_plane,),
    'match_y_plane' : (match_y_plane,),
    'match_z_plane' : (match_z_plane,),
}

from test_periodic_bc_2d import Test
