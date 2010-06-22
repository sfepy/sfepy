# 04.06.2007, c
# last revision: 25.02.2008
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
    'select' : 'nodes in (x < -0.499)',
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.499)',
}
region_3 = {
    'name' : 'Near',
    'select' : 'nodes in (y < -0.499)',
}
region_4 = {
    'name' : 'Far',
    'select' : 'nodes in (y > 0.499)',
}
region_5 = {
    'name' : 'Bottom',
    'select' : 'nodes in (z < -0.499)'
}
region_6 = {
    'name' : 'Top',
    'select' : 'nodes in (z > 0.499)'
}

field_1 = {
    'name' : '3_displacement',
    'dim' : (3,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P2'}
}

field_2 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
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

fe = {
    'chunk_size' : 1000
}

from sfepy.fem.periodic import *

functions = {
    'match_x_plane' : (match_x_plane,),
    'match_y_plane' : (match_y_plane,),
    'match_z_plane' : (match_z_plane,),
}

from test_periodic_bc_2d import Test
