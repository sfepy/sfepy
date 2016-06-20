# 03.10.2007, c
from __future__ import absolute_import
from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_lame

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'

# Whole domain $Y$.
region_1000 = {
    'name' : 'Y',
    'select' : 'all',
}

# Domain $Y_1$.
region_1 = {
    'name' : 'Y1',
    'select' : 'cells of group 1',
}

# Domain $Y_2$.
region_2 = {
    'name' : 'Y2',
    'select' : 'cells of group 2',
}

# Domain $Y_3$.
region_3 = {
    'name' : 'Y3',
    'select' : 'vertices in (x > %f) & (x < %f) & (y > %f) & (y < %f)'\
    % (-0.3, 0.3, -0.48, -0.3),
}

wx = wy = 0.499
region_10 = {
    'name' : 'Bottom',
    'select' : 'vertices in (y < %f)' % -wy,
    'kind' : 'facet',
}

region_11 = {
    'name' : 'Top',
    'select' : 'vertices in (y > %f)' % wy,
    'kind' : 'facet',
}

material_1 = {
    'name' : 'solid',
    'values' : {
        'D' : stiffness_from_lame(2, 1e1, 1e0),
        'density' : 1e-1,
    },
}

field_1 = {
    'name' : '2_displacement',
    'dtype' : 'real',
    'shape' : 'vector',
    'region' : 'Y',
    'approx_order' : 2,
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '2_displacement',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '2_displacement',
    'dual' : 'u',
}

ebc_1 = {
    'name' : 'Fix',
    'region' : 'Bottom',
    'dofs' : {'u.all' : 0.0},
}
ebc_2 = {
    'name' : 'Load',
    'region' : 'Top',
    'dofs' : {'u.all' : 0.2},
}

lcbc_1 = {
    'name' : 'rigid1',
    'region' : 'Y2',
    'dofs' : {'u.all' : None},
    'kind' : 'rigid',
}
lcbc_2 = {
    'name' : 'rigid2',
    'region' : 'Y3',
    'dofs' : {'u.all' : None},
    'kind' : 'rigid',
}

integral_1 = {
    'name' : 'i',
    'order' : 2,
}

equations = {
    'balance' : """dw_lin_elastic.i.Y( solid.D, v, u ) = 0""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
}

from tests_basic import TestLCBC
output_name = 'test_lcbc_2d.vtk'

##
# 03.10.2007, c
class Test( TestLCBC ):
    pass
