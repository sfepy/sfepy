# 05.10.2007, c
# last revision: 25.02.2008
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/special/cube_sphere.mesh'

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

region_10 = {
    'name' : 'Bottom',
    'select' : 'vertices in (z < %f)' % -0.499,
    'kind' : 'facet',
}

region_11 = {
    'name' : 'Top',
    'select' : 'vertices in (z > %f)' % 0.499,
    'kind' : 'facet',
}

material_1 = {
    'name' : 'solid',
    
    'values' : {
        'lam' : 1e1,
        'mu' : 1e0,
        'density' : 1e-1,
    },
}

field_1 = {
    'name' : '3_displacement',
    'dtype' : 'real',
    'shape' : 'vector',
    'region' : 'Y',
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

ebc_1 = {
    'name' : 'Fix',
    'region' : 'Bottom',
    'dofs' : {'u.all' : 0.0},
}
ebc_2 = {
    'name' : 'Load',
    'region' : 'Top',
    'dofs' : {'u.[0,1]' : 0.2, 'u.2' : 0.5},
}

lcbc_1 = {
    'name' : 'rigid1',
    'region' : 'Y2',
    'dofs' : {'u.all' : 'rigid'},
}

integral_1 = {
    'name' : 'i',
    'order' : 2,
}

equations = {
    'balance' : """dw_lin_elastic_iso.i.Y( solid.lam, solid.mu, v, u ) = 0""",
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
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

from testsBasic import TestLCBC
output_name = 'test_lcbc_3d.vtk'

##
# 03.10.2007, c
class Test( TestLCBC ):
    pass
