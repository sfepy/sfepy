# 05.10.2007, c
# last revision: 10.12.2007
fileName_mesh = 'database/phono/cube_sphere.mesh'

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

region_10 = {
    'name' : 'Bottom',
    'select' : 'nodes in (z < %f)' % -0.499,
}

region_11 = {
    'name' : 'Top',
    'select' : 'nodes in (z > %f)' % 0.499,
}

material_1 = {
    'name' : 'solid',
    'mode' : 'here',
    'region' : 'Y',

    'lam' : 1e1,
    'mu' : 1e0,
    'density' : 1e-1,
}

field_1 = {
    'name' : '3_displacement',
    'dim' : (3,1),
    'flags' : (),
    'domain' : 'Y',
    'bases' : {'Y' : '3_4_P1'}
}

variables = {
    'u'  : ('field', 'unknown', '3_displacement', (0, 1, 2), 0),
    'v'  : ('field', 'test', '3_displacement', (0, 1, 2), 'u'),
}

ebc = {
    'Bottom' : (('Fix', (0, 1, 2), 0.0),),
    'Top' : (('Fix', (0, 1), 0.2), ('Fix', (2,), 0.5)),
}

lcbc = {
    'Y2' : ('RigidBody', (0, 1), 'rigid'),
#    'Y3' : ('RigidBody', (0, 1, 2), 'rigid'),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}

equations = {
    'balance' : """dw_sdcc.i1.Y( solid, v, u ) = 0""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'iMax'      : 1,
    'epsA'      : 1e-10,
    'epsR'      : 1.0,
    'macheps'   : 1e-16,
    'linRed'    : 1e-2, # Linear system error < (epsA * linRed).
    'lsRed'     : 0.1,
    'lsRedWarp' : 0.001,
    'lsOn'      : 1.1,
    'lsMin'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'isPlot'    : False,
    'linSolver' : 'umfpack',
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

##
# FE assembling parameters.
fe = {
    'chunkSize' : 1000
}

from testsBasic import TestLCBC
outputName = 'test_lcbc_3d.vtk'

##
# 03.10.2007, c
class Test( TestLCBC ):
    pass
