# last revision: 03.12.2007
fileName_mesh = 'database/simple.vtk'

field_1 = {
    'name' : '3_displacement',
    'dim' : (3,1),
    'flags' : (),
    'domain' : 'Y',
    'bases' : {'Y' : '3_4_P1'}
}

material_1 = {
    'name' : 'solid',
    'mode' : 'here',
    'region' : 'Y',
    'lam' : 1e1, # Cannot use lambda as it is Python keyword...
    'mu' : 1e0,
}

variables = {
    'u' : ('field', 'unknown', '3_displacement', (0, 1, 2), 0),
    'v' : ('field', 'test', '3_displacement', (0, 1, 2), 'u'),
}

# Whole domain $Y$.
region_1000 = {
    'name' : 'Y',
    'select' : 'all',
}

# EBC regions.
region_1 = {
    'name' : 'Left',
    'select' : 'nodes in (x < 0.001)'
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.099)'
}
region_100 = {
    'name' : 'AlsoRight',
    'select' : 'nodes in (x > 0.09)'
}
region_3 = {
    'name' : 'SomewhereTop',
    'select' : 'nodes in (z > 0.017) & (x > 0.01) & (x < 0.08)'
}

ebc = {
    'Left' : (('ZeroSurface', (0,1,2), 0.0 ),),
    'Right' : (('ZeroSurface', (0,1,2), 0.0 ),),
    'SomewhereTop' : (('PerturbedSurface', (2,), 0.01 ),),
    'AlsoRight' : (('ZeroVolume', (0,1,2), 0.0 ),),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'balance_of_forces' : """dw_sdcc.i1.Y( solid, v, u ) = 0""",
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
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

##
# FE assembling parameters.
fe = {
    'chunkSize' : 1000
}

##
# Functions.
from valec import *

##
# Make 'cinc' refer to a cinc_* function according to the mesh file name.
import os.path as op
trunk = op.splitext( op.basename( fileName_mesh ) )[0]
cinc = eval( 'cinc_' + trunk )
print cinc
del op, trunk
