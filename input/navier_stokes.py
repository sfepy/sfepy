# 26.02.2007, c
# last revision: 18.02.2008

fileName_mesh = 'database/pul_klikatak2.mesh'

field_1 = {
    'name' : '3_velocity',
    'dim' : (3,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1B'}
}

field_2 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

# Can use logical operations '&' (and), '|' (or).
region_1000 = {
    'name' : 'Omega',
    'select' : 'elements of group 6',
}

region_0 = {
    'name' : 'Walls',
    'select' : 'nodes of surface -n (r.Outlet +n r.Inlet)',
    'canCells' : False,
}
region_1 = {
    'name' : 'Inlet',
    'select' : 'nodes by cinc( x, y, z, 0 )', # In
    'canCells' : False,
}
region_2 = {
    'name' : 'Outlet',
    'select' : 'nodes by cinc( x, y, z, 1 )', # Out
    'canCells' : False,
}

ebc_1 = {
    'name' : 'Walls',
    'region' : 'Walls',
    'dofs' : (3, 4, 5),
    'value' : 0.0,
}
ebc_2 = {
    'name' : 'Inlet_y',
    'region' : 'Inlet',
    'dofs' : (4,),
    'value' : 1.0,
}
ebc_3 = {
    'name' : 'Inlet_xz',
    'region' : 'Inlet',
    'dofs' : (3, 5),
    'value' : 0.0,
}

material_1 = {
    'name' : 'fluid',
    'mode' : 'here',
    'region' : 'Omega',
    'viscosity' : 1.25e-3,
    'density' : 1e0,
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_velocity',
    'dofs' : (3, 4, 5),
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_velocity',
    'dofs' : (3, 4, 5),
    'dual' : 'u',
}
variable_3 = {
    'name' : 'p',
    'kind' : 'unknown field',
    'field' : 'pressure',
    'dofs' : (9,),
    'order' : 1,
}
variable_4 = {
    'name' : 'q',
    'kind' : 'test field',
    'field' : 'pressure',
    'dofs' : (9,),
    'dual' : 'p',
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}
integral_2 = {
    'name' : 'i2',
    'kind' : 'v',
    'quadrature' : 'gauss_o3_d3',
}

##
# Stationary Navier-Stokes equations.
equations = {
    'namespaces' : {
    'dw_' : ('div', 'grad', 'div_grad', 'convect'),
    },
    'balance' :
    """+ div_grad.i2.Omega( fluid, v, u ) + convect.i2.Omega( v, u )
       - grad.i1.Omega( v, p ) = 0""",
    'incompressibility' :
    """div.i1.Omega( q, u ) = 0""",
}

##
# Stokes equations.
## equations = {
##     'namespaces' : {
##     'dw_' : ('div', 'grad', 'divgrad'),
##     },
##     'balance' :
##     """+ div_grad.i1.Omega( fluid, v, u ) - grad.i1.Omega( v, p ) = 0""",
##     'incompressibility' :
##     """div.i1.Omega( q, u ) = 0""",
## }

##
# FE assembling parameters.
fe = {
    'chunkSize' : 1000
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'iMax'      : 5,
    'epsA'      : 1e-8,
    'epsR'      : 1.0,
    'macheps'   : 1e-16,
    'linRed'    : 1e-2, # Linear system error < (epsA * linRed).
    'lsRed'     : 0.1,
    'lsRedWarp' : 0.001,
    'lsOn'      : 0.99999,
    'lsMin'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'isPlot'    : False,
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

##
# Functions.
from valec import *

##
# Make 'cinc' refer to a cinc_* function according to the mesh file name.
import os.path as op
trunk = op.splitext( op.basename( fileName_mesh ) )[0]
print trunk
cinc = eval( 'cinc_' + trunk )
print cinc
del op, trunk
