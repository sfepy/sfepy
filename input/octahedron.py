# c: 15.02.2008, r: 28.03.2008
fileName_mesh = 'database/octahedron.node'

material_2 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'K' : [[1.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 1.0]],
    'val' : 1.0,
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0,
}

variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 't',
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

def get_line( x, y, z, mode ):
    import numpy as nm
    if mode == 0:
        val = nm.where( (x + y) >= 3.5, 1, 0 )
    elif mode == 1:
        val = nm.where( (x + y) <= -3.5, 1, 0 )
    print mode, val
    return val

region_03 = {
    'name' : 'Gamma_Left',
    'select' : 'nodes by get_line( x, y, z, 0 )',
}
region_4 = {
    'name' : 'Gamma_Right',
    'select' : 'nodes by get_line( x, y, z, 1 )',
}

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Left',
    'dofs' : {'t.0' : 2.0},
}

ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Right',
    'dofs' : {'t.0' : -2.0},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'Temperature' : """dw_diffusion.i1.Omega( coef.K, s, t ) = 0"""
#    'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
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

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

fe = {
    'chunkSize' : 1000
}
