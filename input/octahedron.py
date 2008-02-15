# c: 15.02.2008, r: 15.02.2008
fileName_mesh = 'database/t.1.node'

material_2 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'coef' : 1.0,
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variables = {
    't' : ('field', 'unknown', 'temperature', (30,), 0),
    's' : ('field', 'test', 'temperature', (30,), 't'),
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

ebc = {
    'Gamma_Left' : ('T3', (30,), 2.0),
    'Gamma_Right' : ('T3', (30,), -2.0),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'Temperature' : """dw_laplace.i1.Omega( coef, s, t ) = 0"""
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
