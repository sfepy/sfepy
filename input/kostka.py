# last revision: 18.02.2008
#fileName_mesh = 'database/kostka_big.mesh'
fileName_mesh = 'database/kostka_medium.mesh'

############# Laplace.

material_1 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'coef' : 1.0
}

if fileName_mesh == 'database/kostka_medium.mesh':
    region_1000 = {
        'name' : 'Omega',
        'select' : 'elements of group 0',
    }

    field_1 = {
        'name' : 'temperature',
        'dim' : (1,1),
        'flags' : (),
        'domain' : 'Omega',
        'bases' : {'Omega' : '3_8_Q1'}
    }
    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o1_d3',
    }

elif fileName_mesh == 'database/kostka_big.mesh':
    region_1000 = {
        'name' : 'Omega',
        'select' : 'elements of group 6',
    }

    field_1 = {
        'name' : 'temperature',
        'dim' : (1,1),
        'flags' : (),
        'domain' : 'Omega',
        'bases' : {'Omega' : '3_4_P1'}
    }
    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'custom',
        'vals'    : [[1./3., 1./3.]],
        'weights' : [0.5]
    }

variable_1 = {
    'name' : 'T',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'dofs' : (30,),
    'order' : 0, # order in the global vector of unknowns
}

variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dofs' : (30,),
    'dual' : 'T',
}

region_0 = {
    'name' : 'Surface',
    'select' : 'nodes of surface',
}
region_1 = {
    'name' : 'Bottom',
    'select' : 'nodes in (z < -0.4999999)',
}
region_2 = {
    'name' : 'Top',
    'select' : 'nodes in (z > 0.4999999)',
}
region_03 = {
    'name' : 'Left',
    'select' : 'nodes in (x < -0.4999999)',
}

ebc_1 = {
    'name' : 'T0',
    'region' : 'Surface',
    'dofs' : (30,),
    'value' : -3.0,
}
ebc_4 = {
    'name' : 'T1',
    'region' : 'Top',
    'dofs' : (30,),
    'value' : 1.0,
}
ebc_3 = {
    'name' : 'T2',
    'region' : 'Bottom',
    'dofs' : (30,),
    'value' : -1.0,
}
ebc_2 = {
    'name' : 'T3',
    'region' : 'Left',
    'dofs' : (30,),
    'value' : 2.0,
}

equations = {
    'nice_equation' : """dw_laplace.i1.Omega( coef, s, T ) = 0""",
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
