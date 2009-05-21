#filename_mesh = '../database/kostka_big.mesh'
filename_mesh = '../database/kostka_medium.mesh'

############# Laplace.

material_1 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'val' : 1.0
}

if filename_mesh == '../database/kostka_medium.mesh':
    region_1000 = {
        'name' : 'Omega',
        'select' : 'elements of group 0',
    }

    field_1 = {
        'name' : 'temperature',
        'dim' : (1,1),
        'domain' : 'Omega',
        'bases' : {'Omega' : '3_8_Q1'}
    }
    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o1_d3',
    }
    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
    }

elif filename_mesh == '../database/kostka_big.mesh':
    region_1000 = {
        'name' : 'Omega',
        'select' : 'elements of group 6',
    }

    field_1 = {
        'name' : 'temperature',
        'dim' : (1,1),
        'domain' : 'Omega',
        'bases' : {'Omega' : '3_4_P1'}
    }
    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'custom',
        'vals'    : [[1./3., 1./3., 1./3.]],
        'weights' : [0.5]
    }
    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_iterative',

        'method' : 'cg',
        'i_max'   : 1000,
        'eps_a'   : 1e-12,
    }

variable_1 = {
    'name' : 'T',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0, # order in the global vector of unknowns
}

variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
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
    'dofs' : {'T.0' : -3.0},
}
ebc_4 = {
    'name' : 'T1',
    'region' : 'Top',
    'dofs' : {'T.0' : 1.0},
}
ebc_3 = {
    'name' : 'T2',
    'region' : 'Bottom',
    'dofs' : {'T.0' : -1.0},
}
ebc_2 = {
    'name' : 'T3',
    'region' : 'Left',
    'dofs' : {'T.0' : 2.0},
}

equations = {
    'nice_equation' : """dw_laplace.i1.Omega( coef.val, s, T ) = 0""",
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

##
# FE assembling parameters.
fe = {
    'chunk_size' : 100000
}
