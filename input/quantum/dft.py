def common(mesh='tmp/mesh.vtk', dim=3, n_eigs=5, tau=-1.0):
    assert dim in [2, 3]
    fileName_mesh = mesh
    options = {
        'saveEigVectors' : None,
        'squared' : False,
        'nEigs' : 10,
        'eigenSolver' : 'eigen1',
    }

    # Whole domain $Y$.
    region_1000 = {
        'name' : 'Omega',
        'select' : 'all',
    }

    # Domain $Y_2$.
    region_2 = {
        'name' : 'Surface',
        'select' : 'nodes of surface',
    }

    #def get_sphere( x, y, z, mode ):
    #    import numpy as nm
    #    r = x**2 + y**2 + z**2
    #    val = nm.where( (1 <= r) & (r <= 2), 1, 0 )
    #    return val

    #region_03 = {
    #    'name' : 'sphere',
    #    'select' : 'nodes by get_sphere( x, y, z, 0 )',
    #}

    material_1 = {
        'name' : 'm',
        'mode' : 'here',
        'region' : 'Omega',

        'val' : 0.5,
        'one' : 1.0,
    }

    material_2 = {
        'name' : 'matV',
        'mode' : 'function',
        'region' : 'Omega',

        'function' : 'funV',
        'extraArgs' : {'mode' : 'r^2'},
    }

    if dim == 3:
        base_approx = "3_4_P1"
    else:
        base_approx = "2_3_P1"

    field_0 = {
        'name' : 'field_Psi',
        'dim' : (1,1),
        'flags' : (),
        'domain' : 'Omega',
        'bases' : {'Omega' : base_approx}
    }

    if dim == 3:
        quadr = "gauss_o2_d3"
    else:
        quadr = "gauss_o2_d2"

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : quadr,
    }

    variable_1 = {
        'name' : 'Psi',
        'kind' : 'unknown field',
        'field' : 'field_Psi',
        'order' : 0,
    }
    variable_2 = {
        'name' : 'v',
        'kind' : 'test field',
        'field' : 'field_Psi',
        'dual' : 'Psi',
    }
    variable_3 = {
        'name' : 'V',
        'kind' : 'parameter field',
        'field' : 'field_Psi',
        'like' : 'Psi',
    }
    variable_4 = {
        'name' : 'n',
        'kind' : 'parameter field',
        'field' : 'field_Psi',
        'like' : 'Psi',
    }

    ebc_1 = {
        'name' : 'ZeroSurface',
        'region' : 'Surface',
        'dofs' : {'Psi.0' : 0.0},
    }

    equations = {
        'lhs' : """  dw_laplace.i1.Omega( m.val, v, Psi )
                   + dw_mass_scalar_variable.i1.Omega( matV.V, v, Psi )""",
        'rhs' : """dw_mass_scalar.i1.Omega( v, Psi )""",
        #'sphere' : """ d_volume_integrate.i1.sphere( n )""",
    }

    equations_vh = {
        'poisson' : """dw_laplace.i1.Omega( m.one, v, Psi )
                       = dw_mass_scalar.i1.Omega( v, n )"""
    }

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_iterative',

        'method' : 'cg',
        'iMax'      : 1000,
        'epsA'      : 1e-12,
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'iMax'      : 1,
        'epsA'      : 1e-10,
        'epsR'      : 1.0,
        'macheps'   : 1e-16,
        'linRed'    : 1, # Linear system error < (epsA * linRed).
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

    solver_2 = {
        'name' : 'eigen1',
        'kind' : 'eig.pysparse',

        'tau' : tau,
        'epsA' : 1e-5,
        'iMax' : 150,
        'method' : 'qmrs',
        'verbosity' : 0,
        'strategy' : 1,
    }

    fe = {
        'chunkSize' : 100000
    }
    return locals()
