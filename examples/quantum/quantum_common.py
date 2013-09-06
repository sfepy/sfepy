def common(fun_v, mesh='../../tmp/mesh.vtk', n_eigs=5, tau=0.0):
    filename_mesh = mesh

    options = {
        'save_eig_vectors' : None,
        'n_eigs' : n_eigs,
        'eigen_solver' : 'eigen1',
    }

    # Whole domain $Y$.
    region_1000 = {
        'name' : 'Omega',
        'select' : 'all',
    }

    # Domain $Y_2$.
    region_2 = {
        'name' : 'Surface',
        'select' : 'vertices of surface',
        'kind' : 'facet',
    }

    functions = {
        'fun_v' : (fun_v,),
    }
    
    material_1 = {
        'name' : 'm',

        'values' : {
            'val' : 0.5,
        },
    }

    material_2 = {
        'name' : 'mat_v',

        'function' : 'fun_v',
    }

    field_0 = {
        'name' : 'field_Psi',
        'dtype' : 'real',
        'shape' : 'scalar',
        'region' : 'Omega',
        'approx_order' : 1,
    }

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'order' : 2,
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

    ebc_1 = {
        'name' : 'ZeroSurface',
        'region' : 'Surface',
        'dofs' : {'Psi.0' : 0.0},
    }

    equations = {
        'lhs' : """  dw_laplace.i1.Omega( m.val, v, Psi )
                   + dw_volume_dot.i1.Omega( mat_v.V, v, Psi )""",
        'rhs' : """dw_volume_dot.i1.Omega( v, Psi )""",
    }

    solver_2 = {
        'name' : 'eigen1',
        'kind' : 'eig.pysparse',

        'tau' : tau,
        'eps_a' : 1e-10,
        'i_max' : 150,
        'method' : 'qmrs',
        'verbosity' : 0,
        'strategy' : 1,
    }

    return locals()
