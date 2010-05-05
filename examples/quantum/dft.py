import os
cwd = os.path.split( os.path.join( os.getcwd(), __file__ ) )[0]

def common(mesh='../../tmp/mesh.vtk', dim=3, n_eigs=5, n_electron=5, tau=-1.0):
    assert dim in [2, 3]
    filename_mesh = mesh
    options = {
        'save_eig_vectors' : None,
        'n_eigs' : 10,
        'n_electron' : n_electron,
        'eigen_solver' : 'eigen1',
        'dft_solver' : 'broyden',
        'output_dir' : os.path.join( cwd, 'output/' ),
        'log_filename' : 'log.txt',
        'iter_fig_name' : 'iterations.pdf',
    }

    regions = {
        'Omega' : ('all', {}),
        'Surface' : ('nodes of surface', {'can_cells' : False}),
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
        'region' : 'Omega',

        'values' : {
            'val' : 0.5,
            'one' : 1.0,
        }
    }

    material_2 = {
        'name' : 'mat_v',
        'mode' : 'function',
        'region' : 'Omega',

        'function' : 'fun_v',
        'extra_args' : {'mode' : 'r^2'},
    }

    if dim == 3:
        base_approx = "3_4_P1"
    else:
        base_approx = "2_3_P1"

    field_0 = {
        'name' : 'field_Psi',
        'dim' : (1,1),
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

    ebcs = {
        'ZeroSurface' : ('Surface', {'Psi.0' : 0.0}),
        'VHSurface' : ('Surface', {'Psi.0' : 'core_pot'})
    }

    equations = {
        'lhs' : """  dw_laplace.i1.Omega( m.val, v, Psi )
                   + dw_mass_scalar_w.i1.Omega( mat_v.V, v, Psi )""",
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
        'i_max'      : 1000,
        'eps_a'      : 1e-12,
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 1,
        'eps_a'      : 1e-10,
        'eps_r'      : 1.0,
        'macheps'   : 1e-16,
        'lin_red'    : 1, # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp' : 0.001,
        'ls_on'      : 1.1,
        'ls_min'     : 1e-5,
        'check'     : 0,
        'delta'     : 1e-6,
        'is_plot'    : False,
        'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
    }

    solver_2 = {
        'name' : 'eigen1',
        'kind' : 'eig.pysparse',

        'tau' : -10.0,
        'eps_a' : 1e-5,
        'i_max' : 150,
        'method' : 'qmrs',
        'verbosity' : 0,
        'strategy' : 1,
    }

    fe = {
        'chunk_size' : 100000
    }
    return locals()
