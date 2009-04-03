import numpy as nm

def define():
    filename_mesh = 'database/simple.mesh'

    options = {
        'nls' : 'newton',
        'ls' : 'ls',
    }

    material_2 = {
        'name' : 'ac',
        'mode' : 'here',
        'region' : 'Omega',
        'one' : 1.0,
        'v_n' : 1.0, # m/s
        'rhs' : -1000.0,
        'k' : 1.0,
        'omega' : 1000.0,
        'c' : 343.0, # m/s
        'rho' : 1.55, # kg/m^3
    }

    field_1 = {
        'name' : 'accoustic_pressure',
        'dim' : (1,1),
#        'dtype' : nm.float64,
        'dtype' : nm.complex128,
        'domain' : 'Omega',
        'bases' : {'Omega' : '3_4_P1'}
    }

    variables = {
        'p'   : ('unknown field', 'accoustic_pressure', 0, 'previous'),
        'q'   : ('test field',    'accoustic_pressure', 'p'),
    }

    regions = {
        'Omega' : ('all', {}),
        'Gamma_Left' : ('nodes in (x < 0.00001)', {}),
        'Gamma_Right' : ('nodes in (x > 0.099999)', {}),
    }

    ebcs = {
#        'p_in' : ('Gamma_Left', {'p.0' : 1.0 + 2j}),
        'p_out' : ('Gamma_Right', {'p.0' : -2.0 + 1j}),
    }

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o1_d3',
    }

    integral_3 = {
        'name' : 'isurf',
        'kind' : 's',
        'quadrature' : 'gauss_o2_d2',
    }

    ac = material_2
    equations = {
        'Temperature' :
        """%s * dw_mass_scalar.i1.Omega( q, p )
              - dw_laplace.i1.Omega( ac.one, q, p )
         = %s * dw_surface_integrate.isurf.Gamma_Left( q ) """ \
        % (ac['k'],
           1j * ac['omega'] * ac['rho'] * ac['v_n'])
    }
 #         1j * ac['omega'] / ac['c'],
 #         - %s * dw_surface_dot.isurf.Gamma_Right( q, p )

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.umfpack',
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 1,
        'eps_a'      : 1e-10,
        'eps_r'      : 1.0,
        'macheps'   : 1e-16,
        'lin_red'    : 1e-1, # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp' : 0.001,
        'ls_on'      : 1.1,
        'ls_min'     : 1e-5,
        'check'     : 0,
        'delta'     : 1e-6,
        'is_plot'    : False,
        'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
    }

    fe = {
        'chunk_size' : 100000
    }

    return locals()
