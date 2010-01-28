from sfepy import top_dir

filename_mesh = top_dir + '/meshes/2d/special/two_rectangles.mesh'

v_n = 1.0 # m/s
w = 1000.0
c = 343.0 # m/s
rho = 1.55 # kg/m^3

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

materials = {
    'one' : ('Omega', {'one' : 1.0} ),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_in' : ('nodes in (x < 0.01)', {}),
    'Gamma_out' : ('nodes in (x > 0.99)', {}),
}

fields = {
    'accoustic_pressure' : ((1,1), 'complex', 'Omega', {'Omega' : '2_4_Q1'}),
}

variables = {
    'p'   : ('unknown field', 'accoustic_pressure', 0, 'previous'),
    'q'   : ('test field',    'accoustic_pressure', 'p'),
}

ebcs = {
}

integrals = {
    'ivol' : ('v', 'gauss_o2_d2'),
    'isurf' : ('s', 'gauss_o2_d1'),
}

equations = {
    'Acoustic pressure' :
    """%s * dw_laplace.ivol.Omega( one.one, q, p )
    - %s * dw_mass_scalar.ivol.Omega( q, p )
    - %s * dw_surface_mass_scalar.isurf.Gamma_out( q, p )
    = %s * dw_surface_integrate.isurf.Gamma_in( q )"""
    % (c*c, w*w, 1j*w*c, 1j*w*c*c*rho*v_n)
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {'i_max'      : 1,
                               'eps_a'      : 1e-1,
                               'eps_r'      : 1.0,
                               'macheps'   : 1e-16,
                               # Linear system error < (eps_a * lin_red).
                               'lin_red'    : 1e-1,
                               'ls_red'     : 0.1,
                               'ls_red_warp' : 0.001,
                               'ls_on'      : 1.1,
                               'ls_min'     : 1e-5,
                               'check'     : 0,
                               'delta'     : 1e-6,
                               'is_plot'    : False,
                               # 'nonlinear' or 'linear' (ignore i_max)
                               'problem'   : 'nonlinear',
                               } )
}

fe = {
    'chunk_size' : 100000
}
