# c: 10.10.2008
import numpy as nm

filename_mesh = '../database/kostka_medium.mesh'

regions = {
    'Omega' : ('all', {}),
    'Bottom' : ('nodes in (z < -0.4999999)', {}),
    'Top' : ('nodes in (z > 0.4999999)', {}),
    'Left' : ('nodes in (x < -0.4999999)', {}),
}

field_1 = {
    'name' : 'displacement',
    'dim' : (3,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_8_Q1'}
}

field_2 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_8_Q1'}
}

variables = {
    'u'       : ('unknown field',   'displacement', 0),
    'v'       : ('test field',      'displacement', 'u'),
    'p'       : ('unknown field',   'pressure', 1),
    'q'       : ('test field',      'pressure', 'p'),
}

ebcs = {
    'fix_u' : ('Bottom', {'u.all' : 0.0}),
    'load_u' : ('Top', {'u.2' : 0.2}),
    'load_p' : ('Left', {'p.all' : 1.0}),
}

material_1 = {
    'name' : 'm',
    'region' : 'Omega',
    'values' : {
        'Lame' : {'lambda' : 1.7, 'mu' : 0.3},
        'alpha' : nm.array( [0.132, 0.132, 0.132, 0.092, 0.092, 0.092],
                            dtype = nm.float64 ),
        'K' : nm.array( [[2.0, 0.2, 0.0], [0.2, 1.0, 0.0], [0.0, 0.0, 0.5]],
                        dtype = nm.float64 ),
    }
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}

integral_2 = {
    'name' : 'i2',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}

equations = {
    'eq_1' :
    """dw_lin_elastic_iso.i2.Omega( m.Lame, v, u )
     - dw_biot.i1.Omega( m.alpha, v, p )
       = 0""",
    'eq_2' :
    """dw_biot.i1.Omega( m.alpha, u, q ) + dw_diffusion.i1.Omega( m.K, q, p )
       = 0""",
}

solver_0 = {
    'name' : 'ls_d',
    'kind' : 'ls.scipy_direct',
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

fe = {
    'chunk_size' : 10000
}
