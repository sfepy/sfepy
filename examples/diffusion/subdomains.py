"""
\Delta p_1 = g_1 in Omega_1
\Delta p_2 = g_2 in Omega_2
p_1 = p_2 + c on \Gamma_{12}
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/special/two_rectangles.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

regions = {
    'Omega1' : ('elements of group 1', {}),
    'Omega2' : ('elements of group 2', {}),
    'Gamma1' : ('nodes in (y > 0.29)', {}),
    'Gamma2' : ('nodes in (y < 0.01)', {}),
    'Gamma12' : ('r.Omega1 *n r.Omega2', {}),
    'Gamma12_1' : ('copy r.Gamma12', {'forbid' : 'group 2'}),
    'Gamma12_2' : ('copy r.Gamma12', {'forbid' : 'group 1'}),
}

materials = {
    'one1' : ('Omega1', {'one' : 1.0} ),
    'one2' : ('Omega2', {'one' : 1.0} ),
    'jump1'  : ('Gamma12_1', {'val' : 1.0} ),
    'jump2'  : ('Gamma12_2', {'val' : 1.0} ),
}

fields = {
    'pressure1' : ((1,1), 'real', 'Omega1', {'Omega1' : '2_4_Q1'}),
    'pressure2' : ((1,1), 'real', 'Omega2', {'Omega2' : '2_4_Q1'}),
}

variables = {
    'p1'   : ('unknown field', 'pressure1', 0),
    'q1'   : ('test field',    'pressure1', 'p1'),
    'p2'   : ('unknown field', 'pressure2', 1),
    'q2'   : ('test field',    'pressure2', 'p2'),
}

ebcs = {
    'p1' : ('Gamma1', {'p1.all' : 0.0}),
    'p2' : ('Gamma2', {'p2.0' : 1.0}),
}

integrals = {
    'ivol' : ('v', 'gauss_o2_d2'),
    'isurf' : ('s2', 'gauss_o2_d1'),
}

equations = {
    'eq_1' :
    """dw_laplace.ivol.Omega1( one1.one, q1, p1 )
     + dw_jump.isurf.Gamma12_1( jump1.val, q1, p1, tr(p2) )
     = 0""",
    'eq_2' :
    """dw_laplace.ivol.Omega2( one2.one, q2, p2 )
     + dw_jump.isurf.Gamma12_2( jump2.val, q2, tr(p1), p2 )
     = 0""",
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
