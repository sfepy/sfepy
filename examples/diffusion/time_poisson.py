##
# c: 05.02.2008
from sfepy import top_dir

filename_mesh = top_dir + '/meshes/3d/cylinder.mesh'

t0 = 0.0
t1 = 0.1
n_step = 11

material_2 = {
    'name' : 'coef',
    'region' : 'Omega',
    'values' : {'val' : 0.01},
    'kind' : 'stationary', # 'stationary' or 'time-dependent'
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variable_1 = {
    'name' : 'T',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0,
    'history' : 'previous', # 'previous' or 'full'
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 'T',
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_Left' : ('nodes in (x < 0.00001)', {}),
    'Gamma_Right' : ('nodes in (x > 0.099999)', {}),
}

ebcs = {
    'T1': ('Gamma_Left', {'T.0' : 2.0}),
    'T2': ('Gamma_Right', {'T.0' : -2.0}),
}

def get_ic(coor, ic):
    """Non-constant initial condition."""
    import numpy as nm
    # Normalize x coordinate.
    mi, ma = coor[:,0].min(), coor[:,0].max()
    nx = (coor[:,0] - mi) / (ma - mi)
    return nm.where( (nx > 0.25) & (nx < 0.75 ), 8.0 * (nx - 0.5), 0.0 )

functions = {
    'get_ic' : (get_ic,),
}
    
ics = {
    'ic' : ('Omega', {'T.0' : 'get_ic'}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}

equations = {
    'Temperature' :
    """dw_mass_scalar.i1.Omega( s, dT/dt )
     + dw_laplace.i1.Omega( coef.val, s, T ) = 0"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',

    'presolve' : True,
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
    'problem'   : 'linear', # 'nonlinear' or 'linear' (ignore i_max)
}

solver_2 = {
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'    : t0,
    't1'    : t1,
    'dt'    : None,
    'n_step' : n_step, # has precedence over dt!
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
}

fe = {
    'chunk_size' : 1000
}
