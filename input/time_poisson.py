##
# c: 05.02.2008
filename_mesh = 'database/simple.mesh'

from sfepy.solvers.ts import TimeStepper
t0 = 0.0
t1 = 0.1
n_step = 11
ts = TimeStepper( t0, t1, None, n_step )

material_2 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'val_dt' : 0.01 * ts.dt, # coef * \Delta t.
    'kind' : 'stationary', # 'stationary' or 'time-dependent'
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0,
    'history' : 'previous', # 'previous' or 'full'
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 't',
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}
region_03 = {
    'name' : 'Gamma_Left',
    'select' : 'nodes in (x < 0.00001)',
}
region_4 = {
    'name' : 'Gamma_Right',
    'select' : 'nodes in (x > 0.099999)',
}

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Left',
    'dofs' : {'t.0' : 2.0},
    'value' : 2.0,
}
ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Right',
    'dofs' : {'t.0' : -2.0},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'Temperature' :
    """  dw_laplace.i1.Omega( coef.val_dt, s, t )
       + dw_mass_scalar.i1.Omega( s, t )
       = dw_mass_scalar.i1.Omega( s, t[-1] )"""
}

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
