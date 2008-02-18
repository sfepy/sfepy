##
# c: 05.02.2008, r: 18.02.2008
fileName_mesh = 'database/simple.mesh'

material_2 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'coef' : 0.0001, # coef * \Delta t.
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'dofs' : (30,),
    'order' : 0,
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dofs' : (30,),
    'dual' : 't',
}
variable_10 = {
    'name' : 't0',
    'kind' : 'parameter field',
    'field' : 'temperature',
    'dofs' : (30,),
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
    'dofs' : (30,),
    'value' : 2.0,
}
ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Right',
    'dofs' : (30,),
    'value' : -2.0,
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'Temperature' :
    """  dw_laplace.i1.Omega( coef, s, t )
       + dw_mass_scalar.i1.Omega( s, t ) = dw_mass_scalar_r.i1.Omega( s, t0 )"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'iMax'      : 1,
    'epsA'      : 1e-10,
    'epsR'      : 1.0,
    'macheps'   : 1e-16,
    'linRed'    : 1e-2, # Linear system error < (epsA * linRed).
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
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'    : 0.0,
    't1'    : 1.0,
    'dt'    : None,
    'nStep' : 10, # has precedence over dt!
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'saveSteps' : -1,
    'variableHistory' : {'t' : 't0'},
}

fe = {
    'chunkSize' : 1000
}
