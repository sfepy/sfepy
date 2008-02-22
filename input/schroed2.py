##
# c: 01.02.2008, r: 22.02.2008

#fileName_mesh = 'database/simple.mesh'
#fileName_mesh = 'database/phono/cube_sphere.mesh'
fileName_mesh = 'tmp/t.1.node'

options = {
    'saveEigVectors' : None,
    'squared' : False,
    'nEigs' : 10,
}

if fileName_mesh.find( 'cube_' ) >= 0:
    # Domain $Y_1$.
    region_1 = {
        'name' : 'Y1',
        'select' : 'elements of group 1',
    }

    # Domain $Y_2$.
    region_2 = {
        'name' : 'Omega',
        'select' : 'elements of group 2',
    }

    # Surface of $Y_2$.
    region_100 = {
        'name' : 'Surface',
        'select' : 'r.Y1 *n r.Omega',
        'canCells' : False,
    }
else:
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

field_0 = {
    'name' : 'field_Psi',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}

variable_1 = {
    'name' : 'Psi',
    'kind' : 'unknown field',
    'field' : 'field_Psi',
    'dofs' : (10,),
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : 'field_Psi',
    'dofs' : (10,),
    'dual' : 'Psi',
}
variable_3 = {
    'name' : 'V',
    'kind' : 'parameter field',
    'field' : 'field_Psi',
    'dofs' : (10,),
}
variable_4 = {
    'name' : 'n',
    'kind' : 'parameter field',
    'field' : 'field_Psi',
    'dofs' : (10,),
}

ebc_1 = {
    'name' : 'ZeroSurface',
    'region' : 'Surface',
    'dofs' : (10,),
    'value' : 0.0,
}

equations = {
    'lhs' : """  dw_laplace.i1.Omega( m.val, v, Psi )
               + dw_mass_scalar_variable.i1.Omega( matV.V, v, Psi )""",
    'rhs' : """dw_mass_scalar.i1.Omega( v, Psi )""",
}

equations_vh = {
    'poisson' : """dw_laplace.i1.Omega( m.one, v, Psi )
                   = dw_mass_scalar_r.i1.Omega( v, n )"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_cg',
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

fe = {
    'chunkSize' : 100000
}

##
# c: 01.02.2008, r: 22.02.2008
def funV( ts, coor, region, ig, mode = None ):
    import numpy as nm
    out = {}
    C = 0.5
    r = nm.sqrt( coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2 )
    val = nm.array( - C * 5.0 / r, ndmin = 3 )
    #val = nm.zeros_like( val )
    out['V'] = val
    return out
