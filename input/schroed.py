##
# c: 01.02.2008, r: 25.02.2008

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
               + dw_mass_scalar_variable.i1.Omega( matV.V, v, Psi )""",
    'rhs' : """dw_mass_scalar.i1.Omega( v, Psi )""",
}

fe = {
    'chunkSize' : 100000
}

##
# c: 01.02.2008, r: 01.02.2008
def funV( ts, coor, region, ig, mode = None ):
    import numpy as nm
    out = {}
    C = 0.5
    val = nm.array( C* (coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2), ndmin = 3 )
    #val = nm.zeros_like( val )
    out['V'] = val
    return out
