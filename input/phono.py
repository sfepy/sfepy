# 25.09.2007, c
# last revision: 20.02.2008
"""
u1 is a dummy variable used unly for volume computation.
"""

#fileName_mesh = 'database/phono/cube_sphere.mesh'
#fileName_mesh = 'database/phono/cube_cylinder.mesh'
fileName_mesh = 'database/phono/mesh_circ21.mesh'
#fileName_mesh = 'database/phono/mesh_circ21_small.mesh'

options = {
    'saveEigVectors' : (10, 10),
    'eigRange' : (0, 30), # -> freqRange = eigs[slice(*eigRange)][[0, -1]]
    'freqMargins' : (10, 10), # % of freqRange
    'feps' : 1e-10, # frequency
    'zeps' : 1e-12, # zero finding
    'teps' : 1e-3, # eigenmomentum threshold
    'freqStep' : 0.01, # % of freqRange
#    'eigVectorTransform' : ('selectInPlane', 'z', 1e-1),
#    'plotTranform' : ('clip', (-20, 20)),
    'plotTranform' : ('normalize', (-1, 1)),
    'squared' : False,
}

# Whole domain $Y$.
region_1000 = {
    'name' : 'Y',
    'select' : 'all',
}

# Domain $Y_1$.
region_1 = {
    'name' : 'Y1',
    'select' : 'elements of group 1',
}

# Domain $Y_2$.
region_2 = {
    'name' : 'Y2',
    'select' : 'elements of group 2',
}

# Surface of $Y_2$.
region_100 = {
    'name' : 'Y2_Surface',
    'select' : 'r.Y1 *n r.Y2',
    'canCells' : False,
}

material_1 = {
    'name' : 'matrix',
    'mode' : 'here',
    'region' : 'Y1',

    # aluminium
    'lame' : {'lambda' : 5.898, 'mu' : 2.681}, # in 1e+10 Pa
    'density' : 0.2799, # in 1e4 kg/m3
}

material_2 = {
    'name' : 'inclusion',
    'mode' : 'here',
    'region' : 'Y2',

    # epoxy
    'lame' : {'lambda' : 0.1798, 'mu' : 0.148}, # in 1e+10 Pa
    'density' : 0.1142, # in 1e4 kg/m3
}

if fileName_mesh.find( 'cube_' ) >= 0:
    dim = 3
else:
    dim = 2

if dim == 3:
    field_0 = {
        'name' : '3_displacement_Y1',
        'dim' : (3,1),
        'flags' : (),
        'domain' : 'Y1',
        'bases' : {'Y1' : '3_4_P1'}
    }

    field_1 = {
        'name' : '3_displacement_Y2',
        'dim' : (3,1),
        'flags' : (),
        'domain' : 'Y2',
        'bases' : {'Y2' : '3_4_P1'}
    }

    field_2 = {
        'name' : 'eigenDirection',
        'dim' : (1,1),
        'flags' : (),
        'domain' : 'Y2',
        'bases' : {'Y2' : '3_4_P1'}
    }

    variable_1 = {
        'name' : 'u',
        'kind' : 'unknown field',
        'field' : '3_displacement_Y2',
        'dofs' : (0, 1, 2),
        'order' : 0,
    }
    variable_2 = {
        'name' : 'v',
        'kind' : 'test field',
        'field' : '3_displacement_Y2',
        'dofs' : (0, 1, 2),
        'dual' : 'u',
    }
    variable_3 = {
        'name' : 'u1',
        'kind' : 'parameter field',
        'field' : '3_displacement_Y1',
        'dofs' : (0, 1, 2),
    }
    variable_4 = {
        'name' : 'uc',
        'kind' : 'parameter field',
        'field' : 'eigenDirection',
        'dofs' : (3,),
    }
    variable_5 = {
        'name' : 'd',
        'kind' : 'parameter field',
        'field' : 'eigenDirection',
        'dofs' : (3,),
    }

    ebc_1 = {
        'name' : 'ZeroSurface',
        'region' : 'Y2_Surface',
        'dofs' : (0, 1, 2),
        'value' : 0.0,
    }

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o2_d3',
    }

else:
    field_0 = {
        'name' : '2_displacement_Y1',
        'dim' : (2,1),
        'flags' : (),
        'domain' : 'Y1',
        'bases' : {'Y1' : '2_3_P1'}
    }

    field_1 = {
        'name' : '2_displacement_Y2',
        'dim' : (2,1),
        'flags' : (),
        'domain' : 'Y2',
        'bases' : {'Y2' : '2_3_P1'}
    }

    field_2 = {
        'name' : 'eigenDirection',
        'dim' : (1,1),
        'flags' : (),
        'domain' : 'Y2',
        'bases' : {'Y2' : '2_3_P1'}
    }

    variable_1 = {
        'name' : 'u',
        'kind' : 'unknown field',
        'field' : '2_displacement_Y2',
        'dofs' : (0, 1),
        'order' : 0,
    }
    variable_2 = {
        'name' : 'v',
        'kind' : 'test field',
        'field' : '2_displacement_Y2',
        'dofs' : (0, 1),
        'dual' : 'u',
    }
    variable_3 = {
        'name' : 'u1',
        'kind' : 'parameter field',
        'field' : '2_displacement_Y1',
        'dofs' : (0, 1),
    }
    variable_4 = {
        'name' : 'uc',
        'kind' : 'parameter field',
        'field' : 'eigenDirection',
        'dofs' : (3,),
    }
    variable_5 = {
        'name' : 'd',
        'kind' : 'parameter field',
        'field' : 'eigenDirection',
        'dofs' : (3,),
    }

    ebc_1 = {
        'name' : 'ZeroSurface',
        'region' : 'Y2_Surface',
        'dofs' : (0, 1),
        'value' : 0.0,
    }

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o2_d2',
    }

equations = {
    'lhs' : """dw_sdcc.i1.Y2( inclusion.lame, v, u )""",
    'rhs' : """dw_mass_vector.i1.Y2( inclusion.density, v, u )""",
}

##
# FE assembling parameters.
fe = {
    'chunkSize' : 1000
}

import numpy as nm
def clip( data, plotRange ):
    return nm.clip( data, *plotRange )

def normalize( data, plotRange ):
    aux = nm.arctan( data )
    return clip( aux, plotRange )

##
# 02.10.2007, c
def selectInPlane( vec, shape, normalDirection, eps ):
    nNod, dim = shape
    dirVecs = {2 : {'x': 0, 'y' : 1, 'z' : 1},
               3 : {'x': 0, 'y' : 1, 'z' : 2}}
    ident = nm.eye( dim, dtype = nm.float64 )
    dirVec = ident[:,dirVecs[dim][normalDirection]]

    proj = nm.dot( nm.reshape( vec, (nNod, dim) ), dirVec )
    if nm.any( nm.abs( proj ) > eps ):
        return nm.zeros_like( vec ), True
    else:
        return vec, False
