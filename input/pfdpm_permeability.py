# 24.05.2007, c
# last revision: 25.02.2008
from sfepy.fem.periodic import *

#fileName_mesh = 'database/micro/perf_symm638t.mesh'
fileName_mesh = 'database/micro/perf_symm944t.mesh'
#fileName_mesh = 'database/micro/perf_symm748t.mesh'

if fileName_mesh.find( 'symm' ):
    region_1 = {
        'name' : 'Y1',
        'select' : """elements of group 3""",
    }
    region_2 = {
        'name' : 'Y2',
        'select' : """elements of group 4 +e elements of group 6
                      +e elements of group 8""",
    }
    region_4 = {
        'name' : 'Y1Y2',
        'select' : """r.Y1 +e r.Y2""",
    }
    region_5 = {
        'name' : 'Walls',
        'select' : """r.EBCGamma1 +n r.EBCGamma2""",
    }
    region_310 = {
        'name' : 'EBCGamma1',
        'select' : """(elements of group 1 *n elements of group 3)
                      +n
                      (elements of group 2 *n elements of group 3)
                      """,
    }
    region_320 = {
        'name' : 'EBCGamma2',
        'select' : """(elements of group 5 *n elements of group 4)
                      +n
                      (elements of group 1 *n elements of group 4)
                      +n
                      (elements of group 7 *n elements of group 6)
                      +n
                      (elements of group 2 *n elements of group 6)
                      +n
                      (elements of group 9 *n elements of group 8)
                      +n
                      (elements of group 2 *n elements of group 8)
                      """,
    }


w2 = 0.499
# Sides.
region_20 = {
    'name' : 'Left',
    'select' : 'nodes in (x < %.3f)' % -w2,
}
region_21 = {
    'name' : 'Right',
    'select' : 'nodes in (x > %.3f)' % w2,
}
region_22 = {
    'name' : 'Bottom',
    'select' : 'nodes in (y < %.3f)' % -w2,
}
region_23 = {
    'name' : 'Top',
    'select' : 'nodes in (y > %.3f)' % w2,
}

field_1 = {
    'name' : '2_velocity',
    'dim' : (2,1),
    'flags' : (),
    'domain' : 'Y1Y2',
    'bases' : {'Y1Y2' : '2_3_P2'}
}

field_2 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Y1Y2',
    'bases' : {'Y1Y2' : '2_3_P1'}
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '2_velocity',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '2_velocity',
    'dual' : 'u',
}
variable_3 = {
    'name' : 'p',
    'kind' : 'unknown field',
    'field' : 'pressure',
    'order' : 1,
}
variable_4 = {
    'name' : 'q',
    'kind' : 'test field',
    'field' : 'pressure',
    'dual' : 'p',
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

equations = {
    'namespaces' : {
    'dw_' : ('div', 'grad', 'div_grad', 'convect'),
    },
    'balance' :
    """+ div_grad.i1.Y1Y2( fluid.viscosity, v, u ) - grad.i1.Y1Y2( v, p ) = 0""",
    'incompressibility' :
    """div.i1.Y1Y2( q, u ) = 0""",
}

material_1 = {
    'name' : 'fluid',
    'mode' : 'here',
    'region' : 'Y1Y2',
    'viscosity' : 1.0,
    'density' : 1e0,
}

ebc_1 = {
    'name' : 'walls',
    'region' : 'Walls',
    'dofs' : {'u.all' : 0.0},
}
ebc_2 = {
    'name' : 'top_velocity',
    'region' : 'Top',
    'dofs' : {'u.1' : -1.0, 'u.0' : 0.0},
}
ebc_10 = {
    'name' : 'bottom_pressure',
    'region' : 'Bottom',
    'dofs' : {'p.0' : 0.0},
}

epbc_1 = {
    'name' : 'u_rl',
    'region' : ['Left', 'Right'],
    'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
    'match' : 'matchYLine',
}

##
# FE assembling parameters.
fe = {
    'chunkSize' : 100,
    'cacheOverride' : True,
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'iMax'      : 2,
    'epsA'      : 1e-8,
    'epsR'      : 1e-2,
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

saveFormat = 'hdf5' # 'hdf5' or 'vtk'
