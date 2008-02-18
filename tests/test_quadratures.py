# 13.11.2007, c
# last revision: 18.02.2008
fileName_mesh = 'database/tests/triquad.mesh'

material_1 = {
    'name' : 'm',
    'mode' : 'here',
    'region' : 'Omega',
    'coef' : 1.0,
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}
region_1 = {
    'name' : 'Omega_1',
    'select' : 'elements of group 2',
}
region_2 = {
    'name' : 'Omega_2',
    'select' : 'elements of group 1',
}
region_3 = {
    'name' : 'Gamma_Bottom',
    'select' : 'nodes in (y < 0.00001)',
}
region_4 = {
    'name' : 'Gamma_Top',
    'select' : 'nodes in (y > 0.99999)',
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega_1' : '2_3_P2',
               'Omega_2' : '2_4_Q1'}
##     'bases' : (('Omega_1', '2_3_P1'),
##                ('Omega_2', '2_4_Q1'))
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

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Top',
    'dofs' : (30,),
    'value' : 'ebcSin',
}
ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Bottom',
    'dofs' : (30,),
    'value' : -2.0,
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2', # <quadrature name>
}

import numpy as nm
N = 2
integral_2 = {
    'name' : 'i2',
    'kind' : 'v',
    'quadrature' : 'custom', # <quadrature name>
    'vals'    : zip(nm.linspace( 1e-10, 0.5, N ),
                    nm.linspace( 1e-10, 0.5, N )),
    'weights' : [1./N] * N,
##     'vals'    : [[1./3., 1./3.]],
##     'weights' : [0.5]
}

equations = {
    'Temperature' : """dw_laplace.i2.Omega( m, s, t )
                       = 0"""
##     'Temperature' : """dw_laplace.Omega_1( m, s, t )
##                      + dw_laplace.Omega_2( m, s, t ) = 0"""
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

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

fe = {
    'chunkSize' : 1000
}

import numpy as nm
def ebcSin( bc, ts, coor ):
    val = 2 * nm.sin( coor[:,0] * 5. * nm.pi )
    return val

from sfe.base.testing import TestCommon

##
# 13.11.2007, c
class Test( TestCommon ):

    ##
    # 13.11.2007, c
    def fromConf( conf, options ):

        test = Test( conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # 13.11.2007, c
    def test_problemCreation( self ):
        from sfe.solvers.generic import solveStationary

        problem, vec, data = solveStationary( self.conf )
        ok = True
        return ok
