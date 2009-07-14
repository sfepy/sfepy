# 13.11.2007, c
# last revision: 25.02.2008
filename_mesh = '../database/tests/triquad.mesh'

material_1 = {
    'name' : 'm',
    'region' : 'Omega',
    'values' : {'val' : 1.0},
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
    'order' : 0,
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 't',
}

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Top',
    'dofs' : {'t.0' : 'ebc_sin'},
}
ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Bottom',
    'dofs' : {'t.0' : -2.0},
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
    'Temperature' : """dw_laplace.i2.Omega( m.val, s, t )
                       = 0"""
##     'Temperature' : """dw_laplace.Omega_1( m, s, t )
##                      + dw_laplace.Omega_2( m, s, t ) = 0"""
}

solver_0 = {
    'name' : 'ls',
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

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

fe = {
    'chunk_size' : 1000
}

import numpy as nm
def ebc_sin(ts, coor, bc):
    val = 2 * nm.sin( coor[:,0] * 5. * nm.pi )
    return val

functions = {
    'ebc_sin' : (ebc_sin,),
}


from sfepy.base.testing import TestCommon

##
# 13.11.2007, c
class Test( TestCommon ):

    ##
    # 13.11.2007, c
    def from_conf( conf, options ):

        test = Test( conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # 13.11.2007, c
    def test_problem_creation( self ):
        from sfepy.solvers.generic import solve_stationary

        problem, vec, data = solve_stationary( self.conf )
        ok = True
        return ok
