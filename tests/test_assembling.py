# c: 14.04.2008, r: 14.04.2008
from sfepy.fem.periodic import *

fileName_mesh = 'database/tests/plane.mesh'

material_2 = {
    'name' : 'm',
    'mode' : 'here',
    'region' : 'Omega',
    'K' : [[3.0, 0.1], [0.3, 1.0]],
}

field_1 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P2'}
}

variables = {
    'p' : ('unknown field', 'pressure', 0),
    'q' : ('test field', 'pressure', 'p'),
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

region_1 = {
    'name' : 'Left',
    'select' : 'nodes in (x < -0.499)',
}
region_10 = {
    'name' : 'LeftStrip',
    'select' : 'nodes in (x < -0.499) & (y > -0.199) & (y < 0.199)',
}
region_11 = {
    'name' : 'LeftFix',
    'select' : 'r.Left -n r.LeftStrip',
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.499)',
}
region_20 = {
    'name' : 'RightStrip',
    'select' : 'nodes in (x > 0.499) & (y > -0.199) & (y < 0.199)',
}
region_21 = {
    'name' : 'RightFix',
    'select' : 'r.Right -n r.RightStrip',
}

ebcs = {
    't_left' : ('LeftFix', {'p.0' : 5.0}),
    't_right' : ('RightFix', {'p.0' : 0.0}),
}

epbc_10 = {
    'name' : 'periodic_x',
    'region' : ['LeftStrip', 'RightStrip'],
    'dofs' : {'p.0' : 'p.0'},
    'match' : 'matchYLine',
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

equations = {
    'eq' : """dw_diffusion.i1.Omega( m.K, q, p ) = 0"""
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
    'linSolver' : 'umfpack',
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

fe = {
    'chunkSize' : 1000
}

from sfepy.base.testing import TestCommon

##
# c: 14.04.2008
class Test( TestCommon ):

    ##
    # c: 14.04.2008, r: 14.04.2008
    def fromConf( conf, options ):
        import os.path as op
        from sfepy.solvers.generic import solveStationary

        problem, vec, data = solveStationary( conf )
        name = op.join( options.outDir,
                        op.splitext( op.basename( __file__ ) )[0] + '.vtk' )
        problem.saveState( name, vec )

        test = Test( problem = problem, vec = vec, data = data,
                     conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    
    ##
    # c: 14.04.2008, r: 14.04.2008
    def test_vector_matrix( self ):
        from sfepy.fem.evaluate import evalTermOP
        problem  = self.problem

        state = problem.createStateVector()
        problem.applyEBC( state )

        aux1 = evalTermOP( state, "dw_diffusion.i1.Omega( m.K, q, p )",
                           problem, dwMode = 'vector' )
        aux1g = problem.variables.makeFullVec( aux1 )
        problem.timeUpdate( conf_ebc = {}, conf_epbc = {} )
        mtx = evalTermOP( state, "dw_diffusion.i1.Omega( m.K, q, p )",
                          problem, dwMode = 'matrix' )
        aux2g = mtx * state
        problem.timeUpdate( conf_ebc = self.conf.ebcs,
                            conf_epbc = self.conf.epbcs )
        aux2 = problem.variables.stripStateVector( aux2g, followEPBC = True )

        ret = self.compareVectors( aux1, aux2,
                                   label1 = 'vector mode',
                                   label2 = 'matrix mode' )
        if not ret:
            self.report( 'variable %s: failed' % varName )

        return ret
