# 30.05.2007, c
# last revision: 06.12.2007

fileName_mesh = 'database/tests/plane.mesh'

material_1 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'coef' : 1.0,
}
material_2 = {
    'name' : 'm',
    'mode' : 'here',
    'region' : 'Omega',
    'K' : [[1.0, 0.0], [0.0, 1.0]],
}

field_1 = {
    'name' : 'a_harmonic_field',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P2'}
}

variables = {
    't' : ('field', 'unknown', 'a_harmonic_field', (30,), 0),
    's' : ('field', 'test', 'a_harmonic_field', (30,), 't'),
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

region_1 = {
    'name' : 'Left',
    'select' : 'nodes in (x < -0.499)',
    'canCells' : True,
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.499)',
    'canCells' : True,
}
region_3 = {
    'name' : 'Gamma',
    'select' : 'nodes of surface',
    'canCells' : True,
}

ebc = {
#    'Left' : ('T3', (30,), 'linearY'),
    'Left' : ('T3', (30,), 5.0),
    'Right' : ('T3', (30,), 0.0),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

integral_2 = {
    'name' : 'i2',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d1',
}

equations = {
    'Temperature' : """dw_laplace.i1.Omega( coef, s, t ) = 0"""
#    'Temperature' : """dw_hdpm_d.Omega( m.K, s, t ) = 0"""
}

solution = {
    't' : '- 5.0 * (x - 0.5)',
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

linMin, linMax = 0.0, 2.0

##
# 31.05.2007, c
def linear( bc, ts, coor, which ):
    vals = coor[:,which]
    minVal, maxVal = vals.min(), vals.max()
    vals = (vals - minVal) / (maxVal - minVal) * (linMax - linMin) + linMin
    return vals

##
# 31.05.2007, c
def linearX( bc, ts, coor ):
    return linear( bc, ts, coor, 0 )
def linearY( bc, ts, coor ):
    return linear( bc, ts, coor, 1 )
def linearZ( bc, ts, coor ):
    return linear( bc, ts, coor, 2 )

from sfe.base.testing import TestCommon

##
# 30.05.2007, c
class Test( TestCommon ):

    ##
    # 30.05.2007, c
    def fromConf( conf, options ):
        from sfe.solvers.generic import solveStationary

        problem, vec, data = solveStationary( conf )

        test = Test( problem = problem, vec = vec, data = data,
                     conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # 30.05.2007, c
    def test_solution( self ):
        sol = self.conf.solution
        vec = self.vec
        problem  = self.problem

        ok = True
        for varName, expression in sol.iteritems():
            coor = problem.variables[varName].field.getCoor()
            anaSol = self.evalCoorExpression( expression, coor )
            numSol = problem.variables.getStatePartView( vec, varName )
            ret = self.compareVectors( anaSol, numSol,
                                       label1 = 'analytical %s' % varName,
                                       label2 = 'numerical %s' % varName )
            if not ret:
                self.report( 'variable %s: failed' % varName )

            ok = ok and ret

        return ok

    ##
    # 30.05.2007, c
    # 31.05.2007
    # 02.10.2007
    def test_boundaryFluxes( self ):
        from sfe.base.base import Struct
        from sfe.fem.evaluate import evalTermOP, BasicEvaluator
        problem  = self.problem
        vec = self.vec

        regionNames = ['Left', 'Right', 'Gamma']
        values = [5.0, -5.0, 0.0]

        getState = problem.variables.getStatePartView
        state = vec.copy()
        
        problem.timeUpdate( conf_ebc = {}, conf_epbc = {} )
#        problem.saveEBC( 'aux.vtk' )

        problem.applyEBC( state )
        ev = BasicEvaluator( problem )
        aux = ev.evalResidual( state )[0]

        field = problem.variables['t'].field

        ok = True
        for ii, regionName in enumerate( regionNames ):
            fluxTerm = 'd_hdpm_surfdvel.i2.%s( m.K, t )' % regionName
            val1 = evalTermOP( None, fluxTerm, problem )

            rvec = getState( aux, 't', True )
            reg = problem.domain.regions[regionName]
            nods = reg.getFieldNodes( field, merge = True )
            val2 = rvec[nods].sum() # Assume 1 dof per node.
        
            ok = ok and ((abs( val1 - values[ii] ) < 1e-10) and
                         (abs( val2 - values[ii] ) < 1e-10))
            self.report( '%d. %s: %e == %e == %e'\
                         % (ii, regionName, val1, val2, values[ii]) )
        
        return ok
