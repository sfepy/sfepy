# c: 07.05.2007, r: 08.07.2008
fileName_mesh = 'database/phono/mesh_circ21.mesh'

dim = 2

field_1 = {
    'name' : 'a_harmonic_field',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P1'}
}

variables = {
    't': ('unknown field', 'a_harmonic_field', 0),
    's': ('test field',    'a_harmonic_field', 't'),
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001) & (y < 0.001)', {}),
    'Right' : ('nodes in (x > 0.999)', {}),
    'Gamma' : ('nodes of surface', {}),
}

ebcs = {
    't_left' : ('Gamma', {'t.0' : 'ebc'}),
#    't_right' : ('Right', {'t.0' : 'ebc'}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

material_1 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'val' : 12.0,
    'K' : [[1.0, 0.3], [0.3, 2.0]]
}

material_2 = {
    'name' : 'rhs',
    'mode' : 'function',
    'region' : 'Omega',
    'function' : 'rhs',
    'extraArgs' : {'expression' : None},
}

equations = {
    'Laplace' :
    """2 * dw_laplace.i1.Omega( coef.val, s, t )
    """,
    'Diffusion' :
    """3 * dw_diffusion.i1.Omega( coef.K, s, t )
    """,
}
equations_rhs = {
    'Laplace' :
    """= - dw_volume_lvf.i1.Omega( rhs.val, s )""",
    'Diffusion' :
    """= - dw_volume_lvf.i1.Omega( rhs.val, s )""",
}

val = material_1['val']
solutions = {
    'sincos' : ('t', 'sin( 3.0 * x ) * cos( 4.0 * y )'),
    'poly' : ('t', '(x**2) + (y**2)'),
    'polysin' : ('t', '((x - 0.5)**3) * sin( 5.0 * y )'),
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
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

fe = {
    'chunkSize' : 100000
}

import numpy as nm
try:
    import sympy_operators as sops
except ImportError, exc:
    sops = None
    
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause
outputName = 'test_msm_symbolic_%s.vtk'

##
# c: 07.05.2007, r: 09.05.2008
solution = ['']
def ebc( bc, ts, coor, solution = solution ):
    expression = solution[0]
    val = TestCommon.evalCoorExpression( expression, coor )
    return nm.atleast_1d( val )

##
# c: 07.05.2007, r: 09.05.2008
def rhs( ts, coor, region, ig, expression = None ):
    if expression is None:
        expression = '0.0 * x'

    val = TestCommon.evalCoorExpression( expression, coor )
    return {'val' : nm.atleast_1d( val )}

##
# c: 07.05.2008
class Test( TestCommon ):

    ##
    # c: 07.05.2007, r: 25.06.2008
    def fromConf( conf, options ):
        from sfepy.fem.problemDef import ProblemDefinition

        problem = ProblemDefinition.fromConf( conf, initEquations = False )
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # c: 09.05.2007, r: 08.07.2008
    def _buildRHS( self, equation, sols ):
        from sfepy.fem.equations import buildArgs

        problem  = self.problem
        rhss = {}
        self.report( '%s:' % equation.name )
        self.report( 'evaluating terms, "<=" is solution, "=>" is the rhs:' )
        for term in equation.terms:
            if not hasattr( term, 'symbolic' ):
                self.report( 'term %s has no symbolic description!' % term.name )
                raise ValueError
            expr = term.symbolic['expression']
            argMap = term.symbolic['map']
            self.report( '%s( %s )' %\
                         (term.name, ', '.join( term.argTypes )) )
            self.report( 'multiplicator: %f' % term.sign )
            self.report( '  symbolic:', expr )
            self.report( '  using argument map:', argMap )
            args = buildArgs( term, problem.variables, problem.materials )
            for solName, sol in sols.iteritems():
                rhs = self._evalTerm( sol[1], term, args, sops )
                srhs = "(%s * (%s))" % (term.sign, rhs)
                rhss.setdefault( solName, [] ).append( srhs )

        for key, val in rhss.iteritems():
            rhss[key] = '+'.join( val )

        return rhss

    ##
    # c: 09.05.2007, r: 25.06.2008
    def _evalTerm( self, sol, term, args, sops ):
        """Works for scalar, single unknown terms only!"""
        expr = term.symbolic['expression']
        argMap = term.symbolic['map']
        env = {'x' : sops.Symbol( 'x' ),
               'y' : sops.Symbol( 'y' ),
               'z' : sops.Symbol( 'z' ),
               'dim' : dim}
        for key, val in argMap.iteritems():
            if val == 'state':
                env[key] = sol
            else:
                env[key] = term.getArgs( [val], **args )[0]

            if val[:8] == 'material':
                aux = nm.array( env[key] )
                if nm.prod( aux.shape ) == 1:
                    env[key] = aux.squeeze()
                else:
                    import sympy
                    env[key] = sympy.Matrix( aux )

#        print env

        self.report( '  <= ', sol )
        sops.set_dim( dim )
        val = eval( expr, sops.__dict__, env ).tostr()
        self.report( '   =>', val )
        return val

    ##
    # c: 07.05.2007, r: 30.06.2008
    def _test_msm_symbolic( self, equations ):
        import os.path as op

        if sops is None:
            self.report( 'cannot import sympy, skipping' )
            return True

        problem  = self.problem

        # update data so that buildArgs() works...
        matArgs = {'rhs' : {'expression' : '0 * x'}} 
        problem.updateMaterials( extraMatArgs = matArgs )

        ok = True
        for eqName, equation in equations.iteritems():
            problem.setEquations( {eqName : equation} )
            rhss = self._buildRHS( problem.equations[eqName],
                                   self.conf.solutions )
            erhs = problem.conf.equations_rhs[eqName]  
            problem.setEquations( {eqName : equation + erhs} )
            for solName, sol in problem.conf.solutions.iteritems():
                self.report( 'testing', solName )
                varName, solExpr = sol
                rhsExpr = rhss[solName]

                self.report( 'sol:', solExpr )
                self.report( 'rhs:', rhsExpr )
                matArgs = {'rhs' : {'expression' : rhsExpr}} 
                globals()['solution'][0] = solExpr
                problem.timeUpdate( extraMatArgs = matArgs )
                problem.equations.resetTermCaches()
                vec = problem.solve()
                coor = problem.variables[varName].field.getCoor()
                anaSol = self.evalCoorExpression( solExpr, coor )
                numSol = problem.variables.getStatePartView( vec, varName )

                anaNorm = nm.linalg.norm( anaSol, nm.inf )
                ret = self.compareVectors( anaSol, numSol,
                                           allowedError = anaNorm * 1e-2,
                                           label1 = 'analytical %s' % varName,
                                           label2 = 'numerical %s' % varName,
                                           norm = nm.inf )
                if not ret:
                    self.report( 'variable %s: failed' % varName )

                fname = op.join( self.options.outDir, self.conf.outputName )
                out = {}
                aux = problem.stateToOutput( anaSol )
                out['ana_t'] = aux['t']
                aux = problem.stateToOutput( numSol )
                out['num_t'] = aux['t']

                problem.domain.mesh.write( fname % '_'.join( (solName, eqName) ),
                                           io = 'auto', out = out )

                ok = ok and ret

        return ok

    ##
    # c: 30.06.2008, r: 30.06.2008
    def _getEquations( self, name ):
        """Choose a sub-problem from all equations."""
        return {name : self.problem.conf.equations[name]}
        
    ##
    # c: 30.06.2008, r: 30.06.2008
    def test_msm_symbolic_laplace( self ):
        return self._test_msm_symbolic( self._getEquations( 'Laplace' ) )

    ##
    # c: 30.06.2008, r: 30.06.2008
    def test_msm_symbolic_diffusion( self ):
        return self._test_msm_symbolic( self._getEquations( 'Diffusion' ) )
