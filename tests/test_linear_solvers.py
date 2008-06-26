# c: 02.05.2008, r: 02.05.2008
fileName_mesh = 'database/phono/cube_cylinder.mesh'

material_2 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'val' : 1.0,
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_Left' : ('nodes in (x < 0.0001)', {}),
    'Gamma_Right' : ('nodes in (x > 0.999)', {}),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}

equations = {
    'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
}

solver_100 = {
    'name' : 'dls100',
    'kind' : 'ls.umfpack',
}

solver_101 = {
    'name' : 'ls101',
    'kind' : 'ls.pyamg',

    'method' : 'ruge_stuben_solver',
    'epsA'   : 1e-12,
}

solver_102 = {
    'name' : 'ls102',
    'kind' : 'ls.pyamg',

    'method' : 'smoothed_aggregation_solver',
    'epsA'   : 1e-12,
}

solver_110 = {
    'name' : 'ls110',
    'kind' : 'ls.scipy_iterative',

    'method' : 'cg',
    'iMax'   : 1000,
    'epsA'   : 1e-12,
}

solver_111 = {
    'name' : 'ls111',
    'kind' : 'ls.scipy_iterative',

    'method' : 'bicgstab',
    'iMax'   : 1000,
    'epsA'   : 1e-12,
}

solver_112 = {
    'name' : 'ls112',
    'kind' : 'ls.scipy_iterative',

    'method' : 'qmr',
    'iMax'   : 1000,
    'epsA'   : 1e-12,
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
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

options = {
    'nls' : 'newton',
}

fe = {
    'chunkSize' : 100000
}

from sfepy.base.testing import TestCommon
outputName = 'test_linear_solvers_%s.vtk'

##
# c: 02.05.2008
class Test( TestCommon ):
    canFail = ['ls.pyamg']

    ##
    # c: 02.05.2008, r: 02.05.2008
    def fromConf( conf, options ):
        from sfepy.fem.problemDef import ProblemDefinition

        problem = ProblemDefinition.fromConf( conf )
        problem.timeUpdate()

        test = Test( problem = problem, 
                     conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # c: 02.05.2008, r: 02.05.2008
    def _listLinearSolvers( self, confs ):
        d = []
        for key, val in confs.iteritems():
            if val.kind.find( 'ls.' ) == 0:
                d.append( val )
        d.sort( cmp = lambda a, b: cmp( a.name, b.name ) )

        return d

    ##
    # c: 02.05.2008, r: 07.05.2008
    def test_solvers( self ):
        from sfepy.solvers.generic import solveStationary
        from sfepy.base.base import IndexedStruct
        import os.path as op

        solverConfs = self._listLinearSolvers( self.problem.solverConfs )

        ok = True
        tt = []
        for solverConf in solverConfs:
            if hasattr( solverConf, 'method' ):
                method = solverConf.method
            else:
                method = ''
            name = ' '.join( (solverConf.name, solverConf.kind, method) )
            self.report( name )
            self.report( 'matrix size:', self.problem.mtxA.shape )
            self.report( '        nnz:', self.problem.mtxA.nnz )
            status = IndexedStruct()
            try:
                self.problem.initSolvers( nlsStatus = status,
                                          lsConf = solverConf )
                state = self.problem.solve()
                failed = status.condition != 0
##                 self.problem.mtxA.save( 'mtx_laplace_cube',
##                                         format='%d %d %.12e\n' )
            except Exception, exc:
                failed = True
                status = None

            ok = ok and ((not failed) or (solverConf.kind in self.canFail))

            if status is not None:
                for kv in status.timeStats.iteritems():
                    self.report( '%10s: %7.2f [s]' % kv )
                self.report( 'condition: %d, err0: %.3e, err: %.3e'\
                             % (status.condition, status.err0, status.err) )
                tt.append( [name, status.timeStats['solve']] )

                fname = op.join( self.options.outDir,
                                op.split( self.conf.outputName )[1] ) % name
                self.problem.saveState( fname, state )
            else:
                self.report( 'solver failed:' )
                self.report( exc )
                tt.append( [name, 1e10] )


        tt.sort( cmp = lambda a, b: cmp( a[1], b[1] ) )
        self.report( 'solution times:' )
        for row in tt:
            self.report( '%.2f' % row[1], ':', row[0] )

        return ok
