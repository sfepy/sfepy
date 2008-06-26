from sfepy.base.base import *
from sfepy.solvers.solvers import LinearSolver

import scipy
if scipy.version.version < '0.7.0.dev3861':
    import scipy.linsolve.umfpack as um
else:
    try:
        if scipy.version.version < '0.7.0.dev3998':
            import scipy.splinalg.dsolve.umfpack as um
        else:
            import scipy.sparse.linalg.dsolve.umfpack as um
    except ImportError:
        import scikits.umfpack as um

um.configure( assumeSortedIndices = True )

##
# 10.10.2007, c
class Umfpack( LinearSolver ):
    name = 'ls.umfpack'

    ##
    # c: 10.10.2007, r: 06.03.2008
    def __init__( self, conf, **kwargs ):
        LinearSolver.__init__( self, conf, **kwargs )
        self.umfpack = um.UmfpackContext()

        if self._presolve() and hasattr( self, 'mtx' ):
            if self.mtx is not None:
                self.umfpack.numeric( self.mtx )

    ##
    # 02.12.2005, c
    # 14.04.2006
    # 02.10.2007
    # 10.10.2007, from solve_umfpack()
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        conf = getDefault( conf, self.conf )
        mtx = getDefault( mtx, self.mtx )
        status = getDefault( status, self.status )

##     umfpack.control[um.UMFPACK_PRL] = 4
##     umfpack.control[um.UMFPACK_IRSTEP] = 10
##     umfpack.report_control()
        sol = self.umfpack( um.UMFPACK_A, mtx, rhs, autoTranspose = True )
##     umfpack.report_info()
##    tt = time.clock()
##    vecDX2 = umfpack( um.UMFPACK_At, mtxA, vecR )
##    print "solve = ", time.clock() - tt
##    print nla.norm( vecDX1 - vecDX ), nla.norm( vecDX2 - vecDX )
    
        return sol

    ##
    # c: 10.10.2007, r: 06.03.2008
    def _presolve( self ):
        if hasattr( self, 'presolve' ):
            return self.presolve
        else:
            try:
                return self.conf.presolve
            except:
                return False

##
# c: 22.02.2008
class ScipyIterative( LinearSolver ):
    name = 'ls.scipy_iterative'

    ##
    # c: 22.02.2008, r: 23.06.2008
    def __init__( self, conf, **kwargs ):
        if scipy.version.version < '0.7.0.dev3861':
            import scipy.linalg as la
        else:
            if scipy.version.version < '0.7.0.dev3998':
                import scipy.splinalg.isolve as la
            else:
                import scipy.sparse.linalg.isolve as la
        solver = getattr( la, conf.method )
        LinearSolver.__init__( self, conf, solver = solver, **kwargs )
        
    ##
    # c: 22.02.2008, r: 22.02.2008
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        conf = getDefault( conf, self.conf )
        mtx = getDefault( mtx, self.mtx )
        status = getDefault( status, self.status )

        sol, info = self.solver( mtx, rhs, tol = conf.epsA, maxiter = conf.iMax )
        
        return sol

##
# c: 02.05.2008, r: 02.05.2008
class PyAMGSolver( LinearSolver ):
    name = 'ls.pyamg'

    ##
    # c: 02.05.2008, r: 02.05.2008
    def __init__( self, conf, **kwargs ):
        try:
            import pyamg
        except ImportError:
            output( 'cannot import pyamg!' )
            raise

        try:
            solver = getattr( pyamg, conf.method )
        except AttributeError:
            output( 'pyamg.%s does not exist!' % conf.method )
            output( 'using pyamg.smoothed_aggregation_solver instead' )
            solver = pyamg.smoothed_aggregation_solver

        LinearSolver.__init__( self, conf,
                               solver = solver, mg = None,
                               **kwargs )

        if hasattr( self, 'mtx' ):
            if self.mtx is not None:
                self.mg = self.solver( self.mtx )
        
    ##
    # c: 02.05.2008, r: 02.05.2008
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        conf = getDefault( conf, self.conf )
        mtx = getDefault( mtx, self.mtx )
        status = getDefault( status, self.status )

        if (self.mg is None) or (mtx is not self.mtx):
            self.mg = self.solver( mtx )

        sol = self.mg.solve( rhs, tol = conf.epsA )
        
        return sol
