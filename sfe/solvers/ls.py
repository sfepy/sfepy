from sfe.base.base import *
from sfe.solvers.solvers import LinearSolver

import scipy
if scipy.version.version == "0.6.0":
    import scipy.linsolve.umfpack as um
else:
    try:
        import scipy.splinalg.dsolve.umfpack as um
    except ImportError:
        import scikits.umfpack as um

um.configure( assumeSortedIndices = True )

##
# 10.10.2007, c
class Umfpack( LinearSolver ):
    name = 'ls.umfpack'

    ##
    # 10.10.2007, c
    def __init__( self, conf, **kwargs ):
        LinearSolver.__init__( self, conf, **kwargs )
        self.umfpack = um.UmfpackContext()

        if self._presolve() and hasattr( self, 'mtx' ):
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

##    tt = time.clock()
##    vecDX = sp.solve( mtxA, vecR )
##    print "solve = ", time.clock() - tt

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
    # 10.10.2007, c
    def _presolve( self ):
        try:
            return self.conf.presolve
        except:
            return False
