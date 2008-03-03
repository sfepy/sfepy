from sfe.base.base import *
from sfe.solvers.solvers import EigenvalueSolver

##
# c: 03.03.2008, r: 03.03.2008
class SymeigEigenvalueSolver( EigenvalueSolver ):
    name = 'eig.symeig'
    
    ##
    # c: 03.03.2008, r: 03.03.2008
    def __init__( self, conf, **kwargs ):
        EigenvalueSolver.__init__( self, conf, **kwargs )
        try:
            import symeig
            self.symeig = symeig.symeig
        except:
            self.symeig = None
            output( 'cannot import symeig, required by %s solver' % self.name )
            raise

    ##
    # c: 03.03.2008, r: 03.03.2008
    def __call__( self, mtxA, mtxB = None, nEigs = None,
                  eigenvectors = None, status = None, conf = None ):
        conf = getDefault( conf, self.conf )
        mtxA = getDefault( mtxA, self.mtxA )
        mtxB = getDefault( mtxB, self.mtxB )
        nEigs = getDefault( nEigs, self.nEigs )
        eigenvectors = getDefault( eigenvectors, self.eigenvectors )
        status = getDefault( status, self.status )

        if nEigs is None:
            rng = None
        else:
            rng = (1, nEigs)

        tt = time.clock()
        mtxA = mtxA.toarray()
        if mtxB is not None:
            mtxB = mtxB.toarray()
        out = self.symeig( mtxA, mtxB, range = rng, eigenvectors = eigenvectors )
        if status is not None:
            status['time'] = time.clock() - tt

        return out
##
# c: 03.03.2008, r: 03.03.2008
class ScipyEigenvalueSolver( EigenvalueSolver ):
    name = 'eig.scipy'
    
    ##
    # c: 03.03.2008, r: 03.03.2008
    def __init__( self, conf, **kwargs ):
        EigenvalueSolver.__init__( self, conf, **kwargs )

    ##
    # c: 03.03.2008, r: 03.03.2008
    def __call__( self, mtxA, mtxB = None, nEigs = None,
                  eigenvectors = None, status = None, conf = None ):
        conf = getDefault( conf, self.conf )
        mtxA = getDefault( mtxA, self.mtxA )
        mtxB = getDefault( mtxB, self.mtxB )
        nEigs = getDefault( nEigs, self.nEigs )
        eigenvectors = getDefault( eigenvectors, self.eigenvectors )
        status = getDefault( status, self.status )

        tt = time.clock()
        if nEigs is None:
            mtxA = mtxA.toarray()
            if mtxB is not None:
                mtxB = mtxB.toarray()
            out = nla.eig( mtxA, mtxB, right = eigenvectors )
            eigs = out[0]
            ii = nm.argsort( eigs )
            if eigenvectors:
                mtxEV = out[1][:,ii]
                out = (eigs[ii], mtxEV)
            else:
                out = (eigs,)
        else:
            out = sc.splinalg.eigen_symmetric( mtxA, k = nEigs, M = mtxB )

        if status is not None:
            status['time'] = time.clock() - tt

        return out

##
# c: 03.03.2008, r: 03.03.2008
class PysparseEigenvalueSolver( EigenvalueSolver ):
    name = 'eig.pysparse'

    def _convert_mat(mtx):
        from pysparse import spmatrix
        A = spmatrix.ll_mat(*mtx.shape)
        for i in xrange( mtx.indptr.shape[0] - 1 ):
            ii = slice( mtx.indptr[i], mtx.indptr[i+1] )
            nInRow = ii.stop - ii.start
            A.update_add_at( mtx.data[ii], [i] * nInRow, mtx.indices[ii] )
        return A
    _convert_mat = staticmethod( _convert_mat )

    ##
    # c: 03.03.2008, r: 03.03.2008
    def __init__( self, conf, **kwargs ):
        EigenvalueSolver.__init__( self, conf, **kwargs )

    ##
    # c: 03.03.2008, r: 03.03.2008
    def __call__( self, mtxA, mtxB = None, nEigs = None,
                  eigenvectors = None, status = None, conf = None ):
        from pysparse import jdsym, itsolvers, precon
        conf = getDefault( conf, self.conf )
        mtxA = getDefault( mtxA, self.mtxA )
        mtxB = getDefault( mtxB, self.mtxB )
        nEigs = getDefault( nEigs, self.nEigs )
        eigenvectors = getDefault( eigenvectors, self.eigenvectors )
        status = getDefault( status, self.status )

        output( "loading..." )
        A = self._convert_mat( mtxA )
        output( "...done" )
        if mtxB is not None:
            M = self._convert_mat( mtxB )

        output( "solving..." )
        tt = time.clock()
        Atau=A.copy()
        Atau.shift(-conf.tau,M)
        K=precon.jacobi(Atau)
        A=A.to_sss();
        if mtxB is not None:
            M=M.to_sss();

        method = getattr( itsolvers, conf.method )
        kconv, lmbd, Q, it, it_in = jdsym.jdsym( A, M, K, nEigs, conf.tau,
                                                 conf.epsA, conf.iMax, 
                                                 method,
                                                 clvl = conf.verbosity,
                                                 strategy = conf.strategy )

        output( "number of converged eigenvalues:", kconv )
        ttt = time.clock() - tt
        output( '...done in %.2f s' % ttt )

        if status is not None:
            status['time'] = ttt
            status['Q'] = Q
            status['it'] = it
            status['it_in'] = it_in

        return lmbd, Q
