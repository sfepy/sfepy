from sfepy.base.base import *
from sfepy.solvers.solvers import EigenvalueSolver

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
    # c: 03.03.2008, r: 08.04.2008
    def __call__( self, mtx_a, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None, conf = None ):
        conf = get_default( conf, self.conf )
        mtx_a = get_default( mtx_a, self.mtx_a )
        mtx_b = get_default( mtx_b, self.mtx_b )
        n_eigs = get_default( n_eigs, self.n_eigs )
        eigenvectors = get_default( eigenvectors, self.eigenvectors )
        status = get_default( status, self.status )

        if n_eigs is None:
            rng = None
        else:
            rng = (1, n_eigs)

        tt = time.clock()
        mtx_a, mtx_b = self._to_array( mtx_a, mtx_b )
        out = self.symeig( mtx_a, mtx_b, range = rng, eigenvectors = eigenvectors )
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
    # c: 03.03.2008, r: 08.04.2008
    def __call__( self, mtx_a, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None, conf = None ):
        conf = get_default( conf, self.conf )
        mtx_a = get_default( mtx_a, self.mtx_a )
        mtx_b = get_default( mtx_b, self.mtx_b )
        n_eigs = get_default( n_eigs, self.n_eigs )
        eigenvectors = get_default( eigenvectors, self.eigenvectors )
        status = get_default( status, self.status )

        tt = time.clock()
        if n_eigs is None:
            mtx_a, mtx_b = self._to_array( mtx_a, mtx_b )
            out = nla.eig( mtx_a, mtx_b, right = eigenvectors )
            if eigenvectors:
                eigs = out[0]
            else:
                eigs = out
            ii = nm.argsort( eigs )
            if eigenvectors:
                mtx_ev = out[1][:,ii]
                out = (eigs[ii], mtx_ev)
            else:
                out = (eigs,)
        else:
            out = sc.splinalg.eigen_symmetric( mtx_a, k = n_eigs, M = mtx_b )

        if status is not None:
            status['time'] = time.clock() - tt

        return out

##
# c: 08.04..2008, r: 08.04..2008
class ScipySGEigenvalueSolver( ScipyEigenvalueSolver ):
    name = 'eig.sgscipy'
    
    ##
    # c: 08.04..2008, r: 08.04..2008
    def __call__( self, mtx_a, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None, conf = None ):
        """eigenvectors arg ignored, computes them always"""
        import scipy.lib.lapack as ll
        conf = get_default( conf, self.conf )
        mtx_a = get_default( mtx_a, self.mtx_a )
        mtx_b = get_default( mtx_b, self.mtx_b )
        n_eigs = get_default( n_eigs, self.n_eigs )
        eigenvectors = get_default( eigenvectors, self.eigenvectors )
        status = get_default( status, self.status )

        tt = time.clock()
        if n_eigs is None:
            mtx_a, mtx_b = self._to_array( mtx_a, mtx_b )
            if nm.iscomplexobj( mtx_a ):
                if mtx_b is None:
                    fun = ll.get_lapack_funcs( ['heev'], arrays = (mtx_a,) )[0]
                else:
                    fun = ll.get_lapack_funcs( ['hegv'], arrays = (mtx_a,) )[0]
            else:
                if mtx_b is None:
                    fun = ll.get_lapack_funcs( ['syev'], arrays = (mtx_a,) )[0]
                else:
                    fun = ll.get_lapack_funcs( ['sygv'], arrays = (mtx_a,) )[0]
    ##         print fun
            if mtx_b is None:
                out = fun( mtx_a )
            else:
                out = fun( mtx_a, mtx_b )

            if not eigenvectors:
                out = out[0]
            else:
                out = out[:-1]
            
        else:
            out = ScipyEigenvalueSolver.__call__( self, mtx_a, mtx_b, n_eigs,
                  eigenvectors, status = None )

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
            n_in_row = ii.stop - ii.start
            A.update_add_at( mtx.data[ii], [i] * n_in_row, mtx.indices[ii] )
        return A
    _convert_mat = staticmethod( _convert_mat )

    ##
    # c: 03.03.2008, r: 03.03.2008
    def __init__( self, conf, **kwargs ):
        EigenvalueSolver.__init__( self, conf, **kwargs )

    ##
    # c: 03.03.2008, r: 03.03.2008
    def __call__( self, mtx_a, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None, conf = None ):
        from pysparse import jdsym, itsolvers, precon
        conf = get_default( conf, self.conf )
        mtx_a = get_default( mtx_a, self.mtx_a )
        mtx_b = get_default( mtx_b, self.mtx_b )
        n_eigs = get_default( n_eigs, self.n_eigs )
        eigenvectors = get_default( eigenvectors, self.eigenvectors )
        status = get_default( status, self.status )

        output( "loading..." )
        A = self._convert_mat( mtx_a )
        output( "...done" )
        if mtx_b is not None:
            M = self._convert_mat( mtx_b )

        output( "solving..." )
        tt = time.clock()
        Atau=A.copy()
        Atau.shift(-conf.tau,M)
        K=precon.jacobi(Atau)
        A=A.to_sss();
        if mtx_b is not None:
            M=M.to_sss();

        method = getattr( itsolvers, conf.method )
        kconv, lmbd, Q, it, it_in = jdsym.jdsym( A, M, K, n_eigs, conf.tau,
                                                 conf.eps_a, conf.i_max, 
                                                 method,
                                                 clvl = conf.verbosity,
                                                 strategy = conf.strategy )

        output( "number of converged eigenvalues:", kconv )
        ttt = time.clock() - tt
        output( '...done in %.2f s' % ttt )

        if status is not None:
            status['time'] = ttt
            status['q'] = Q
            status['it'] = it
            status['it_in'] = it_in

        return lmbd, Q
