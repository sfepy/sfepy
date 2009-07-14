from sfepy.base.base import *
from sfepy.solvers.solvers import EigenvalueSolver

##
# c: 25.09.2007, r: 08.04.2008
def eig( mtx_a, mtx_b = None, num = None, eigenvectors = True,
         return_time = None, method = 'eig.scipy', **ckwargs ):

    kwargs = {'name' : 'aux', 'kind' : method}
    kwargs.update( ckwargs )
    conf = Struct( **kwargs )
    solver = Solver.any_from_conf( conf )

    status = {}
    out = solver( mtx_a, mtx_b, num, eigenvectors, status )
    if return_time is not None:
        return_time[0] = status['time']
        
    return out

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
            try:
                from scipy.splinalg import eigen_symmetric
            except ImportError:
                eigen_symmetric = None

            try:
                from scipy.sparse.linalg.eigen.arpack import eigen_symmetric
            except ImportError:
                eigen_symmetric = None

            if eigen_symmetric is None:
                raise ImportError('cannot import eigen_symmetric!')
                
            out = eigen_symmetric( mtx_a, k = n_eigs, M = mtx_b )

        if status is not None:
            status['time'] = time.clock() - tt

        return out

##
# c: 08.04..2008, r: 08.04..2008
class ScipySGEigenvalueSolver( ScipyEigenvalueSolver ):
    """Solver for symmetric problems."""
    name = 'eig.sgscipy'
    
    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
       solver_20 = {
            'name' : 'eigen2',
            'kind' : 'eig.sgscipy',

            'force_n_eigs' : True,
        }
        """
        get = conf.get_default_attr

        force_n_eigs = get( 'force_n_eigs', False )

        common = EigenvalueSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

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
        if (n_eigs is None) or (conf.force_n_eigs):
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
                if n_eigs is None:
                    out = out[0]
                else:
                    out = out[0][:n_eigs]
            else:
                if n_eigs is None:
                    out = out[:-1]
                else:
                    out = (out[0][:n_eigs], out[1][:,:n_eigs])
                    
        else:
            out = ScipyEigenvalueSolver.__call__( self, mtx_a, mtx_b, n_eigs,
                  eigenvectors, status = None )

        if status is not None:
            status['time'] = time.clock() - tt

        return out


class LOBPCGEigenvalueSolver( EigenvalueSolver ):
    name = 'eig.scipy_lobpcg'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_2 = {
            'name' : 'lobpcg',
            'kind' : 'eig.scipy_lobpcg',

            'i_max' : 20,
            'n_eigs' : 5,
            'eps_a' : None,
            'largest' : True,
            'precond' : None,
            'verbosity' : 0,
        }
        """
        get = conf.get_default_attr

        i_max = get( 'i_max', 20 )
        n_eigs = get( 'n_eigs', None )
        eps_a = get( 'eps_a', None )
        largest = get( 'largest', True )
        precond = get( 'precond', None )
        verbosity = get( 'verbosity', 0 )

        common = EigenvalueSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

    def __init__( self, conf, **kwargs ):
        EigenvalueSolver.__init__( self, conf, **kwargs )

        from scipy.sparse.linalg.eigen import lobpcg
        self.lobpcg = lobpcg

    def __call__( self, mtx_a, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None, conf = None ):
        conf = get_default( conf, self.conf )
        mtx_a = get_default( mtx_a, self.mtx_a )
        mtx_b = get_default( mtx_b, self.mtx_b )
        n_eigs = get_default( n_eigs, self.n_eigs )
        eigenvectors = get_default( eigenvectors, self.eigenvectors )
        status = get_default( status, self.status )

        if n_eigs is None:
            n_eigs = mtx_a.shape[0]
        else:
            n_eigs = min(n_eigs, mtx_a.shape[0])

        tt = time.clock()

        x = nm.zeros( (mtx_a.shape[0], n_eigs), dtype = nm.float64 )
        x[:n_eigs] = nm.eye( n_eigs, dtype = nm.float64 )
        out = self.lobpcg( mtx_a, x, mtx_b,
                           M = conf.precond,
                           tol = conf.eps_a, maxiter = conf.i_max,
                           largest = conf.largest,
                           verbosityLevel = conf.verbosity )

        if status is not None:
            status['time'] = time.clock() - tt

        if not eigenvectors:
            out = out[0]

        return out

##
# c: 03.03.2008, r: 03.03.2008
class PysparseEigenvalueSolver( EigenvalueSolver ):
    name = 'eig.pysparse'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_2 = {
            'name' : 'eigen1',
            'kind' : 'eig.pysparse',

            'i_max' : 150,
            'eps_a' : 1e-5,
            'tau' : -10.0,
            'method' : 'qmrs',
            'verbosity' : 0,
            'strategy' : 1,
        }
        """
        get = conf.get_default_attr

        i_max = get( 'i_max', 100 )
        eps_a = get( 'eps_a', 1e-5 )
        tau = get( 'tau', 0.0 )
        method = get( 'method', 'qmrs' )
        verbosity = get( 'verbosity', 0 )
        strategy = get( 'strategy', 1 )

        common = EigenvalueSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

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
