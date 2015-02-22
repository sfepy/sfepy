import time

import numpy as nm

from sfepy.base.base import output, get_default, try_imports, Struct
from sfepy.solvers.solvers import SolverMeta, Solver, EigenvalueSolver

def eig(mtx_a, mtx_b=None, n_eigs=None, eigenvectors=True,
        return_time=None, method='eig.scipy', **ckwargs):
    """
    Utility function that constructs an eigenvalue solver given by
    `method`, calls it and returns solution.
    """
    kwargs = {'name' : 'aux', 'kind' : method}
    kwargs.update(ckwargs)
    conf = Struct(**kwargs)
    solver = Solver.any_from_conf(conf)

    status = {}
    out = solver(mtx_a, mtx_b, n_eigs, eigenvectors, status)
    if return_time is not None:
        return_time[0] = status['time']

    return out

def standard_call(call):
    """
    Decorator handling argument preparation and timing for eigensolvers.
    """
    def _standard_call(self, mtx_a, mtx_b=None, n_eigs=None,
                       eigenvectors=None, status=None, conf=None, **kwargs):
        tt = time.clock()

        conf = get_default(conf, self.conf)
        mtx_a = get_default(mtx_a, self.mtx_a)
        mtx_b = get_default(mtx_b, self.mtx_b)
        n_eigs = get_default(n_eigs, self.n_eigs)
        eigenvectors = get_default(eigenvectors, self.eigenvectors)
        status = get_default(status, self.status)

        result = call(self, mtx_a, mtx_b, n_eigs, eigenvectors, status, conf,
                      **kwargs)

        ttt = time.clock() - tt
        if status is not None:
            status['time'] = ttt

        return result

    return _standard_call

class ScipyEigenvalueSolver(EigenvalueSolver):
    """
    SciPy-based solver for both dense and sparse problems (if `n_eigs`
    is given).
    """
    name = 'eig.scipy'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', "{'eig', 'eigh'}", 'eig', False,
         """The method for solving general or symmetric eigenvalue problems:
            for dense problems :func:`eig()` or :func:`eigh()` are used, for
            sparse problems (if `n_eigs` is given) :func:`eigs()` or
            :func:`eigsh()` are used."""),
        ('*', '*', None, False,
         'Additional parameters supported by the method.'),
    ]

    def __init__(self, conf, **kwargs):
        EigenvalueSolver.__init__(self, conf, **kwargs)

        import scipy.linalg as sla
        self.sla = sla

        aux = try_imports(['import scipy.sparse.linalg as ssla'],
                          'cannot import scipy sparse eigenvalue solvers!')
        self.ssla = aux['ssla']

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None):
        kwargs = self.build_solver_kwargs(conf)

        if n_eigs is None:
            mtx_a, mtx_b = self._to_array(mtx_a, mtx_b)
            eig = self.sla.eig if conf.method == 'eig' else self.sla.eigh
            out = eig(mtx_a, mtx_b, right=eigenvectors, **kwargs)

        else:
            eig = self.ssla.eigs if conf.method == 'eig' else self.ssla.eigsh
            out = eig(mtx_a, M=mtx_b, k=n_eigs, which='SM',
                      return_eigenvectors=eigenvectors, **kwargs)

        if eigenvectors:
            eigs = out[0]

        else:
            eigs = out

        ii = nm.argsort(eigs)

        if eigenvectors:
            mtx_ev = out[1][:, ii]
            out = (eigs[ii], mtx_ev)

        else:
            out = eigs[ii]

        return out

class ScipySGEigenvalueSolver(EigenvalueSolver):
    """
    SciPy-based solver for dense symmetric problems.
    """
    name = 'eig.sgscipy'

    __metaclass__ = SolverMeta

    def __init__(self, conf, **kwargs):
        EigenvalueSolver.__init__(self, conf, **kwargs)

        import scipy.lib.lapack as llapack
        self.llapack = llapack

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None):
        """
        Notes
        -----
        Eigenvectors argument is ignored, as they are computed always.
        """
        ll = self.llapack

        mtx_a, mtx_b = self._to_array(mtx_a, mtx_b)
        if nm.iscomplexobj(mtx_a):
            if mtx_b is None:
                fun = ll.get_lapack_funcs(['heev'], arrays=(mtx_a,))[0]
            else:
                fun = ll.get_lapack_funcs(['hegv'], arrays=(mtx_a,))[0]
        else:
            if mtx_b is None:
                fun = ll.get_lapack_funcs(['syev'], arrays=(mtx_a,))[0]
            else:
                fun = ll.get_lapack_funcs(['sygv'], arrays=(mtx_a,))[0]

        if mtx_b is None:
            out = fun(mtx_a)
        else:
            out = fun(mtx_a, mtx_b)

        if not eigenvectors:
            if n_eigs is None:
                out = out[0]
            else:
                out = out[0][:n_eigs]
        else:
            if n_eigs is None:
                out = out[:-1]
            else:
                out = (out[0][:n_eigs], out[1][:, :n_eigs])

        return out

class LOBPCGEigenvalueSolver(EigenvalueSolver):
    """
    SciPy-based LOBPCG solver for sparse symmetric problems.
    """
    name = 'eig.scipy_lobpcg'

    __metaclass__ = SolverMeta

    _parameters = [
        ('i_max', 'int', 20, False,
         'The maximum number of iterations.'),
        ('eps_a', 'float', None, False,
         'The absolute tolerance for the convergence.'),
        ('largest', 'bool', True, False,
         'If True, solve for the largest eigenvalues, otherwise the smallest.'),
        ('precond', '{dense matrix, sparse matrix, LinearOperator}',
         None, False,
         'The preconditioner.'),
    ]

    def __init__(self, conf, **kwargs):
        EigenvalueSolver.__init__(self, conf, **kwargs)

        from scipy.sparse.linalg.eigen import lobpcg
        self.lobpcg = lobpcg

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None):

        if n_eigs is None:
            n_eigs = mtx_a.shape[0]
        else:
            n_eigs = min(n_eigs, mtx_a.shape[0])

        x = nm.zeros((mtx_a.shape[0], n_eigs), dtype=nm.float64)
        x[:n_eigs] = nm.eye(n_eigs, dtype=nm.float64)

        out = self.lobpcg(mtx_a, x, mtx_b,
                          M=conf.precond,
                          tol=conf.eps_a, maxiter=conf.i_max,
                          largest=conf.largest,
                          verbosityLevel=conf.verbose)

        if not eigenvectors:
            out = out[0]

        return out

class PysparseEigenvalueSolver(EigenvalueSolver):
    """
    Pysparse-based eigenvalue solver for sparse symmetric problems.
    """
    name = 'eig.pysparse'

    __metaclass__ = SolverMeta

    _parameters = [
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_a', 'float', 1e-5, False,
         'The absolute tolerance for the convergence.'),
        ('tau', 'float', 0.0, False,
         'The target value.'),
        ('method', "{'cgs', 'qmrs'}", 'qmrs', False,
         'The linear iterative solver supported by ``pysparse``.'),
        ('verbosity', 'int', 0, False,
         'The ``pysparse`` verbosity level.'),
        ('strategy', '{0, 1}', 1, False,
         """The shifting and sorting strategy of JDSYM: strategy=0 enables the
            default JDSYM algorithm, strategy=1 enables JDSYM to avoid
            convergence to eigenvalues smaller than `tau`."""),
    ]

    @staticmethod
    def _convert_mat(mtx):
        from pysparse import spmatrix
        A = spmatrix.ll_mat(*mtx.shape)
        for i in xrange(mtx.indptr.shape[0] - 1):
            ii = slice(mtx.indptr[i], mtx.indptr[i+1])
            n_in_row = ii.stop - ii.start
            A.update_add_at(mtx.data[ii], [i] * n_in_row, mtx.indices[ii])
        return A

    def __init__(self, conf, **kwargs):
        EigenvalueSolver.__init__(self, conf, **kwargs)

        imp = try_imports(['from pysparse import jdsym, itsolvers, precon',
                           'from pysparse.eigen import jdsym;'
                           ' from pysparse import itsolvers, precon'],
                          'cannot import pysparse eigensolvers!')
        self.imp = imp

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None,
                 eigenvectors=None, status=None, conf=None):
        jdsym = self.imp['jdsym']
        itsolvers = self.imp['itsolvers']
        precon = self.imp['precon']

        output('solving...', verbose=conf.verbose)

        A = self._convert_mat(mtx_a)
        Atau = A.copy()

        if mtx_b is not None:
            M = self._convert_mat(mtx_b)
            Atau.shift(-conf.tau, M)

        K = precon.jacobi(Atau)
        A = A.to_sss()

        if mtx_b is not None:
            M = M.to_sss()

        else:
            M = None

        method = getattr(itsolvers, conf.method)
        kconv, lmbd, Q, it, it_in = jdsym.jdsym(A, M, K, n_eigs, conf.tau,
                                                conf.eps_a, conf.i_max,
                                                method,
                                                clvl=conf.verbosity,
                                                strategy=conf.strategy)

        output('number of converged eigenvalues:', kconv, verbose=conf.verbose)

        output('...done', verbose=conf.verbose)

        if status is not None:
            status['q'] = Q
            status['it'] = it
            status['it_in'] = it_in

        ii = nm.argsort(lmbd)

        if eigenvectors:
            out = (lmbd[ii], Q[:, ii])

        else:
            out = lmbd[ii]

        return out
