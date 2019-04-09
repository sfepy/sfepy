from __future__ import absolute_import
import time
import hashlib

import numpy as nm
import warnings

import scipy.sparse as sps
import six
from six.moves import range

warnings.simplefilter('ignore', sps.SparseEfficiencyWarning)

from sfepy.base.base import output, get_default, assert_, try_imports
from sfepy.solvers.solvers import SolverMeta, LinearSolver

def solve(mtx, rhs, solver_class=None, solver_conf=None):
    """
    Solve the linear system with the matrix `mtx` and the right-hand side
    `rhs`.

    Convenience wrapper around the linear solver classes below.
    """
    solver_class = get_default(solver_class, ScipyDirect)
    solver_conf = get_default(solver_conf, {})

    solver = solver_class(solver_conf, mtx=mtx)
    solution = solver(rhs)

    return solution

def _get_cs_matrix_hash(mtx, chunk_size=100000):
    def _gen_array_chunks(arr):
        ii = 0
        while len(arr[ii:]):
            yield arr[ii:ii+chunk_size].tobytes()
            ii += chunk_size

    sha1 = hashlib.sha1()
    for chunk in _gen_array_chunks(mtx.indptr):
        sha1.update(chunk)
    for chunk in _gen_array_chunks(mtx.indices):
        sha1.update(chunk)
    for chunk in _gen_array_chunks(mtx.data):
        sha1.update(chunk)

    digest = sha1.hexdigest()
    return digest

def _is_new_matrix(mtx, mtx_digest, force_reuse=False):
    if not isinstance(mtx, sps.csr_matrix):
        return True, mtx_digest

    if force_reuse:
        return False, mtx_digest

    id0, digest0 = mtx_digest
    id1 = id(mtx)
    digest1 = _get_cs_matrix_hash(mtx)
    if (id1 == id0) and (digest1 == digest0):
        return False, (id1, digest1)

    return True, (id1, digest1)

def standard_call(call):
    """
    Decorator handling argument preparation and timing for linear solvers.
    """
    def _standard_call(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                       i_max=None, mtx=None, status=None, context=None,
                       **kwargs):
        tt = time.clock()

        conf = get_default(conf, self.conf)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)
        context = get_default(context, self.context)

        assert_(mtx.shape[0] == mtx.shape[1] == rhs.shape[0])
        if x0 is not None:
            assert_(x0.shape[0] == rhs.shape[0])

        result = call(self, rhs, x0, conf, eps_a, eps_r, i_max, mtx, status,
                      context=context, **kwargs)
        if isinstance(result, tuple):
            result, n_iter = result

        else:
            n_iter = -1 # Number of iterations is undefined/unavailable.

        ttt = time.clock() - tt
        if status is not None:
            status['time'] = ttt
            status['n_iter'] = n_iter

        return result

    return _standard_call

def petsc_call(call):
    """
    Decorator handling argument preparation and timing for PETSc-based linear
    solvers.
    """
    def _petsc_call(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                    i_max=None, mtx=None, status=None, comm=None,
                    context=None, **kwargs):
        tt = time.clock()

        conf = get_default(conf, self.conf)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)
        context = get_default(context, self.context)
        comm = get_default(comm, self.comm)

        mshape = mtx.size if isinstance(mtx, self.petsc.Mat) else mtx.shape
        rshape = [rhs.size] if isinstance(rhs, self.petsc.Vec) else rhs.shape

        assert_(mshape[0] == mshape[1] == rshape[0])
        if x0 is not None:
            xshape = [x0.size] if isinstance(x0, self.petsc.Vec) else x0.shape
            assert_(xshape[0] == rshape[0])

        result = call(self, rhs, x0, conf, eps_a, eps_r, i_max, mtx, status,
                      comm, context=context, **kwargs)

        ttt = time.clock() - tt
        if status is not None:
            status['time'] = ttt
            status['n_iter'] = self.ksp.getIterationNumber()

        return result

    return _petsc_call


class ScipyDirect(LinearSolver):
    """
    Direct sparse solver from SciPy.
    """
    name = 'ls.scipy_direct'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', "{'auto', 'umfpack', 'superlu'}", 'auto', False,
         'The actual solver to use.'),
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
    ]

    def __init__(self, conf, method=None, **kwargs):
        LinearSolver.__init__(self, conf, solve=None, **kwargs)
        um = self.sls = None
        if method is None:
            method = self.conf.method

        aux = try_imports(['import scipy.linsolve as sls',
                           'import scipy.splinalg.dsolve as sls',
                           'import scipy.sparse.linalg.dsolve as sls'],
                          'cannot import scipy sparse direct solvers!')
        if 'sls' in aux:
            self.sls = aux['sls']
        else:
            raise ValueError('SuperLU not available!')

        if method in ['auto', 'umfpack']:
            aux = try_imports([
                'import scipy.linsolve.umfpack as um',
                'import scipy.splinalg.dsolve.umfpack as um',
                'import scipy.sparse.linalg.dsolve.umfpack as um',
                'import scikits.umfpack as um'])

            is_umfpack = True if 'um' in aux\
                and hasattr(aux['um'], 'UMFPACK_OK') else False
            if method == 'umfpack' and not is_umfpack:
                raise ValueError('UMFPACK not available!')
        elif method == 'superlu':
            is_umfpack = False
        else:
            raise ValueError('uknown solution method! (%s)' % method)

        if is_umfpack:
            self.sls.use_solver(useUmfpack=True,
                                assumeSortedIndices=True)
        else:
            self.sls.use_solver(useUmfpack=False)

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):

        if conf.use_presolve:
            self.presolve(mtx)

        if self.solve is not None:
            # Matrix is already prefactorized.
            return self.solve(rhs)
        else:
            return self.sls.spsolve(mtx, rhs)

    def presolve(self, mtx):
        is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest)
        if is_new:
            self.solve = self.sls.factorized(mtx)
            self.mtx_digest = mtx_digest


class ScipySuperLU(ScipyDirect):
    """
    SuperLU - direct sparse solver from SciPy.
    """
    name = 'ls.scipy_superlu'

    _parameters = [
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
    ]

    def __init__(self, conf, **kwargs):
        ScipyDirect.__init__(self, conf, method='superlu', **kwargs)


class ScipyUmfpack(ScipyDirect):
    """
    UMFPACK - direct sparse solver from SciPy.
    """
    name = 'ls.scipy_umfpack'

    _parameters = [
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
    ]

    def __init__(self, conf, **kwargs):
        ScipyDirect.__init__(self, conf, method='umfpack', **kwargs)


class ScipyIterative(LinearSolver):
    """
    Interface to SciPy iterative solvers.

    The `eps_r` tolerance is both absolute and relative - the solvers
    stop when either the relative or the absolute residual is below it.
    """
    name = 'ls.scipy_iterative'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'cg', False,
         'The actual solver to use.'),
        ('setup_precond', 'callable', lambda mtx, context: None, False,
         """User-supplied function for the preconditioner initialization/setup.
            It is called as setup_precond(mtx, context), where mtx is the
            matrix, context is a user-supplied context, and should return one
            of {sparse matrix, dense matrix, LinearOperator}.
         """),
        ('callback', 'callable', None, False,
         """User-supplied function to call after each iteration. It is called
            as callback(xk), where xk is the current solution vector, except
            the gmres method, where the argument is the residual.
         """),
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_a', 'float', 1e-8, False,
         'The absolute tolerance for the residual.'),
        ('eps_r', 'float', 1e-8, False,
         'The relative tolerance for the residual.'),
        ('*', '*', None, False,
         'Additional parameters supported by the method.'),
    ]

    # All iterative solvers in scipy.sparse.linalg pass a solution vector into
    # a callback except those below, that take a residual vector.
    _callbacks_res = ['gmres']

    def __init__(self, conf, context=None, **kwargs):
        import scipy.sparse.linalg.isolve as la

        LinearSolver.__init__(self, conf, context=context, **kwargs)

        try:
            solver = getattr(la, self.conf.method)
        except AttributeError:
            output('scipy solver %s does not exist!' % self.conf.method)
            output('using cg instead')
            solver = la.cg
        self.solver = solver
        self.converged_reasons = {
            0 : 'successful exit',
            1 : 'number of iterations',
            -1 : 'illegal input or breakdown',
        }

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, context=None, **kwargs):
        solver_kwargs = self.build_solver_kwargs(conf)

        eps_a = get_default(eps_a, self.conf.eps_a)
        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)

        setup_precond = get_default(kwargs.get('setup_precond', None),
                                    self.conf.setup_precond)
        callback = get_default(kwargs.get('callback', lambda sol: None),
                               self.conf.callback)

        self.iter = 0
        def iter_callback(sol):
            self.iter += 1
            msg = '%s: iteration %d' % (self.conf.name, self.iter)
            if conf.verbose > 2:
                if conf.method not in self._callbacks_res:
                    res = mtx * sol - rhs

                else:
                    res = sol

                rnorm = nm.linalg.norm(res)
                msg += ': |Ax-b| = %e' % rnorm
            output(msg, verbose=conf.verbose > 1)

            # Call an optional user-defined callback.
            callback(sol)

        precond = setup_precond(mtx, context)

        if conf.method == 'qmr':
            prec_args = {'M1' : precond, 'M2' : precond}

        else:
            prec_args = {'M' : precond}

        solver_kwargs.update(prec_args)

        try:
            sol, info = self.solver(mtx, rhs, x0=x0, atol=eps_a, tol=eps_r,
                                    maxiter=i_max, callback=iter_callback,
                                    **solver_kwargs)
        except TypeError:
            sol, info = self.solver(mtx, rhs, x0=x0, tol=eps_r,
                                    maxiter=i_max, callback=iter_callback,
                                    **solver_kwargs)

        output('%s: %s convergence: %s (%s, %d iterations)'
               % (self.conf.name, self.conf.method,
                  info, self.converged_reasons[nm.sign(info)], self.iter),
               verbose=conf.verbose)

        return sol, self.iter

class PyAMGSolver(LinearSolver):
    """
    Interface to PyAMG solvers.

    The `method` parameter can be one of: 'smoothed_aggregation_solver',
    'ruge_stuben_solver'. The `accel` parameter specifies the Krylov
    solver name, that is used as an accelerator for the multigrid solver.
    """
    name = 'ls.pyamg'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'smoothed_aggregation_solver', False,
         'The actual solver to use.'),
        ('accel', 'str', None, False,
         'The accelerator.'),
        ('callback', 'callable', None, False,
         """User-supplied function to call after each iteration. It is called
            as callback(xk), where xk is the current solution vector, except
            the gmres accelerator, where the argument is the residual norm.
         """),
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_r', 'float', 1e-8, False,
         'The relative tolerance for the residual.'),
        ('force_reuse', 'bool', False, False,
         """If True, skip the check whether the MG solver object corresponds
            to the `mtx` argument: it is always reused."""),
        ('*', '*', None, False,
         """Additional parameters supported by the method. Use the 'method:'
            prefix for arguments of the method construction function
            (e.g. 'method:max_levels' : 5), and the 'solve:' prefix for
            the subsequent solver call."""),
    ]

    # All iterative solvers in pyamg.krylov pass a solution vector into
    # a callback except those below, that take a residual vector norm.
    _callbacks_res = ['gmres']

    def __init__(self, conf, **kwargs):
        try:
            import pyamg
        except ImportError:
            msg =  'cannot import pyamg!'
            raise ImportError(msg)

        LinearSolver.__init__(self, conf, mg=None, **kwargs)

        try:
            solver = getattr(pyamg, self.conf.method)
        except AttributeError:
            output('pyamg.%s does not exist!' % self.conf.method)
            output('using pyamg.smoothed_aggregation_solver instead')
            solver = pyamg.smoothed_aggregation_solver
        self.solver = solver

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        solver_kwargs = self.build_solver_kwargs(conf)

        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)

        callback = get_default(kwargs.get('callback', lambda sol: None),
                               self.conf.callback)

        self.iter = 0
        def iter_callback(sol):
            self.iter += 1
            msg = '%s: iteration %d' % (self.conf.name, self.iter)
            if conf.verbose > 2:
                if conf.accel not in self._callbacks_res:
                    res = mtx * sol - rhs

                else:
                    res = sol

                rnorm = nm.linalg.norm(res)
                msg += ': |Ax-b| = %e' % rnorm
            output(msg, verbose=conf.verbose > 1)

            # Call an optional user-defined callback.
            callback(sol)

        is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest,
                                            force_reuse=conf.force_reuse)
        if is_new or (self.mg is None):
            _kwargs = {key[7:] : val
                       for key, val in six.iteritems(solver_kwargs)
                       if key.startswith('method:')}
            self.mg = self.solver(mtx, **_kwargs)
            self.mtx_digest = mtx_digest

        _kwargs = {key[6:] : val
                   for key, val in six.iteritems(solver_kwargs)
                   if key.startswith('solve:')}
        sol = self.mg.solve(rhs, x0=x0, accel=conf.accel, tol=eps_r,
                            maxiter=i_max, callback=iter_callback,
                            **_kwargs)

        return sol, self.iter

class PyAMGKrylovSolver(LinearSolver):
    """
    Interface to PyAMG Krylov solvers.
    """
    name = 'ls.pyamg_krylov'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'cg', False,
         'The actual solver to use.'),
        ('setup_precond', 'callable', lambda mtx, context: None, False,
         """User-supplied function for the preconditioner initialization/setup.
            It is called as setup_precond(mtx, context), where mtx is the
            matrix, context is a user-supplied context, and should return one
            of {sparse matrix, dense matrix, LinearOperator}.
         """),
        ('callback', 'callable', None, False,
         """User-supplied function to call after each iteration. It is called
            as callback(xk), where xk is the current solution vector, except
            the gmres method, where the argument is the residual norm.
         """),
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_r', 'float', 1e-8, False,
         'The relative tolerance for the residual.'),
        ('*', '*', None, False,
         'Additional parameters supported by the method.'),
    ]

    # All iterative solvers in pyamg.krylov pass a solution vector into
    # a callback except those below, that take a residual vector norm.
    _callbacks_res = ['gmres']

    def __init__(self, conf, context=None, **kwargs):
        try:
            import pyamg.krylov as krylov
        except ImportError:
            msg =  'cannot import pyamg.krylov!'
            raise ImportError(msg)

        LinearSolver.__init__(self, conf, mg=None,
                              context=context, **kwargs)

        try:
            solver = getattr(krylov, self.conf.method)
        except AttributeError:
            output('pyamg.krylov.%s does not exist!' % self.conf.method)
            raise

        self.solver = solver
        self.converged_reasons = {
            0 : 'successful exit',
            1 : 'number of iterations',
            -1 : 'illegal input or breakdown',
        }

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, context=None, **kwargs):
        solver_kwargs = self.build_solver_kwargs(conf)

        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)

        setup_precond = get_default(kwargs.get('setup_precond', None),
                                    self.conf.setup_precond)
        callback = get_default(kwargs.get('callback', lambda sol: None),
                               self.conf.callback)

        self.iter = 0
        def iter_callback(sol):
            self.iter += 1
            msg = '%s: iteration %d' % (self.conf.name, self.iter)
            if conf.verbose > 2:
                if conf.method not in self._callbacks_res:
                    res = mtx * sol - rhs

                else:
                    res = sol

                rnorm = nm.linalg.norm(res)
                msg += ': |Ax-b| = %e' % rnorm
            output(msg, verbose=conf.verbose > 1)

            # Call an optional user-defined callback.
            callback(sol)

        precond = setup_precond(mtx, context)

        sol, info = self.solver(mtx, rhs, x0=x0, tol=eps_r, maxiter=i_max,
                                M=precond, callback=iter_callback,
                                **solver_kwargs)

        output('%s: %s convergence: %s (%s, %d iterations)'
               % (self.conf.name, self.conf.method,
                  info, self.converged_reasons[nm.sign(info)], self.iter),
               verbose=conf.verbose)

        return sol, self.iter

class PETScKrylovSolver(LinearSolver):
    """
    PETSc Krylov subspace solver.

    The solver supports parallel use with a given MPI communicator (see `comm`
    argument of :func:`PETScKrylovSolver.__init__()`) and allows passing in
    PETSc matrices and vectors. Returns a (global) PETSc solution vector
    instead of a (local) numpy array, when given a PETSc right-hand side
    vector.

    The solver and preconditioner types are set upon the solver object
    creation. Tolerances can be overridden when called by passing a `conf`
    object.

    Convergence is reached when `rnorm < max(eps_r * rnorm_0, eps_a)`,
    where, in PETSc, `rnorm` is by default the norm of *preconditioned*
    residual.
    """
    name = 'ls.petsc'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'cg', False,
         'The actual solver to use.'),
        ('setup_precond', 'callable', None, False,
         """User-supplied function for the preconditioner initialization/setup.
            It is called as setup_precond(mtx, context), where mtx is the
            matrix, context is a user-supplied context, and should return an
            object with `setUp(self, pc)` and `apply(self, pc, x, y)` methods.

            Has precedence over the `precond`/`sub_precond` parameters.
         """),
        ('precond', 'str', 'icc', False,
         'The preconditioner.'),
        ('sub_precond', 'str', 'none', False,
         'The preconditioner for matrix blocks (in parallel runs).'),
        ('precond_side', "{'left', 'right', 'symmetric', None}", None, False,
         'The preconditioner side.'),
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_a', 'float', 1e-8, False,
         'The absolute tolerance for the residual.'),
        ('eps_r', 'float', 1e-8, False,
         'The relative tolerance for the residual.'),
        ('eps_d', 'float', 1e5, False,
         'The divergence tolerance for the residual.'),
        ('force_reuse', 'bool', False, False,
         """If True, skip the check whether the KSP solver object corresponds
            to the `mtx` argument: it is always reused."""),
        ('*', '*', None, False,
         """Additional parameters supported by the method. Can be used to pass
            all PETSc options supported by :func:`petsc.Options()`."""),
    ]

    _precond_sides = {None : None, 'left' : 0, 'right' : 1, 'symmetric' : 2}

    def __init__(self, conf, comm=None, context=None, **kwargs):
        if comm is None:
            from sfepy.parallel.parallel import init_petsc_args; init_petsc_args

        from petsc4py import PETSc as petsc

        converged_reasons = {}
        for key, val in six.iteritems(petsc.KSP.ConvergedReason.__dict__):
            if isinstance(val, int):
                converged_reasons[val] = key

        LinearSolver.__init__(self, conf, petsc=petsc, comm=comm,
                              converged_reasons=converged_reasons,
                              fields=None, ksp=None, pmtx=None,
                              context=context, **kwargs)

    def set_field_split(self, field_ranges, comm=None):
        """
        Setup local PETSc ranges for fields to be used with 'fieldsplit'
        preconditioner.

        This function must be called before solving the linear system.
        """
        comm = get_default(comm, self.comm)

        self.fields = []
        for key, rng in six.iteritems(field_ranges):
            if isinstance(rng, slice):
                rng = rng.start, rng.stop

            size = rng[1] - rng[0]
            field_is = self.petsc.IS().createStride(size, first=rng[0], step=1,
                                                    comm=comm)
            self.fields.append((key, field_is))

    def create_ksp(self, options=None, comm=None):
        optDB = self.petsc.Options()

        optDB['sub_pc_type'] = self.conf.sub_precond
        if options is not None:
            for key, val in six.iteritems(options):
                optDB[key] = val

        ksp = self.petsc.KSP()
        ksp.create(comm)

        ksp.setType(self.conf.method)
        pc = ksp.getPC()
        if self.conf.setup_precond is None:
            pc.setType(self.conf.precond)

        else:
            pc.setType(pc.Type.PYTHON)
        ksp.setFromOptions()

        if (pc.type == 'fieldsplit'):
            if self.fields is not None:
                pc.setFieldSplitIS(*self.fields)

            else:
                msg = 'PETScKrylovSolver.set_field_split() has to be called!'
                raise ValueError(msg)

        side = self._precond_sides[self.conf.precond_side]
        if side is not None:
            ksp.setPCSide(side)

        return ksp

    def create_petsc_matrix(self, mtx, comm=None):
        if isinstance(mtx, self.petsc.Mat):
            pmtx = mtx

        else:
            mtx = sps.csr_matrix(mtx)

            pmtx = self.petsc.Mat()
            pmtx.createAIJ(mtx.shape, csr=(mtx.indptr, mtx.indices, mtx.data),
                           comm=comm)

        return pmtx

    @petsc_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, comm=None, context=None,
                 **kwargs):
        solver_kwargs = self.build_solver_kwargs(conf)

        eps_a = get_default(eps_a, self.conf.eps_a)
        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)
        eps_d = self.conf.eps_d

        is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest,
                                            force_reuse=conf.force_reuse)
        if (not is_new) and self.ksp is not None:
            ksp = self.ksp
            pmtx = self.pmtx

        else:
            pmtx = self.create_petsc_matrix(mtx, comm=comm)

            ksp = self.create_ksp(options=solver_kwargs, comm=comm)
            ksp.setOperators(pmtx)
            ksp.setTolerances(atol=eps_a, rtol=eps_r, divtol=eps_d,
                              max_it=i_max)

            setup_precond = self.conf.setup_precond
            if setup_precond is not None:
                ksp.pc.setPythonContext(setup_precond(mtx, context))

            ksp.setFromOptions()
            self.mtx_digest = mtx_digest
            self.ksp = ksp
            self.pmtx = pmtx

        if isinstance(rhs, self.petsc.Vec):
            prhs = rhs

        else:
            prhs = pmtx.getVecLeft()
            prhs[...] = rhs

        if x0 is not None:
            if isinstance(x0, self.petsc.Vec):
                psol = x0

            else:
                psol = pmtx.getVecRight()
                psol[...] = x0

            ksp.setInitialGuessNonzero(True)

        else:
            psol = pmtx.getVecRight()

            ksp.setInitialGuessNonzero(False)

        ksp.solve(prhs, psol)
        output('%s(%s, %s/proc) convergence: %s (%s, %d iterations)'
               % (ksp.getType(), ksp.getPC().getType(), self.conf.sub_precond,
                  ksp.reason, self.converged_reasons[ksp.reason],
                  ksp.getIterationNumber()),
               verbose=conf.verbose)

        if isinstance(rhs, self.petsc.Vec):
            sol = psol

        else:
            sol = psol[...].copy()

        return sol


class MUMPSSolver(LinearSolver):
    """
    Interface to MUMPS solver.
    """
    name = 'ls.mumps'

    __metaclass__ = SolverMeta

    _parameters = [
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
    ]

    def __init__(self, conf, **kwargs):
        import sfepy.solvers.ls_mumps as mumps

        self.mumps_ls = None
        mumps.load_mumps_libraries()  # try to load MUMPS libraries

        LinearSolver.__init__(self, conf, mumps=mumps, mumps_ls=None,
                              mumps_presolved=False, **kwargs)

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        if not self.mumps_presolved:
            self.presolve(mtx, presolve_flag=conf.use_presolve)

        out = rhs.copy()
        self.mumps_ls.set_rhs(out)
        self.mumps_ls(3)  # solve

        return out

    def presolve(self, mtx, presolve_flag=False):
        is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest)
        if not isinstance(mtx, sps.coo_matrix):
            mtx = mtx.tocoo()
        if self.mumps_ls is None:
            system = 'complex' if mtx.dtype.name.startswith('complex')\
                else 'real'
            is_sym = self.mumps.coo_is_symmetric(mtx)
            self.mumps_ls = self.mumps.MumpsSolver(system=system,
                                                   is_sym=is_sym)

        if is_new:
            if self.conf.verbose:
                self.mumps_ls.set_verbose()

            self.mumps_ls.set_mtx_centralized(mtx)
            self.mumps_ls(4)  # analyze + factorize
            if presolve_flag:
                self.mumps_presolved = True
            self.mtx_digest = mtx_digest

    def __del__(self):
        if self.mumps_ls is not None:
            del(self.mumps_ls)


class MUMPSParallelSolver(LinearSolver):
    """
    Interface to MUMPS parallel solver.
    """
    name = 'ls.mumps_par'

    __metaclass__ = SolverMeta

    _parameters = []

    def __init__(self, conf, **kwargs):
        import multiprocessing
        import sfepy.solvers.ls_mumps as mumps

        mumps.load_mumps_libraries()  # try to load MUMPS libraries

        LinearSolver.__init__(self, conf, mumps=mumps, mumps_ls=None,
                              number_of_cpu=multiprocessing.cpu_count(),
                              mumps_presolved=False, **kwargs)

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        from mpi4py import MPI
        import sys
        from sfepy import data_dir
        import os.path as op
        from tempfile import gettempdir

        def tmpfile(fname):
            return op.join(gettempdir(), fname)

        if not isinstance(mtx, sps.coo_matrix):
            mtx = mtx.tocoo()

        is_sym = self.mumps.coo_is_symmetric(mtx)
        rr, cc, data = mtx.row + 1, mtx.col + 1, mtx.data
        if is_sym:
            idxs = nm.where(cc >= rr)[0]  # upper triangular matrix
            rr, cc, data = rr[idxs], cc[idxs], data[idxs]

        n = mtx.shape[0]
        nz = rr.shape[0]
        flags = nm.memmap(tmpfile('vals_flags.array'), dtype='int32',
                          mode='w+', shape=(4,))
        flags[0] = n
        flags[1] = 1 if data.dtype.name.startswith('complex') else 0
        flags[2] = int(is_sym)
        flags[3] = int(self.conf.verbose)

        idxs = nm.memmap(tmpfile('idxs.array'), dtype='int32',
                         mode='w+', shape=(2, nz))
        idxs[0, :] = rr
        idxs[1, :] = cc

        dtype = {0: 'float64', 1: 'complex128'}[flags[1]]
        vals_mtx = nm.memmap(tmpfile('vals_mtx.array'), dtype=dtype,
                             mode='w+', shape=(nz,))
        vals_rhs = nm.memmap(tmpfile('vals_rhs.array'), dtype=dtype,
                             mode='w+', shape=(n,))
        vals_mtx[:] = data
        vals_rhs[:] = rhs

        mumps_call = op.join(data_dir, 'sfepy', 'solvers',
                             'ls_mumps_parallel.py')
        comm = MPI.COMM_SELF.Spawn(sys.executable, args=[mumps_call],
                                   maxprocs=self.number_of_cpu)
        comm.Disconnect()

        out = nm.memmap(tmpfile('vals_x.array'), dtype=dtype, mode='r')

        return out


class SchurMumps(MUMPSSolver):
    r"""
    Mumps Schur complement solver.
    """
    name = 'ls.schur_mumps'

    __metaclass__ = SolverMeta

    _parameters = [
        ('schur_variables', 'list', None, True,
         'The list of Schur variables.'),
    ]

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        import scipy.linalg as sla

        if not isinstance(mtx, sps.coo_matrix):
            mtx = mtx.tocoo()

        system = 'complex' if mtx.dtype.name.startswith('complex') else 'real'
        self.mumps_ls = self.mumps.MumpsSolver(system=system)

        if self.conf.verbose:
            self.mumps_ls.set_verbose()

        schur_list = []
        for schur_var in conf.schur_variables:
            slc = self.context.equations.variables.di.indx[schur_var]
            schur_list.append(nm.arange(slc.start, slc.stop, slc.step, dtype='i') + 1)

        self.mumps_ls.set_mtx_centralized(mtx)
        out = rhs.copy()
        self.mumps_ls.set_rhs(out)

        S, y2 = self.mumps_ls.get_schur(nm.hstack(schur_list))
        x2 = sla.solve(S.T, y2)  # solve the dense Schur system using scipy.linalg

        return self.mumps_ls.expand_schur(x2)


class MultiProblem(ScipyDirect):
    r"""
    Conjugate multiple problems.

    Allows to define conjugate multiple problems.
    """
    name = 'ls.cm_pb'

    __metaclass__ = SolverMeta

    _parameters = ScipyDirect._parameters + [
        ('others', 'list', None, True,
         'The list of auxiliary problem definition files.'),
        ('coupling_variables', 'list', None, True,
         'The list of coupling variables.'),
    ]

    def __init__(self, conf, context=None, **kwargs):
        ScipyDirect.__init__(self, conf, context=context, **kwargs)

    def init_subproblems(self, conf, **kwargs):
        from sfepy.discrete.state import State
        from sfepy.discrete import Problem
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from scipy.spatial import cKDTree as KDTree

        # init subproblems
        problem = self.context
        pb_vars = problem.get_variables()
        # get "master" DofInfo and last index
        pb_adi_indx = problem.equations.variables.adi.indx
        self.adi_indx = pb_adi_indx.copy()
        last_indx = -1
        for ii in six.itervalues(self.adi_indx):
            last_indx = nm.max([last_indx, ii.stop])

        # coupling variables
        self.cvars_to_pb = {}
        for jj in conf.coupling_variables:
            self.cvars_to_pb[jj] = [None, None]
            if jj in pb_vars.names:
                if pb_vars[jj].dual_var_name is not None:
                    self.cvars_to_pb[jj][0] = -1

                else:
                    self.cvars_to_pb[jj][1] = -1

        # init subproblems
        self.subpb = []
        required, other = get_standard_keywords()
        master_prefix = output.get_output_prefix()
        for ii, ifname in enumerate(conf.others):
            sub_prefix = master_prefix[:-1] + '-sub%d:' % (ii + 1)
            output.set_output_prefix(sub_prefix)
            kwargs['master_problem'] = problem
            confi = ProblemConf.from_file(ifname, required, other,
                                          define_args=kwargs)
            pbi = Problem.from_conf(confi, init_equations=True)
            sti = State(pbi.equations.variables)
            pbi.equations.set_data(None, ignore_unknown=True)
            pbi.time_update()
            pbi.update_materials()
            sti.apply_ebc()
            pbi_vars = pbi.get_variables()
            output.set_output_prefix(master_prefix)
            self.subpb.append([pbi, sti, None])

            # append "slave" DofInfo
            for jj in pbi_vars.names:
                if not(pbi_vars[jj].is_state()):
                    continue

                didx = pbi.equations.variables.adi.indx[jj]
                ndof = didx.stop - didx.start
                if jj in self.adi_indx:
                    if ndof != \
                      (self.adi_indx[jj].stop - self.adi_indx[jj].start):
                        raise ValueError('DOFs do not match!')

                else:
                    self.adi_indx.update({
                        jj: slice(last_indx, last_indx + ndof, None)})
                    last_indx += ndof

            for jj in conf.coupling_variables:
                if jj in pbi_vars.names:
                    if pbi_vars[jj].dual_var_name is not None:
                        self.cvars_to_pb[jj][0] = ii

                    else:
                        self.cvars_to_pb[jj][1] = ii

        self.subpb.append([problem, None, None])

        self.cvars_to_pb_map = {}
        for varname, pbs in six.iteritems(self.cvars_to_pb):
            # match field nodes
            coors = []
            for ii in pbs:
                pbi = self.subpb[ii][0]
                pbi_vars = pbi.get_variables()
                fcoors = pbi_vars[varname].field.coors
                dc = nm.abs(nm.max(fcoors, axis=0)\
                            - nm.min(fcoors, axis=0))
                ax = nm.where(dc > 1e-9)[0]
                coors.append(fcoors[:,ax])

            if len(coors[0]) != len(coors[1]):
                raise ValueError('number of nodes does not match!')

            kdtree = KDTree(coors[0])
            map_12 = kdtree.query(coors[1])[1]

            pbi1 = self.subpb[pbs[0]][0]
            pbi1_vars = pbi1.get_variables()
            eq_map_1 = pbi1_vars[varname].eq_map

            pbi2 = self.subpb[pbs[1]][0]
            pbi2_vars = pbi2.get_variables()
            eq_map_2 = pbi2_vars[varname].eq_map

            dpn = eq_map_2.dpn
            nnd = map_12.shape[0]

            map_12_nd = nm.zeros((nnd * dpn,), dtype=nm.int32)
            if dpn > 1:
                for ii in range(dpn):
                    map_12_nd[ii::dpn] = map_12 * dpn + ii
            else:
                map_12_nd = map_12

            idx = nm.where(eq_map_2.eq >= 0)[0]
            self.cvars_to_pb_map[varname] = eq_map_1.eq[map_12[idx]]

    def sparse_submat(self, Ad, Ar, Ac, gr, gc, S):
        """
        A[gr,gc] = S
        """

        if type(gr) is slice:
            gr = nm.arange(gr.start, gr.stop)

        if type(gc) is slice:
            gc = nm.arange(gc.start, gc.stop)

        for ii, lrow in enumerate(S):
            m = lrow.indices.shape[0]
            idxrow = nm.ones((m,), dtype=nm.int32) * gr[ii]
            Ar = nm.hstack([Ar, idxrow])
            Ac = nm.hstack([Ac, gc[lrow.indices]])
            Ad = nm.hstack([Ad, lrow.data])

        return Ad, Ar, Ac

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        self.init_subproblems(self.conf, **kwargs)

        max_indx = 0
        hst = nm.hstack
        for ii in six.itervalues(self.adi_indx):
            max_indx = nm.max([max_indx, ii.stop])

        new_rhs = nm.zeros((max_indx,), dtype=rhs.dtype)
        new_rhs[:rhs.shape[0]] = rhs

        # copy "master" matrices
        pbi = self.subpb[-1][0]
        adi_indxi = pbi.equations.variables.adi.indx
        mtxc = mtx.tocsc()
        aux_data = nm.array([], dtype=mtxc.dtype)
        aux_rows = nm.array([], dtype=nm.int32)
        aux_cols = nm.array([], dtype=nm.int32)

        for jk, jv in six.iteritems(adi_indxi):
            if jk in self.cvars_to_pb:
                if not(self.cvars_to_pb[jk][0] == -1):
                    continue

            gjv = self.adi_indx[jk]
            ii = gjv.start
            for jj in nm.arange(jv.start, jv.stop):
                ptr = mtxc.indptr[jj]
                nn = mtxc.indptr[jj + 1] - ptr
                sl = slice(ptr, ptr + nn, None)
                aux_data = hst([aux_data, mtxc.data[sl]])
                aux_rows = hst([aux_rows, mtxc.indices[sl]])
                aux_cols = hst([aux_cols, nm.ones((nn,), dtype=nm.int32) * ii])
                ii += 1

        # copy "slave" (sub)matricies
        mtxs = []
        for kk, (pbi, sti0, _) in enumerate(self.subpb[:-1]):
            x0i = sti0.get_reduced()
            evi = pbi.get_evaluator()
            mtxi = evi.eval_tangent_matrix(x0i, mtx=pbi.mtx_a)
            rhsi = evi.eval_residual(x0i)
            mtxs.append(mtxi)

            adi_indxi = pbi.equations.variables.adi.indx
            for ik, iv in six.iteritems(adi_indxi):
                if ik in self.cvars_to_pb:
                    if not(self.cvars_to_pb[ik][0] == kk):
                        continue

                giv = self.adi_indx[ik]
                for jk, jv in six.iteritems(adi_indxi):
                    gjv = self.adi_indx[jk]
                    if jk in self.cvars_to_pb:
                        if not(self.cvars_to_pb[jk][0] == kk):
                            continue

                    aux_data, aux_rows, aux_cols =\
                        self.sparse_submat(aux_data, aux_rows, aux_cols,
                                           giv, gjv, mtxi[iv, jv])

                new_rhs[giv] = rhsi[iv]

        mtxs.append(mtx)
        # copy "coupling" (sub)matricies
        for varname, pbs in six.iteritems(self.cvars_to_pb):
            idx = pbs[1]
            pbi = self.subpb[idx][0]
            mtxi = mtxs[idx]
            gjv = self.adi_indx[varname]
            jv = pbi.equations.variables.adi.indx[varname]
            adi_indxi = pbi.equations.variables.adi.indx
            for ik, iv in six.iteritems(adi_indxi):
                if ik == varname:
                    continue

                giv = self.adi_indx[ik]
                aux_mtx = mtxi[iv,:].tocsc()
                for ll, jj in enumerate(nm.arange(jv.start, jv.stop)):
                    ptr = aux_mtx.indptr[jj]
                    nn = aux_mtx.indptr[jj + 1] - ptr
                    if nn < 1:
                        continue
                    sl = slice(ptr, ptr + nn, None)
                    aux_data = hst([aux_data, aux_mtx.data[sl]])
                    aux_rows = hst([aux_rows, aux_mtx.indices[sl] + giv.start])
                    jjr = gjv.start + self.cvars_to_pb_map[varname][ll]
                    aux_cols = hst([aux_cols,
                                    nm.ones((nn,), dtype=nm.int32) * jjr])

        # create new matrix
        new_mtx = sps.coo_matrix((aux_data, (aux_rows, aux_cols))).tocsr()

        res0 = ScipyDirect.__call__(self, new_rhs, mtx=new_mtx)

        res = []
        for kk, (pbi, sti0, _) in enumerate(self.subpb):
            adi_indxi = pbi.equations.variables.adi.indx
            max_indx = 0
            for ii in six.itervalues(adi_indxi):
                max_indx = nm.max([max_indx, ii.stop])

            resi = nm.zeros((max_indx,), dtype=res0.dtype)
            for ik, iv in six.iteritems(adi_indxi):
                giv = self.adi_indx[ik]
                if ik in self.cvars_to_pb:
                    if pbi is self.subpb[self.cvars_to_pb[ik][1]][0]:
                        giv = self.cvars_to_pb_map[ik] + giv.start

                resi[iv] = res0[giv]

            if sti0 is not None:
                sti = sti0.copy()
                sti.set_reduced(-resi)
                pbi.setup_default_output()
                pbi.save_state(pbi.get_output_name(), sti)
                self.subpb[kk][-1] = sti

            res.append(resi)

        return res[-1]
