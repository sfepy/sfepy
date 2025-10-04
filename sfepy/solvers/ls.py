import hashlib

import numpy as nm
import warnings

import scipy.sparse as sps

warnings.simplefilter('ignore', sps.SparseEfficiencyWarning)

from sfepy.base.base import output, get_default, assert_, try_imports
from sfepy.base.timing import Timer
from sfepy.solvers.solvers import LinearSolver

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
        timer = Timer(start=True)

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

        elapsed = timer.stop()
        if status is not None:
            status['time'] = elapsed
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
        timer = Timer(start=True)

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

        elapsed = timer.stop()
        if status is not None:
            status['time'] = elapsed
            status['n_iter'] = self.ksp.getIterationNumber()

        return result

    return _petsc_call


class ScipyDirect(LinearSolver):
    """
    Direct sparse solver from SciPy.
    """
    name = 'ls.scipy_direct'

    _parameters = [
        ('method', "{'auto', 'umfpack', 'superlu'}", 'auto', False,
         'The actual solver to use.'),
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
        ('use_mtx_digest', 'bool', True, False,
         """If True, determine automatically a reused matrix using its
            SHA1 digest. If False, .clear() has to be called
            manually whenever the matrix changes - expert use only!"""),
    ]

    def __init__(self, conf, method=None, **kwargs):
        LinearSolver.__init__(self, conf, solve=None, **kwargs)
        self.sls = None
        if method is None:
            method = self.conf.method

        aux = try_imports(['import scipy.sparse.linalg as sls',
                           'import scipy.sparse.linalg.dsolve as sls'],
                          'cannot import scipy sparse direct solvers!')
        if 'sls' in aux:
            self.sls = aux['sls']
        else:
            raise ValueError('SuperLU not available!')

        if method in ['auto', 'umfpack']:
            aux = try_imports(['import scikits.umfpack as um'])

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

        self.clear()

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        if not conf.use_presolve:
            self.clear()

        if conf.use_presolve:
            self.presolve(mtx, use_mtx_digest=conf.use_mtx_digest)

            # Matrix is already prefactorized.
            return self.solve(rhs)

        else:
            return self.sls.spsolve(mtx, rhs)

    def clear(self):
        if self.solve is not None:
            del self.solve

        self.solve = None

    def presolve(self, mtx, use_mtx_digest=True):
        if use_mtx_digest:
            is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest)

        else:
            is_new, mtx_digest = False, None

        if is_new or (self.solve is None):
            self.solve = self.sls.factorized(mtx.tocsc())
            self.mtx_digest = mtx_digest


class ScipySuperLU(ScipyDirect):
    """
    SuperLU - direct sparse solver from SciPy.
    """
    name = 'ls.scipy_superlu'

    _parameters = [
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
        ('use_mtx_digest', 'bool', True, False,
         """If True, determine automatically a reused matrix using its
            SHA1 digest. If False, .clear() has to be called
            manually whenever the matrix changes - expert use only!"""),
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
        ('use_mtx_digest', 'bool', True, False,
         """If True, determine automatically a reused matrix using its
            SHA1 digest. If False, .clear() has to be called
            manually whenever the matrix changes - expert use only!"""),
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
        import scipy.sparse.linalg as la

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
        if conf.method == 'gmres':
            from packaging import version
            import scipy as sp

            if version.parse(sp.__version__) >= version.parse('1.4.0'):
                solver_kwargs.update({'callback_type' : 'legacy'})

        try:
            sol, info = self.solver(mtx, rhs, x0=x0, atol=eps_a, rtol=eps_r,
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

class PyPardisoSolver(LinearSolver):
    """
    PyPardiso direct solver.

    PyPardiso (https://github.com/haasad/PyPardiso) is a python package to
    solve large sparse linear systems of equations with the Intel oneAPI Math
    Kernel Library PARDISO solver, a shared-memory multiprocessing parallel
    direct sparse solver. PyPardiso is not a python interface to the PARDISO
    Solver Project (https://panua.ch/pardiso).
    """
    name = 'ls.pypardiso'

    _parameters = [
        ('use_presolve', 'bool', False, False,
         """If True, pre-factorize the matrix. It is not needed for performance
            here, as it just calls pypardiso.spsolve() on zeros."""),
    ]

    def __init__(self, conf, method=None, **kwargs):
        aux = try_imports(['import pypardiso as pp'],
                          'cannot import PyPardiso!')
        LinearSolver.__init__(self, conf, pp=aux['pp'], **kwargs)

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        return self.pp.spsolve(mtx, rhs)

    def presolve(self, mtx, use_mtx_digest=True):
        # PyPardiso does its own digest.
        self.pp.spsolve(mtx, nm.zeros(mtx.shape[0], dtype=mtx.dtype))

class PyAMGSolver(LinearSolver):
    """
    Interface to PyAMG solvers.

    The `method` parameter can be one of: 'smoothed_aggregation_solver',
    'ruge_stuben_solver'. The `accel` parameter specifies the Krylov
    solver name, that is used as an accelerator for the multigrid solver.
    """
    name = 'ls.pyamg'

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
                       for key, val in solver_kwargs.items()
                       if key.startswith('method:')}
            self.mg = self.solver(mtx, **_kwargs)
            self.mtx_digest = mtx_digest

        _kwargs = {key[6:] : val
                   for key, val in solver_kwargs.items()
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
        ('create_matrix', 'callable', None, False,
         """User-defined function returning the linear system matrix and
            optionally a precoditioning matrix in a PETSc format. It is called
            as create_matrix(self, mtx, context, comm), where self is the
            solver instance, mtx is the matrix passed to self.__call__(),
            context is a user-supplied context and comm the PETSc communicator.
            """),
        ('block_size', 'int', None, False,
         'The number of degrees of freedom per node.'),
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
        for key, val in petsc.KSP.ConvergedReason.__dict__.items():
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
        for key, rng in field_ranges.items():
            if isinstance(rng, slice):
                rng = rng.start, rng.stop

            size = rng[1] - rng[0]
            field_is = self.petsc.IS().createStride(size, first=rng[0], step=1,
                                                    comm=comm)
            self.fields.append((key, field_is))

    def create_ksp(self, options=None, comm=None):
        optDB = self.petsc.Options()

        if self.conf.sub_precond != 'none':
            optDB['sub_pc_type'] = self.conf.sub_precond

        if options is not None:
            for key, val in options.items():
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

    def init_ksp(self, conf=None, comm=None):
        if conf is None:
            conf = self.conf

        solver_kwargs = self.build_solver_kwargs(conf)
        self.ksp = self.create_ksp(options=solver_kwargs, comm=comm)
        return self.ksp

    @petsc_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, comm=None, context=None,
                 **kwargs):
        eps_a = get_default(eps_a, conf.eps_a)
        eps_r = get_default(eps_r, conf.eps_r)
        i_max = get_default(i_max, conf.i_max)
        eps_d = conf.eps_d

        is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest,
                                            force_reuse=conf.force_reuse)
        ppmtx = None
        if (not is_new) and (self.pmtx is not None):
            pmtx = self.pmtx
        else:
            is_new = True

            if conf.create_matrix is None:
                pmtx = self.create_petsc_matrix(mtx, comm=comm)

            else:
                pmtx = conf.create_matrix(self, mtx, context, comm)
                if not isinstance(pmtx, self.petsc.Mat):
                    pmtx, ppmtx = pmtx # matrix and preconditioning matrix.

                # Convert to Mat if necessary.
                pmtx = self.create_petsc_matrix(pmtx, comm=comm)
                if ppmtx is not None:
                    ppmtx = self.create_petsc_matrix(ppmtx, comm=comm)

            if conf.block_size is not None:
                pmtx.setBlockSize(conf.block_size)

            self.mtx_digest = mtx_digest
            self.pmtx = pmtx

        if self.ksp is not None:
            ksp = self.ksp
            if is_new:
                ksp.setOperators(A=pmtx, P=ppmtx)

        else:
            ksp = self.init_ksp(conf=conf, comm=comm)
            ksp.setOperators(A=pmtx, P=ppmtx)

        ksp.setTolerances(atol=eps_a, rtol=eps_r, divtol=eps_d,
                          max_it=i_max)

        if conf.setup_precond is not None:
            ksp.pc.setPythonContext(conf.setup_precond(mtx, context))

        ksp.setFromOptions()

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
               % (ksp.getType(), ksp.getPC().getType(), conf.sub_precond,
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

    _parameters = [
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
        ('use_mtx_digest', 'bool', True, False,
         """If True, determine automatically a reused matrix using its
            SHA1 digest. If False, .clear() has to be called
            manually whenever the matrix changes - expert use only!"""),
        ('memory_relaxation', 'int', 20, False,
         'The percentage increase in the estimated working space.'),
    ]

    @staticmethod
    def coo_is_symmetric(mtx, tol=1e-9):
        """Symmetry check of the sparse matrix."""
        row, col, data = mtx.row, mtx.col, mtx.data

        out_of_diag = (row != col)
        row, col, data = row[out_of_diag], col[out_of_diag], data[out_of_diag]

        idxs_u = nm.where(col > row)[0]
        idxs_l = nm.where(col < row)[0]
        if idxs_l.shape[0] != idxs_u.shape[0]:
            return False

        iu = nm.lexsort((row[idxs_u], col[idxs_u]))
        il = nm.lexsort((col[idxs_l], row[idxs_l]))

        err = nm.abs(data[idxs_u[iu]] - data[idxs_l[il]])
        if nm.any(err/nm.abs(data[idxs_u].max()) >= tol):
            return False

        return True

    def __init__(self, conf, **kwargs):
        aux = try_imports(['import mumpspy as mumps', 'import mumps'],
                          'cannot import MUMPS!')

        LinearSolver.__init__(self, conf, mumps=aux['mumps'], mumps_ls=None,
                              mumps_presolved=False, **kwargs)
        self.clear()

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        if not conf.use_presolve:
            self.clear()

        self.presolve(mtx, use_mtx_digest=conf.use_mtx_digest)

        return self.mumps_ls.solve(rhs)

    def clear(self):
        if self.mumps_ls is not None:
            del self.mumps_ls

        self.mumps_ls = None

    def presolve(self, mtx, use_mtx_digest=True, factorize=True):
        if use_mtx_digest:
            is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest)

        else:
            is_new, mtx_digest = False, None

        if is_new or (self.mumps_ls is None):
            if not isinstance(mtx, sps.coo_matrix):
                mtx = mtx.tocoo()

            is_sym = self.coo_is_symmetric(mtx)

            if self.mumps_ls is None:
                if self.mumps.__name__ == 'mumpspy':
                    system = 'complex' if mtx.dtype.name.startswith('complex')\
                        else 'real'
                    mem_relax = self.conf.memory_relaxation
                    self.mumps_ls = self.mumps.MumpsSolver(system=system,
                                                           is_sym=is_sym,
                                                           mem_relax=mem_relax)
                    if self.conf.verbose:
                        self.mumps_ls.set_verbose()
                else:
                    self.mumps_ls = self.mumps.Context(self.conf.verbose)

            if self.mumps.__name__ == 'mumpspy':
                self.mumps_ls.set_mtx(mtx, factorize=factorize)
            else:
                self.mumps_ls.set_matrix(mtx, symmetric=is_sym)
                if factorize:
                    self.mumps_ls.factor()

            self.mtx_digest = mtx_digest

    def __del__(self):
        self.clear()


class CholeskySolver(ScipyDirect):
    """
    Interface to scikit-sparse.cholesky solver.
    """
    name = 'ls.cholesky'

    _parameters = [
        ('use_presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
        ('use_mtx_digest', 'bool', True, False,
         """If True, determine automatically a reused matrix using its
            SHA1 digest. If False, .clear() has to be called
            manually whenever the matrix changes - expert use only!"""),
    ]

    def __init__(self, conf, **kwargs):
        LinearSolver.__init__(self, conf, solve=None, **kwargs)
        self.sls = None

        aux = try_imports(['from sksparse.cholmod import cholesky'],
                          'cannot import cholesky sparse solver!')
        if 'cholesky' in aux:
            self.sls = aux['cholesky']
        else:
            raise ValueError('cholesky not available!')

        self.clear()

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        if not conf.use_presolve:
            self.clear()

        self.presolve(mtx, use_mtx_digest=conf.use_mtx_digest)

        return self.solve(rhs)

    def presolve(self, mtx, use_mtx_digest=True):
        if use_mtx_digest:
            is_new, mtx_digest = _is_new_matrix(mtx, self.mtx_digest)

        else:
            is_new, mtx_digest = False, None

        if is_new or (self.solve is None):
            self.solve = self.sls(mtx.tocsc())
            self.mtx_digest = mtx_digest


class SchurMumps(MUMPSSolver):
    r"""
    Mumps Schur complement solver.
    """
    name = 'ls.schur_mumps'

    _parameters = MUMPSSolver._parameters + [
        ('schur_variables', 'list', None, True,
         'The list of Schur variables.'),
    ]

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        if not conf.use_presolve:
            self.clear()

        schur_list = []
        for schur_var in conf.schur_variables:
            slc = self.context.equations.variables.adi.indx[schur_var]
            schur_list.append(nm.arange(slc.start, slc.stop, slc.step, dtype='i'))

        schur_list = nm.hstack(schur_list)

        self.presolve(mtx, use_mtx_digest=conf.use_mtx_digest, factorize=False)

        if self.mumps.__name__ == 'mumpspy':
            # shur_list indexing starts from 1!
            return self.mumps_ls.schur_solve(schur_list + 1, rhs)
        else:
            self.mumps_ls.schur(schur_list)
            return self.mumps_ls.solve_schur(rhs)


class MultiProblem(ScipyDirect):
    r"""
    Conjugate multiple problems.

    Allows to define conjugate multiple problems.
    """
    name = 'ls.cm_pb'

    _parameters = ScipyDirect._parameters + [
        ('others', 'list', None, True,
         'The list of auxiliary problem definition files.'),
        ('coupling_variables', 'list', None, True,
         'The list of coupling variables.'),
    ]

    def __init__(self, conf, context=None, **kwargs):
        ScipyDirect.__init__(self, conf, context=context, **kwargs)

    def init_subproblems(self, conf, **kwargs):
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
        for ii in self.adi_indx.values():
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
            pbi.setup_default_output(conf=confi)
            pbi.set_output_dir(problem.output_dir)
            pbi.equations.set_data(None, ignore_unknown=True)
            pbi.time_update()
            pbi.update_materials()
            sti = pbi.get_initial_state()
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
        for varname, pbs in self.cvars_to_pb.items():
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
        for ii in self.adi_indx.values():
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

        for jk, jv in adi_indxi.items():
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
            x0i = sti0.get_state(pbi.active_only)
            evi = pbi.get_evaluator()
            mtxi = evi.eval_tangent_matrix(x0i, mtx=pbi.mtx_a)
            rhsi = evi.eval_residual(x0i)
            mtxs.append(mtxi)

            adi_indxi = pbi.equations.variables.adi.indx
            for ik, iv in adi_indxi.items():
                if ik in self.cvars_to_pb:
                    if not(self.cvars_to_pb[ik][0] == kk):
                        continue

                giv = self.adi_indx[ik]
                for jk, jv in adi_indxi.items():
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
        for varname, pbs in self.cvars_to_pb.items():
            idx = pbs[1]
            pbi = self.subpb[idx][0]
            mtxi = mtxs[idx]
            gjv = self.adi_indx[varname]
            jv = pbi.equations.variables.adi.indx[varname]
            adi_indxi = pbi.equations.variables.adi.indx
            for ik, iv in adi_indxi.items():
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
            for ii in adi_indxi.values():
                max_indx = nm.max([max_indx, ii.stop])

            resi = nm.zeros((max_indx,), dtype=res0.dtype)
            for ik, iv in adi_indxi.items():
                giv = self.adi_indx[ik]
                if ik in self.cvars_to_pb:
                    if pbi is self.subpb[self.cvars_to_pb[ik][1]][0]:
                        giv = self.cvars_to_pb_map[ik] + giv.start

                resi[iv] = res0[giv]

            if sti0 is not None:
                sti = sti0.copy()
                sti.set_state(-resi, pbi.active_only)
                pbi.save_state(pbi.get_output_name(), sti)
                self.subpb[kk][-1] = sti

            res.append(resi)

        return res[-1]


class RMMSolver(LinearSolver):
    """
    Special solver for explicit transient elastodynamics.

    The solver uses the reciprocal mass matrix algorithm [1]_, [2]_ to directly
    construct a sparse inverse mass matrix. Instead of solving a linear system,
    calling the solver simply performs a sparse matrix multiplication.

    Limitations:

    - Assumes that the density is constant in time.
    - Uses the direct EBC application, i.e., no EBC projection matrix.

    .. [1] González, J.A., Kolman, R., Cho, S.S., Felippa, C.A., Park, K.C.,
           2018. Inverse mass matrix via the method of localized Lagrange
           multipliers. International Journal for Numerical Methods in
           Engineering 113, 277–295. https://doi.org/10.1002/nme.5613

    .. [2] González, J.A., Kopačka, J., Kolman, R., Cho, S.S., Park, K.C.,
           2019. Inverse mass matrix for isogeometric explicit transient
           analysis via the method of localized Lagrange multipliers.
           International Journal for Numerical Methods in Engineering 117,
           939–966. https://doi.org/10.1002/nme.5986
    """
    name = 'ls.rmm'

    _parameters = [
        ('rmm_term', 'str', None, True,
         """The RMM term definition, see
         :class:`MassTerm <sfepy.terms.terms_mass.MassTerm>`."""),
        ('debug', 'bool', False, False,
         'If True, run in debug mode.'),
    ]

    def __init__(self, conf, context=None, **kwargs):
        LinearSolver.__init__(self, conf, context=context, mtx_im=None, a0=None,
                              **kwargs)

    def init_rmm(self, mtx):
        from sfepy.discrete.evaluate import eval_equations, apply_ebc_to_matrix

        problem = self.context
        equations, variables = problem.create_evaluable(
            self.conf.rmm_term, preserve_caches=True,
            copy_materials=False, mode='weak',
            active_only=problem.active_only,
        )
        vu = next(variables.iter_state())

        mtx_a = eval_equations(equations, variables, preserve_caches=True,
                               mode='weak', dw_mode='matrix', term_mode='DPM',
                               active_only=problem.active_only)
        if not problem.active_only:
            apply_ebc_to_matrix(mtx_a, vu.eq_map.eq_ebc,
                                (vu.eq_map.master, vu.eq_map.slave))
        mtx_a.eliminate_zeros()
        mtx_ia = mtx_a.copy()
        mtx_ia.setdiag(1.0 / mtx_a.diagonal())

        mtx_c = eval_equations(equations, variables, preserve_caches=True,
                               mode='weak', dw_mode='matrix', term_mode='RMM',
                               active_only=problem.active_only)
        mtx_c.eliminate_zeros()

        mtx_im = mtx_ia @ (mtx_c @ mtx_ia)
        if not problem.active_only:
            apply_ebc_to_matrix(mtx_im, vu.eq_map.eq_ebc,
                                (vu.eq_map.master, vu.eq_map.slave))

        if self.conf.debug:
            mtx_m = eval_equations(
                equations, variables, preserve_caches=True,
                mode='weak', dw_mode='matrix', term_mode=None,
                active_only=problem.active_only,
            )

            if not problem.active_only:
                mtx_r = vu.eq_map.get_operator()
                mtx_imr = mtx_r.T @ mtx_im @ mtx_r
                mtx_mr = mtx_r.T @ mtx_m @ mtx_r
                mtx_mor = mtx_r.T @ mtx @ mtx_r

            else:
                mtx_imr = mtx_im
                mtx_mr = mtx_m
                mtx_mor = mtx

            dim = problem.domain.shape.dim
            output('total mass check: AMM:', mtx_mr.sum() / dim,
                   'RMM:', nm.linalg.inv(mtx_imr.toarray()).sum() / dim,
                   'M:', mtx_mor.sum() / dim)

        return mtx_im

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):

        if self.mtx_im is None:
            self.mtx_im = self.init_rmm(mtx)

        sol = self.mtx_im @ rhs
        if self.a0 is not None:
            # To make RMMSolver work with the standard Newton solver, a0 has to
            # be set to the previous acceleration and M term has to be
            # nullified (use dw_zero). This option is not used in the current
            # implementation of CentralDifferenceTS.
            sol += self.a0

        return sol
