import time

import numpy as nm
import warnings

import scipy.sparse as sps

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

def standard_call(call):
    """
    Decorator handling argument preparation and timing for linear solvers.
    """
    def _standard_call(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                       i_max=None, mtx=None, status=None, **kwargs):
        tt = time.clock()

        conf = get_default(conf, self.conf)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        assert_(mtx.shape[0] == mtx.shape[1] == rhs.shape[0])
        if x0 is not None:
            assert_(x0.shape[0] == rhs.shape[0])

        result = call(self, rhs, x0, conf, eps_a, eps_r, i_max, mtx, status,
                      **kwargs)

        ttt = time.clock() - tt
        if status is not None:
            status['time'] = ttt

        return result

    return _standard_call

class ScipyDirect(LinearSolver):
    """
    Direct sparse solver from SciPy.
    """
    name = 'ls.scipy_direct'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', "{'auto', 'umfpack', 'superlu'}", 'auto', False,
         'The actual solver to use.'),
        ('presolve', 'bool', False, False,
         'If True, pre-factorize the matrix.'),
        ('warn', 'bool', True, False,
         'If True, allow warnings.'),
    ]

    def __init__(self, conf, **kwargs):
        LinearSolver.__init__(self, conf, **kwargs)
        um = self.sls = None

        aux = try_imports(['import scipy.linsolve as sls',
                           'import scipy.splinalg.dsolve as sls',
                           'import scipy.sparse.linalg.dsolve as sls'],
                          'cannot import scipy sparse direct solvers!')
        self.sls = aux['sls']
        aux = try_imports(['import scipy.linsolve.umfpack as um',
                           'import scipy.splinalg.dsolve.umfpack as um',
                           'import scipy.sparse.linalg.dsolve.umfpack as um',
                           'import scikits.umfpack as um'])
        if 'um' in aux:
            um = aux['um']

        if um is not None:
            is_umfpack = hasattr(um, 'UMFPACK_OK')
        else:
            is_umfpack = False

        method = self.conf.method
        if method == 'superlu':
            self.sls.use_solver(useUmfpack=False)
        elif method == 'umfpack':
            if not is_umfpack and self.conf.warn:
                output('umfpack not available, using superlu!')
        elif method != 'auto':
            raise ValueError('uknown solution method! (%s)' % method)

        if method != 'superlu' and is_umfpack:
            self.sls.use_solver(useUmfpack=True,
                                assumeSortedIndices=True)

        self.solve = None
        if self._presolve() and hasattr(self, 'mtx'):
            if self.mtx is not None:
                self.solve = self.sls.factorized(self.mtx)

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):

        if self.solve is not None:
            # Matrix is already prefactorized.
            return self.solve(rhs)
        else:
            return self.sls.spsolve(mtx, rhs)

    def _presolve(self):
        if hasattr(self, 'presolve'):
            return self.presolve
        else:
            return self.conf.presolve

class ScipyIterative(LinearSolver):
    """
    Interface to SciPy iterative solvers.

    The `eps_r` tolerance is both absolute and relative - the solvers
    stop when either the relative or the absolute residual is below it.

    A preconditioner can be anything that the SciPy solvers accept (sparse
    matrix, dense matrix, LinearOperator).
    """
    name = 'ls.scipy_iterative'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'cg', False,
         'The actual solver to use.'),
        ('precond', '{sparse matrix, dense matrix, LinearOperator}',
         None, False,
         'The preconditioner.'),
        ('callback', 'function', None, False,
         """User-supplied function to call after each iteration. It is called
            as callback(xk), where xk is the current solution vector."""),
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_r', 'float', 1e-8, False,
         'The relative or absolute tolerance for the residual.'),
    ]

    def __init__(self, conf, **kwargs):
        import scipy.sparse.linalg.isolve as la

        LinearSolver.__init__(self, conf, **kwargs)

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
                 i_max=None, mtx=None, status=None, **kwargs):

        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)

        precond = get_default(kwargs.get('precond', None), self.conf.precond)
        callback = get_default(kwargs.get('callback', None), self.conf.callback)

        if conf.method == 'qmr':
            prec_args = {'M1' : precond, 'M2' : precond}

        else:
            prec_args = {'M' : precond}

        sol, info = self.solver(mtx, rhs, x0=x0, tol=eps_r, maxiter=i_max,
                                callback=callback, **prec_args)
        output('%s convergence: %s (%s)'
               % (self.conf.method,
                  info, self.converged_reasons[nm.sign(info)]))

        return sol

class PyAMGSolver(LinearSolver):
    """
    Interface to PyAMG solvers.
    """
    name = 'ls.pyamg'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'smoothed_aggregation_solver', False,
         'The actual solver to use.'),
        ('accel', 'str', None, False,
         'The accelerator.'),
        ('i_max', 'int', 100, False,
         'The maximum number of iterations.'),
        ('eps_r', 'float', 1e-8, False,
         'The relative tolerance for the residual.'),
    ]

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

        if hasattr(self, 'mtx'):
            if self.mtx is not None:
                self.mg = self.solver(self.mtx)

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):

        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)

        if (self.mg is None) or (mtx is not self.mtx):
            self.mg = self.solver(mtx)
            self.mtx = mtx

        sol = self.mg.solve(rhs, x0=x0, accel=conf.accel, tol=eps_r,
                            maxiter=i_max)

        return sol

class PETScKrylovSolver(LinearSolver):
    """
    PETSc Krylov subspace solver.

    The solver and preconditioner types are set upon the solver object
    creation. Tolerances can be overriden when called by passing a `conf`
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
        ('precond', 'str', 'icc', False,
         'The preconditioner.'),
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
    ]

    _precond_sides = {None : None, 'left' : 0, 'right' : 1, 'symmetric' : 2}

    def __init__(self, conf, **kwargs):
        try:
            import petsc4py
            petsc4py.init([])
            from petsc4py import PETSc
        except ImportError:
            msg = 'cannot import petsc4py!'
            raise ImportError(msg)

        LinearSolver.__init__(self, conf, petsc=PETSc, pmtx=None, **kwargs)

        ksp = PETSc.KSP().create()

        ksp.setType(self.conf.method)
        ksp.getPC().setType(self.conf.precond)
        side = self._precond_sides[self.conf.precond_side]
        if side is not None:
            ksp.setPCSide(side)
        self.ksp = ksp

        self.converged_reasons = {}
        for key, val in ksp.ConvergedReason.__dict__.iteritems():
            if isinstance(val, int):
                self.converged_reasons[val] = key

    def set_matrix(self, mtx):
        mtx = sps.csr_matrix(mtx)

        pmtx = self.petsc.Mat().createAIJ(mtx.shape,
                                          csr=(mtx.indptr,
                                               mtx.indices,
                                               mtx.data))
        sol, rhs = pmtx.getVecs()
        return pmtx, sol, rhs

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):

        eps_a = get_default(eps_a, self.conf.eps_a)
        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)
        eps_d = self.conf.eps_d

        # There is no use in caching matrix in the solver - always set as new.
        pmtx, psol, prhs = self.set_matrix(mtx)

        ksp = self.ksp
        ksp.setOperators(pmtx)
        ksp.setFromOptions() # PETSc.Options() not used yet...
        ksp.setTolerances(atol=eps_a, rtol=eps_r, divtol=eps_d, max_it=i_max)

        # Set PETSc rhs, solve, get solution from PETSc solution.
        if x0 is not None:
            psol[...] = x0
            ksp.setInitialGuessNonzero(True)
        prhs[...] = rhs
        ksp.solve(prhs, psol)
        sol = psol[...].copy()
        output('%s(%s) convergence: %s (%s)'
               % (self.conf.method, self.conf.precond,
                  ksp.reason, self.converged_reasons[ksp.reason]))

        return sol

class PETScParallelKrylovSolver(PETScKrylovSolver):
    """
    PETSc Krylov subspace solver able to run in parallel by storing the
    system to disk and running a separate script via `mpiexec`.

    The solver and preconditioner types are set upon the solver object
    creation. Tolerances can be overriden when called by passing a `conf`
    object.

    Convergence is reached when `rnorm < max(eps_r * rnorm_0, eps_a)`,
    where, in PETSc, `rnorm` is by default the norm of *preconditioned*
    residual.
    """
    name = 'ls.petsc_parallel'

    __metaclass__ = SolverMeta

    _parameters = PETScKrylovSolver._parameters + [
        ('log_dir', 'str', '.', False,
         'The directory for storing logs.'),
        ('n_proc', 'int', 1, False,
         'The number of processes.'),
        ('sub_precond', 'str', 'icc', False,
         'The preconditioner for matrix blocks.'),
    ]

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        import os, sys, shutil, tempfile
        from sfepy import base_dir
        from sfepy.base.ioutils import ensure_path

        eps_a = get_default(eps_a, self.conf.eps_a)
        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)
        eps_d = self.conf.eps_d

        petsc = self.petsc

        # There is no use in caching matrix in the solver - always set as new.
        pmtx, psol, prhs = self.set_matrix(mtx)

        ksp = self.ksp
        ksp.setOperators(pmtx)
        ksp.setFromOptions() # PETSc.Options() not used yet...
        ksp.setTolerances(atol=eps_a, rtol=eps_r, divtol=eps_d, max_it=i_max)

        output_dir = tempfile.mkdtemp()

        # Set PETSc rhs, solve, get solution from PETSc solution.
        if x0 is not None:
            psol[...] = x0
            sol0_filename = os.path.join(output_dir, 'sol0.dat')

        else:
            sol0_filename = ''

        prhs[...] = rhs

        script_filename = os.path.join(base_dir, 'solvers/petsc_worker.py')

        mtx_filename = os.path.join(output_dir, 'mtx.dat')
        rhs_filename = os.path.join(output_dir, 'rhs.dat')
        sol_filename = os.path.join(output_dir, 'sol.dat')
        status_filename = os.path.join(output_dir, 'status.txt')

        log_filename = os.path.join(self.conf.log_dir, 'sol.log')
        ensure_path(log_filename)

        output('storing system to %s...' % output_dir)
        tt = time.clock()
        view_mtx = petsc.Viewer().createBinary(mtx_filename, mode='w')
        view_rhs = petsc.Viewer().createBinary(rhs_filename, mode='w')
        pmtx.view(view_mtx)
        prhs.view(view_rhs)
        if sol0_filename:
            view_sol0 = petsc.Viewer().createBinary(sol0_filename, mode='w')
            psol.view(view_sol0)
        output('...done in %.2f s' % (time.clock() - tt))

        command = [
            'mpiexec -n %d' % self.conf.n_proc,
            sys.executable, script_filename,
            '-mtx %s' % mtx_filename, '-rhs %s' % rhs_filename,
            '-sol0 %s' % sol0_filename, '-sol %s' % sol_filename,
            '-status %s' % status_filename,
            '-ksp_type %s' % self.conf.method,
            '-pc_type %s' % self.conf.precond,
            '-sub_pc_type %s' % self.conf.sub_precond,
            '-ksp_atol %.3e' % self.conf.eps_a,
            '-ksp_rtol %.3e' % self.conf.eps_r,
            '-ksp_max_it %d' % self.conf.i_max,
            '-ksp_monitor %s' % log_filename,
            '-ksp_view %s' % log_filename,
        ]
        if self.conf.precond_side is not None:
            command.append('-ksp_pc_side %s' % self.conf.precond_side)

        out = os.system(" ".join(command))
        assert_(out == 0)

        output('reading solution...')
        tt = time.clock()
        view_sol = self.petsc.Viewer().createBinary(sol_filename, mode='r')
        psol = petsc.Vec().load(view_sol)

        fd = open(status_filename, 'r')
        line = fd.readline().split()
        reason = int(line[0])
        elapsed = float(line[1])
        fd.close()
        output('...done in %.2f s' % (time.clock() - tt))

        sol = psol[...].copy()
        output('%s(%s, %s/proc) convergence: %s (%s)'
               % (self.conf.method, self.conf.precond, self.conf.sub_precond,
                  reason, self.converged_reasons[reason]))
        output('elapsed: %.2f [s]' % elapsed)

        shutil.rmtree(output_dir)

        return sol

class SchurGeneralized(ScipyDirect):
    r"""
    Generalized Schur complement.

    Defines the matrix blocks and calls user defined function.
    """
    name = 'ls.schur_generalized'

    __metaclass__ = SolverMeta

    _parameters = ScipyDirect._parameters + [
        ('blocks', 'dict', None, True,
         """The description of blocks: ``{block_name1 : [variable_name1, ...],
            ...}``."""),
        ('function', 'callable', None, True,
         'The user defined function.'),
    ]

    def __init__(self, conf, problem, **kwargs):
        from sfepy.discrete.state import State

        ScipyDirect.__init__(self, conf, **kwargs)

        equations = problem.equations
        aux_state = State(equations.variables)

        conf.idxs = {}
        for bk, bv in conf.blocks.iteritems():
            aux_state.fill(0.0)
            for jj in bv:
                idx = equations.variables.di.indx[jj]
                aux_state.vec[idx] = nm.nan

            aux_state.apply_ebc()
            vec0 = aux_state.get_reduced()
            conf.idxs[bk] = nm.where(nm.isnan(vec0))[0]

    @standard_call
    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):

        mtxi= self.orig_conf.idxs
        mtxslc_s = {}
        mtxslc_f = {}
        nn = {}

        for ik, iv in mtxi.iteritems():
            ptr = 0
            nn[ik] = len(iv)
            mtxslc_s[ik] = []
            mtxslc_f[ik] = []
            while ptr < nn[ik]:
                idx0 = iv[ptr:]
                idxrange = nm.arange(idx0[0], idx0[0] + len(idx0))
                aux = nm.where(idx0 == idxrange)[0]
                mtxslc_s[ik].append(slice(ptr + aux[0], ptr + aux[-1] + 1))
                mtxslc_f[ik].append(slice(idx0[aux][0], idx0[aux][-1] + 1))
                ptr += aux[-1] + 1

        mtxs = {}
        rhss = {}
        ress = {}
        get_sub = mtx._get_submatrix
        for ir in mtxi.iterkeys():
            rhss[ir] = nm.zeros((nn[ir],), dtype=nm.float64)
            ress[ir] = nm.zeros((nn[ir],), dtype=nm.float64)
            for jr, idxr in enumerate(mtxslc_f[ir]):
                rhss[ir][mtxslc_s[ir][jr]] = rhs[idxr]

            for ic in mtxi.iterkeys():
                mtxid = '%s%s' % (ir, ic)
                mtxs[mtxid] = nm.zeros((nn[ir], nn[ic]), dtype=nm.float64)
                for jr, idxr in enumerate(mtxslc_f[ir]):
                    for jc, idxc in enumerate(mtxslc_f[ic]):
                        iir = mtxslc_s[ir][jr]
                        iic = mtxslc_s[ic][jc]
                        mtxs[mtxid][iir, iic] = get_sub(idxr, idxc).todense()

        self.orig_conf.function(ress, mtxs, rhss, nn)

        res = nm.zeros_like(rhs)
        for ir in mtxi.iterkeys():
            for jr, idxr in enumerate(mtxslc_f[ir]):
                res[idxr] = ress[ir][mtxslc_s[ir][jr]]

        return res

    def _presolve(self):
        if hasattr(self, 'presolve'):
            return self.presolve
        else:
            return self.conf.presolve

class SchurComplement(SchurGeneralized):
    r"""
    Schur complement.

    Solution of the linear system

    .. math::
       \left[ \begin{array}{cc}
       A & B \\
       C & D \end{array} \right]
       \cdot
       \left[ \begin{array}{c}
       u \\
       v \end{array} \right]
       =
       \left[ \begin{array}{c}
       f \\
       g \end{array} \right]

    is obtained by solving the following equation:

    .. math::
       (D - C A^{-1} B) \cdot v = g - C A^{-1} f

    variable(s) :math:`u` are specified in "eliminate" list,
    variable(s) :math:`v` are specified in "keep" list,

    See: http://en.wikipedia.org/wiki/Schur_complement
    """
    name = 'ls.schur_complement'

    __metaclass__ = SolverMeta

    _parameters = ScipyDirect._parameters + [
        ('eliminate', 'list', None, True,
         'The list of variables to eliminate.'),
        ('keep', 'list', None, True,
         'The list of variables to keep.'),
    ]

    @staticmethod
    def schur_fun(res, mtx, rhs, nn):
        import scipy.sparse as scs
        import scipy.sparse.linalg as sls

        invA = sls.splu(scs.csc_matrix(mtx['11']))
        invAB = nm.zeros_like(mtx['12'])
        for j, b in enumerate(mtx['12'].T):
            invAB[:,j] = invA.solve(b)

        invAf = invA.solve(rhs['1'])

        spC = scs.csc_matrix(mtx['21'])
        k_rhs = rhs['2'] - spC * invAf
        res['2'] = sls.spsolve(scs.csc_matrix(mtx['22'] - spC * invAB), k_rhs)
        res['1'] = invAf - nm.dot(invAB, res['2'])

    def __init__(self, conf, **kwargs):
        get = self.conf.get
        conf.blocks = {'1': get('eliminate', None,
                                'missing "eliminate" in options!'),
                       '2': get('keep', None,
                                'missing "keep" in options!'),}
        conf.function = SchurComplement.schur_fun

        SchurGeneralized.__init__(self, conf, **kwargs)

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

    def __init__(self, conf, problem, **kwargs):
        from sfepy.discrete.state import State
        from sfepy.discrete import Problem
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from scipy.spatial import cKDTree as KDTree

        ScipyDirect.__init__(self, conf, **kwargs)

        # init subproblems
        pb_vars = problem.get_variables()
        # get "master" DofInfo and last index
        pb_adi_indx = problem.equations.variables.adi.indx
        self.adi_indx = pb_adi_indx.copy()
        last_indx = -1
        for ii in self.adi_indx.itervalues():
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
        for varname, pbs in self.cvars_to_pb.iteritems():
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

        max_indx = 0
        hst = nm.hstack
        for ii in self.adi_indx.itervalues():
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

        for jk, jv in adi_indxi.iteritems():
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
            for ik, iv in adi_indxi.iteritems():
                if ik in self.cvars_to_pb:
                    if not(self.cvars_to_pb[ik][0] == kk):
                        continue

                giv = self.adi_indx[ik]
                for jk, jv in adi_indxi.iteritems():
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
        for varname, pbs in self.cvars_to_pb.iteritems():
            idx = pbs[1]
            pbi = self.subpb[idx][0]
            mtxi = mtxs[idx]
            gjv = self.adi_indx[varname]
            jv = pbi.equations.variables.adi.indx[varname]
            adi_indxi = pbi.equations.variables.adi.indx
            for ik, iv in adi_indxi.iteritems():
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
            for ii in adi_indxi.itervalues():
                max_indx = nm.max([max_indx, ii.stop])

            resi = nm.zeros((max_indx,), dtype=res0.dtype)
            for ik, iv in adi_indxi.iteritems():
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

    def _presolve(self):
        if hasattr(self, 'presolve'):
            return self.presolve
        else:
            return self.conf.presolve
