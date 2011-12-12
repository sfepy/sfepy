import numpy as nm
import warnings

import scipy.sparse as sps

warnings.simplefilter('ignore', sps.SparseEfficiencyWarning)

from sfepy.base.base import output, get_default, assert_, Struct
from sfepy.solvers.solvers import make_get_conf, LinearSolver

def try_imports(imports, fail_msg=None):
    for imp in imports:
        try:
            exec imp
            break
        except:
            pass
        else:
            if fail_msg is not None:
                raise ValueError(fail_msg)
    return locals()

class ScipyDirect(LinearSolver):
    name = 'ls.scipy_direct'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are set to default values.

        Example configuration, all items::

            solver_1100 = {
                'name' : 'dls1100',
                'kind' : 'ls.scipy_direct',

                'method' : 'superlu',
                'presolve' : False,
                'warn' : True,
            }
        """
        get = make_get_conf(conf, kwargs)
        common = LinearSolver.process_conf(conf)

        return Struct(method=get('method', 'auto'),
                      presolve=get('presolve', False),
                      warn=get('warn', True),
                      i_max=None, eps_a=None, eps_r=None) + common

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

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        conf = get_default(conf, self.conf)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

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

class Umfpack(ScipyDirect):
    """This class stays for compatability with old input files. Use ScipyDirect
    isntead."""
    name = 'ls.umfpack'

    def __init__(self, conf, **kwargs):
        conf.method = 'umfpack'
        ScipyDirect.__init__(self, conf, **kwargs)

##
# c: 22.02.2008
class ScipyIterative( LinearSolver ):
    """
    Interface to scipy iterative solvers.

    Notes
    -----
    The `eps_r` tolerance is both absolute and relative - the solvers
    stop when either the relative or the absolute residual is below it.
    """
    name = 'ls.scipy_iterative'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are set to default values.

        Example configuration, all items::

            solver_110 = {
                'name' : 'ls110',
                'kind' : 'ls.scipy_iterative',

                'method' : 'cg',
                'i_max' : 1000,
                'eps_r' : 1e-12,
            }
        """
        get = make_get_conf(conf, kwargs)
        common = LinearSolver.process_conf(conf)

        return Struct(method=get('method', 'cg'),
                      i_max=get('i_max', 100),
                      eps_a=None,
                      eps_r=get('eps_r', 1e-8)) + common

    def __init__(self, conf, **kwargs):
        import scipy.sparse.linalg.isolve as la

        LinearSolver.__init__(self, conf, **kwargs)

        try:
            solver = getattr( la, self.conf.method )
        except AttributeError:
            output( 'scipy solver %s does not exist!' % self.conf.method )
            output( 'using cg instead' )
            solver = la.cg
        self.solver = solver
        self.converged_reasons = {
            0 : 'successful exit',
            1 : 'number of iterations',
            -1 : 'illegal input or breakdown',
        }

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        conf = get_default(conf, self.conf)
        eps_r = get_default(eps_r, self.conf.eps_r)
        i_max = get_default(i_max, self.conf.i_max)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        sol, info = self.solver(mtx, rhs, x0=x0, tol=eps_r, maxiter=i_max)
        output('%s convergence: %s (%s)'
               % (self.conf.method,
                  info, self.converged_reasons[nm.sign(info)]))

        return sol

##
# c: 02.05.2008, r: 02.05.2008
class PyAMGSolver( LinearSolver ):
    """
    Interface to PyAMG solvers.

    Notes
    -----
    Uses relative convergence tolerance, i.e. eps_r is scaled by `||b||`.
    """
    name = 'ls.pyamg'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are set to default values.

        Example configuration, all items::

            solver_102 = {
                'name' : 'ls102',
                'kind' : 'ls.pyamg',

                'method' : 'smoothed_aggregation_solver',
                'accel' : 'cg'
                'eps_r' : 1e-12,
            }
        """
        get = make_get_conf(conf, kwargs)
        common = LinearSolver.process_conf(conf)

        return Struct(method=get('method', 'smoothed_aggregation_solver'),
                      accel = get('accel', None),
                      i_max=None, eps_a=None,
                      eps_r=get('eps_r', 1e-8)) + common

    ##
    # c: 02.05.2008, r: 02.05.2008
    def __init__( self, conf, **kwargs ):
        try:
            import pyamg
        except ImportError:
            msg =  'cannot import pyamg!'
            raise ImportError( msg )

        LinearSolver.__init__(self, conf, eps_r=conf.eps_r, mg=None, **kwargs)

        try:
            solver = getattr( pyamg, self.conf.method )
        except AttributeError:
            output( 'pyamg.%s does not exist!' % self.conf.method )
            output( 'using pyamg.smoothed_aggregation_solver instead' )
            solver = pyamg.smoothed_aggregation_solver
        self.solver = solver

        if hasattr( self, 'mtx' ):
            if self.mtx is not None:
                self.mg = self.solver( self.mtx )

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        conf = get_default(conf, self.conf)
        eps_r = get_default(eps_r, self.eps_r)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        if (self.mg is None) or (mtx is not self.mtx):
            self.mg = self.solver(mtx)
            self.mtx = mtx

        sol = self.mg.solve(rhs, x0=x0, accel=conf.accel, tol=eps_r)

        return sol

class PETScKrylovSolver( LinearSolver ):
    """
    PETSc Krylov subspace solver.

    The solver and preconditioner types are set upon the solver object
    creation. Tolerances can be overriden when called by passing a `conf`
    object.

    Notes
    -----
    Convergence is reached when `rnorm < max(eps_r * rnorm_0, eps_a)`,
    where, in PETSc, `rnorm` is by default the norm of *preconditioned*
    residual.
    """
    name = 'ls.petsc'

    _precond_sides = {None : None, 'left' : 0, 'right' : 1, 'symmetric' : 2}

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are set to default values.

        Example configuration, all items::

            solver_120 = {
                'name' : 'ls120',
                'kind' : 'ls.petsc',

                'method' : 'cg', # ksp_type
                'precond' : 'icc', # pc_type
                'precond_side' : 'left', # ksp_pc_side
                'eps_a' : 1e-12, # abstol
                'eps_r' : 1e-12, # rtol
                'i_max' : 1000, # maxits
            }
        """
        get = make_get_conf(conf, kwargs)
        common = LinearSolver.process_conf(conf)

        return Struct(method=get('method', 'cg'),
                      precond=get('precond', 'icc'),
                      precond_side=get('precond_side', None),
                      i_max=get('i_max', 100),
                      eps_a=get('eps_a', 1e-8),
                      eps_r=get('eps_r', 1e-8)) + common

    def __init__( self, conf, **kwargs ):
        try:
            import petsc4py
            petsc4py.init([])
            from petsc4py import PETSc
        except ImportError:
            msg = 'cannot import petsc4py!'
            raise ImportError( msg )

        LinearSolver.__init__(self, conf, eps_a=conf.eps_a, eps_r=conf.eps_r,
                              petsc=PETSc, pmtx=None, **kwargs)

        ksp = PETSc.KSP().create()

        ksp.setType( self.conf.method )
        ksp.getPC().setType( self.conf.precond )
        side = self._precond_sides[self.conf.precond_side]
        if side is not None:
            ksp.setPCSide(side)
        self.ksp = ksp

        self.converged_reasons = {}
        for key, val in ksp.ConvergedReason.__dict__.iteritems():
            if isinstance(val, int):
                self.converged_reasons[val] = key

    def set_matrix( self, mtx ):
        mtx = sps.csr_matrix(mtx)

        pmtx = self.petsc.Mat().createAIJ( mtx.shape,
                                           csr = (mtx.indptr,
                                                  mtx.indices,
                                                  mtx.data) )
        sol, rhs = pmtx.getVecs()
        return pmtx, sol, rhs

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        conf = get_default(conf, self.conf)
        eps_a = get_default(eps_a, self.eps_a)
        eps_r = get_default(eps_r, self.eps_r)
        i_max = get_default(i_max, self.conf.i_max)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        # There is no use in caching matrix in the solver - always set as new.
        pmtx, psol, prhs = self.set_matrix(mtx)

        ksp = self.ksp
        ksp.setOperators(pmtx)
        ksp.setFromOptions() # PETSc.Options() not used yet...
        ksp.setTolerances(atol=eps_a, rtol=eps_r, max_it=i_max)

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

    Notes
    -----
    Convergence is reached when `rnorm < max(eps_r * rnorm_0, eps_a)`,
    where, in PETSc, `rnorm` is by default the norm of *preconditioned*
    residual.
    """
    name = 'ls.petsc_parallel'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are set to default values.

        Example configuration, all items::

            solver_1 = {
                'name' : 'ls',
                'kind' : 'ls.petsc_parallel',

                'n_proc' : 5, # Number of processes to run.

                'method' : 'cg', # ksp_type
                'precond' : 'bjacobi', # pc_type
                'sub_precond' : 'icc', # sub_pc_type
                'eps_a' : 1e-12, # abstol
                'eps_r' : 1e-12, # rtol
                'i_max' : 1000, # maxits
            }
        """
        get = make_get_conf(conf, kwargs)
        common = PETScKrylovSolver.process_conf(conf, kwargs)

        return Struct(n_proc=get('n_proc', 1),
                      sub_precond=get('sub_precond', 'icc')) + common

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        import os, sys, shutil, tempfile, time
        from sfepy import base_dir, data_dir
        from sfepy.base.ioutils import ensure_path

        conf = get_default(conf, self.conf)
        eps_a = get_default(eps_a, self.eps_a)
        eps_r = get_default(eps_r, self.eps_r)
        i_max = get_default(i_max, self.conf.i_max)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        petsc = self.petsc

        # There is no use in caching matrix in the solver - always set as new.
        pmtx, psol, prhs = self.set_matrix(mtx)

        ksp = self.ksp
        ksp.setOperators(pmtx)
        ksp.setFromOptions() # PETSc.Options() not used yet...
        ksp.setTolerances(atol=eps_a, rtol=eps_r, max_it=i_max)

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

        log_filename = os.path.join(data_dir, 'tmp/sol.log')
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

class SchurComplement(ScipyDirect):
    r"""
    Schur complement.

    Solution of the linear system
    ..math::
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
    ..math::
       (D - C A^{-1} B) \cdot v = g - C A^{-1} f

    variable(s) :math:`u` are specified in "eliminate" list,
    variable(s) :math:`v` are specified in "keep" list,

    See: http://en.wikipedia.org/wiki/Schur_complement
    """
    name = 'ls.schur_complement'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Example configuration::

            solvers = {
                'ls': ('ls.schur_complement', {
                    'eliminate': ['displacement'],
                    'keep': ['pressure'],
                    'needs_problem_instance': True,
            }
        """
        get = make_get_conf(conf, kwargs)
        common = LinearSolver.process_conf(conf)

        return Struct(method=get('method', 'auto'),
                      presolve=False,
                      warn=get('warn', True),
                      i_max=None, eps_a=None, eps_r=None) + common

    def __init__(self, conf, **kwargs):
        from sfepy.fem.state import State

        ScipyDirect.__init__(self, conf, **kwargs)

        equations = self.problem.equations
        aux_state = State(equations.variables)
        bl = {1: conf.eliminate, 2: conf.keep}
        conf.idxs = {}
        for ii in [1, 2]:
            aux_state.fill(0.0)
            for jj in bl[ii]:
                idx = equations.variables.di.indx[jj]
                aux_state.vec[idx] = nm.nan

            aux_state.apply_ebc()
            vec0 = aux_state.get_reduced()
            conf.idxs[ii] = nm.where(nm.isnan(vec0))[0]

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        import scipy.sparse as scs

        conf = get_default(conf, self.conf)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        mtxi= self.orig_conf.idxs
        mtxslc_s = {}
        mtxslc_f = {}
        nn = {}

        for ii in [1, 2]:
            ptr = 0
            nn[ii] = len(mtxi[ii])
            mtxslc_s[ii] = []
            mtxslc_f[ii] = []
            while ptr < nn[ii]:
                idx0 = mtxi[ii][ptr:]
                idxrange = nm.arange(idx0[0], idx0[0] + len(idx0))
                aux = nm.where(idx0 == idxrange)[0]
                mtxslc_s[ii].append(slice(ptr + aux[0], ptr + aux[-1] + 1))
                mtxslc_f[ii].append(slice(idx0[aux][0], idx0[aux][-1] + 1))
                ptr += aux[-1] + 1

        mtxs = {}
        rhss = {}
        ress = {}
        for ir in [1, 2]:
            rhss[str(ir)] = nm.zeros((nn[ir],), dtype=nm.float64)
            ress[str(ir)] = nm.zeros((nn[ir],), dtype=nm.float64)
            for jr, idxr in enumerate(mtxslc_f[ir]):
                rhss[str(ir)][mtxslc_s[ir][jr]] = rhs[idxr]

            for ic in [1, 2]:
                mtxid = '%d%d' % (ir, ic)
                mtxs[mtxid] = nm.zeros((nn[ir], nn[ic]), dtype=nm.float64)
                for jr, idxr in enumerate(mtxslc_f[ir]):
                    for jc, idxc in enumerate(mtxslc_f[ic]):
                        mtxs[mtxid][mtxslc_s[ir][jr], mtxslc_s[ic][jc]] = mtx._get_submatrix(idxr, idxc).todense()

        invA = self.sls.splu(scs.csc_matrix(mtxs['11']))
        spC = scs.csc_matrix(mtxs['21'])
        invAB = nm.zeros((nn[1], nn[2]), dtype=nm.float64)
        for ic in xrange(nn[2]):
            invAB[:,ic] = invA.solve(mtxs['12'][:,ic])
        invAf = invA.solve(rhss['1'])
        k_rhs = rhss['2'] - spC * invAf
        ress['2'] = self.sls.spsolve(scs.csc_matrix(mtxs['22'] - spC * invAB), k_rhs)
        ress['1'] = invAf - nm.dot(invAB, ress['2'])

        res = nm.zeros_like(rhs)
        for ir in [1, 2]:
            for jr, idxr in enumerate(mtxslc_f[ir]):
                res[idxr] = ress[str(ir)][mtxslc_s[ir][jr]]

        return res

    def _presolve(self):
        if hasattr(self, 'presolve'):
            return self.presolve
        else:
            return self.conf.presolve
