from __future__ import absolute_import

import numpy as nm
import scipy.sparse as sps

from sfepy.base.base import get_default, try_imports, Struct
from sfepy.base.timing import Timer
from sfepy.solvers.solvers import Solver, EigenvalueSolver
import six
from six.moves import range

def eig(mtx_a, mtx_b=None, n_eigs=None, eigenvectors=True,
        return_time=None, solver_kind='eig.scipy', **ckwargs):
    """
    Utility function that constructs an eigenvalue solver given by
    `solver_kind`, calls it and returns its output.
    """
    kwargs = {'name' : 'aux', 'kind' : solver_kind}
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
        timer = Timer(start=True)

        conf = get_default(conf, self.conf)
        mtx_a = get_default(mtx_a, self.mtx_a)
        mtx_b = get_default(mtx_b, self.mtx_b)
        n_eigs = get_default(n_eigs, self.n_eigs)
        eigenvectors = get_default(eigenvectors, self.eigenvectors)
        status = get_default(status, self.status)

        if n_eigs == 0:
            result =  self._ret_zero(mtx_a, eigenvectors=eigenvectors)

        else:
            result = call(self, mtx_a, mtx_b, n_eigs, eigenvectors, status,
                          conf, **kwargs)

        elapsed = timer.stop()
        if status is not None:
            status['time'] = elapsed

        return result

    return _standard_call

class ScipyEigenvalueSolver(EigenvalueSolver):
    """
    SciPy-based solver for both dense and sparse problems.

    The problem is consirered sparse if `n_eigs` argument is not None.
    """
    name = 'eig.scipy'

    _parameters = [
        ('method', "{'eig', 'eigh', 'eigs', 'eigsh'}", 'eigs', False,
         """The method for solving general or symmetric eigenvalue problems:
            for dense problems :func:`eig()` or :func:`eigh()` can be used, for
            sparse problems :func:`eigs()` or :func:`eigsh()` should be
            used."""),
        ('which', "'LM' | 'SM' | 'LR' | 'SR' | 'LI' | 'SI'", 'SM', False,
         """Which eigenvectors and eigenvalues to find,
            see :func:`scipy.sparse.linalg.eigs()`
            or :func:`scipy.sparse.linalg.eigsh()`. For dense problmes,
            only 'LM' and 'SM' can be used"""),
        ('linear_solver', "({'ls.cholesky', 'ls.mumps', ...}, ls_conf)",
         None, False,
         """The method used to construct an inverse linear operator. If None, the
            eigenvalue solver will solve the linear system internally."""),
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

        if conf.method in ('eig', 'eigh'):
            mtx_a, mtx_b = self._to_array(mtx_a, mtx_b)
            if conf.method == 'eig':
                out = self.sla.eig(mtx_a, mtx_b, right=eigenvectors, **kwargs)

            else:
                out = self.sla.eigh(mtx_a, mtx_b,
                                    eigvals_only=not eigenvectors, **kwargs)

        else:
            eig = self.ssla.eigs if conf.method == 'eigs' else self.ssla.eigsh
            if conf.linear_solver is not None:
                import sfepy.solvers.ls as ls
                from sfepy.solvers import solver_table
                import scipy.sparse as sps
                from sfepy.base.base import Struct

                ls_solvers = {}
                for k, v in solver_table.items():
                    if not k.startswith('ls'):
                        continue

                    aux = [par[0] == 'use_presolve' for par in v._parameters]
                    if nm.any(aux):
                        ls_solvers[k] = v

                fake_mtx_a = sps.csc_matrix(mtx_a.shape)
                solver_name, solver_conf = conf.linear_solver
                ls_conf = Struct(use_presolve=True, **solver_conf)
                ls = ls_solvers[solver_name](ls_conf)
                ls.mtx = mtx_a
                ls.presolve(mtx_a)
                matvec = ls.__call__

                inv_op_a = self.ssla.LinearOperator(mtx_a.shape, matvec=matvec)
                out = eig(fake_mtx_a, M=mtx_b, k=n_eigs, which=conf.which,
                          OPinv=inv_op_a,
                          return_eigenvectors=eigenvectors, **kwargs)
            else:
                out = eig(mtx_a, M=mtx_b, k=n_eigs, which=conf.which,
                          return_eigenvectors=eigenvectors, **kwargs)

        if eigenvectors:
            eigs = out[0]

        else:
            eigs = out

        if nm.iscomplexobj(eigs):
            ii = nm.argsort(nm.linalg.norm(eigs[:, None], axis=1))

        else:
            ii = nm.argsort(eigs)

        if n_eigs is not None and (conf.method in ('eig', 'eigh')):
            if conf.which == 'SM':
                ii = ii[:n_eigs]

            elif conf.which == 'LM':
                ii = ii[:-n_eigs-1:-1]

            else:
                raise ValueError("only 'LM' or 'SM' can be used with dense"
                                 " problems! (%s)" % conf.which)

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

    def __init__(self, conf, **kwargs):
        EigenvalueSolver.__init__(self, conf, **kwargs)

        try:
            import scipy.linalg.lapack as llapack
        except ImportError:
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

        # Fix output order of scipy.linalg.lapack functions.
        if out[0].ndim == 2:
            out = (out[1], out[0]) + out[2:]

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

        from scipy.sparse.linalg import lobpcg
        self.lobpcg = lobpcg

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None):

        if n_eigs is None:
            n_eigs = mtx_a.shape[0]
        else:
            n_eigs = min(n_eigs, mtx_a.shape[0])

        rng = nm.random.default_rng(12345)
        x = rng.normal(size=(mtx_a.shape[0], n_eigs))

        out = self.lobpcg(mtx_a, x, mtx_b,
                          M=conf.precond,
                          tol=conf.eps_a, maxiter=conf.i_max,
                          largest=conf.largest,
                          verbosityLevel=conf.verbose)

        if not eigenvectors:
            out = out[0]

        return out

def init_slepc_args():
    try:
        import sys, slepc4py

    except ImportError:
        return

    argv = [arg for arg in sys.argv if arg not in ['-h', '--help']]
    slepc4py.init(argv)

class SLEPcEigenvalueSolver(EigenvalueSolver):
    """
    General SLEPc eigenvalue problem solver.
    """
    name = 'eig.slepc'

    _parameters = [
        ('method', 'str', 'krylovschur', False,
         'The actual solver to use.'),
        ('problem', 'str', 'gnhep', False,
         """The problem type: Hermitian (hep), non-Hermitian (nhep), generalized
         Hermitian (ghep), generalized non-Hermitian (gnhep), generalized
         non-Hermitian with positive semi-definite B (pgnhep), and generalized
         Hermitian-indefinite (ghiep)."""),
        ('i_max', 'int', 20, False,
         'The maximum number of iterations.'),
        ('eps', 'float', None, False,
         'The convergence tolerance.'),
        ('conv_test', '{"abs", "rel", "norm", "user"}, ', 'abs', False,
         'The type of convergence test.'),
        ('which', """{'largest_magnitude', 'smallest_magnitude',
        'largest_real', 'smallest_real',
        'largest_imaginary', 'smallest_imaginary', 'target_magnitude',
        'target_real', 'target_imaginary', 'all', 'which_user'}""",
         'largest_magnitude', False,
         'Which eigenvectors and eigenvalues to find.'),
        ('*', '*', None, False,
         'Additional parameters supported by the method.'),
    ]

    def __init__(self, conf, comm=None, context=None, **kwargs):
        if comm is None:
            init_slepc_args()

        from petsc4py import PETSc as petsc
        from slepc4py import SLEPc as slepc

        EigenvalueSolver.__init__(self, conf, petsc=petsc, slepc=slepc,
                                  comm=comm, context=context, **kwargs)

    def create_eps(self, options=None, comm=None):
        optDB = self.petsc.Options()

        if options is not None:
            for key, val in six.iteritems(options):
                optDB[key] = val

        es = self.slepc.EPS()
        es.create(comm)

        return es

    def create_petsc_matrix(self, mtx, comm=None):
        if mtx is None or isinstance(mtx, self.petsc.Mat):
            pmtx = mtx

        else:
            mtx = sps.csr_matrix(mtx)

            pmtx = self.petsc.Mat()
            pmtx.createAIJ(mtx.shape, csr=(mtx.indptr, mtx.indices, mtx.data),
                           comm=comm)

        return pmtx

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None, comm=None, context=None):
        solver_kwargs = self.build_solver_kwargs(conf)

        pmtx_a = self.create_petsc_matrix(mtx_a, comm=comm)
        pmtx_b = self.create_petsc_matrix(mtx_b, comm=comm)

        es = self.create_eps(options=solver_kwargs, comm=comm)
        es.setType(conf.method)
        es.setProblemType(getattr(es.ProblemType, conf.problem.upper()))
        es.setDimensions(nev=n_eigs)
        es.setTolerances(tol=conf.eps, max_it=conf.i_max)
        es.setOperators(pmtx_a, pmtx_b)
        es.setConvergenceTest(getattr(es.Conv, conf.conv_test.upper()))
        es.setWhichEigenpairs(getattr(es.Which, conf.which.upper()))
        es.setFromOptions()

        es.solve()

        n_converged = es.getConverged()
        if status is not None:
            status['n_iter'] = es.getIterationNumber()
            status['n_converged'] = n_converged

        if not eigenvectors:
            out = nm.array([es.getEigenvalue(ii) for ii in range(n_converged)])

        else:
            vr, vi = pmtx_a.createVecs()
            eigs = []
            vrs, vis = [], []
            is_real = True
            for ii in range(n_converged):
                val = es.getEigenpair(ii, vr, vi)
                eigs.append(val if val.imag != 0 else val.real)
                vrs.append(vr.getArray().copy())
                vis.append(vi.getArray().copy())
                if is_real and nm.sum(nm.abs(vis[-1])) > 0.0:
                    is_real = False

            eigs = nm.array(eigs)
            vecs = nm.array(vrs).T
            if not is_real:
                vecs += 1j * nm.array(vis)

            out = (eigs, vecs)

        return out

class MatlabEigenvalueSolver(EigenvalueSolver):
    """
    Matlab eigenvalue problem solver.
    """
    name = 'eig.matlab'

    _parameters = [
        ('method', """{'eig', 'eigs', None}""", 'eigs', False,
         """The solution method. Note that eig() function cannot be used for
         all inputs. If `n_eigs` is not None, eigs() is used regardless of
         this parameter."""),
        ('balance', """{'balance', 'nobalance'}""", 'balance', False,
         'The balance option for eig().'),
        ('algorithm', """{'chol', 'qz'}""", 'chol', False,
         'The algorithm option for eig().'),
        ('which',
         """{'lm', 'sm', 'la', 'sa', 'be' 'lr', 'sr', 'li', 'si', sigma}""",
         'lm', False,
         'Which eigenvectors and eigenvalues to find with eigs().'),
        ('*', '*', None, False,
         'Additional parameters supported by eigs().'),
    ]

    def __init__(self, conf, comm=None, context=None, **kwargs):
        import matlab.engine as me

        EigenvalueSolver.__init__(self, conf, me=me, context=context,
                                  **kwargs)

    def solver_call(self, mtx_filename, eigs_filename, which=None):
        import os
        import scipy.io as sio

        eng = self.me.start_matlab()
        eng.cd(os.path.dirname(__file__))
        eng.matlab_eig(mtx_filename, eigs_filename)
        eng.quit()

        return sio.loadmat(eigs_filename)

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None, comm=None, context=None):
        import os
        import shutil
        import tempfile
        import scipy.io as sio

        solver_kwargs = self.build_solver_kwargs(conf)

        dirname = tempfile.mkdtemp()
        mtx_filename = os.path.join(dirname, 'matrices.mat')
        eigs_filename = os.path.join(dirname, 'eigs.mat')
        sio.savemat(mtx_filename, {
            'A' : mtx_a,
            'B' : mtx_b if mtx_b is not None else 'None',
            'n_eigs' : n_eigs if n_eigs is not None else 'None',
            'eigenvectors' : (eigenvectors
                              if eigenvectors is not None else False),
            'method' : conf.method,
            'balance' : conf.balance,
            'algorithm' : conf.algorithm,
            'which' : conf.which,
            'verbose' : conf.verbose,
            'eigs_options' : solver_kwargs,
        })

        evp = self.solver_call(mtx_filename, eigs_filename, conf.which)

        shutil.rmtree(dirname)

        out = evp['vals'][:, 0]
        if eigenvectors:
            out =  (out, evp['vecs'])

        return out


class OctaveEigenvalueSolver(MatlabEigenvalueSolver):
    """
    Octave eigenvalue problem solver.
    """
    name = 'eig.octave'

    def __init__(self, conf, comm=None, context=None, **kwargs):
        from oct2py import octave as oe

        EigenvalueSolver.__init__(self, conf, oe=oe, context=context,
                                  **kwargs)

    def solver_call(self, mtx_filename, eigs_filename, which):
        import os
        import scipy.io as sio

        self.oe.addpath(os.path.dirname(__file__))
        self.oe.matlab_eig(mtx_filename, eigs_filename)
        self.oe.exit()

        evp = sio.loadmat(eigs_filename)
        evals = evp['vals']

        which = which.lower()
        sort_funs = {
            'lm': (nm.abs, -1),
            'sm': (nm.abs, 1),
            'la': (nm.real, -1),
            'sa': (nm.real, 1),
            'lr': (nm.real, -1),
            'sr': (nm.real, 1),
            'li': (nm.imag, -1),
            'si': (nm.imag, 1),
            'be': (nm.real, 1),
        }

        sfun, order = sort_funs[which]
        idxs = nm.argsort(sfun(evals), axis=0)
        idxs = idxs.ravel()[::order]
        if which == 'be':
            idxs = nm.hstack([idxs[::-1], idxs])
            k2 = len(idxs) // 2
            idxs = idxs[k2:(k2 + len(idxs))]

        out = {'vals': evals[idxs, :]}
        if 'vecs' in evp:
            out['vecs'] = evp['vecs'][:, idxs]

        return out


class PrimmeEigenvalueSolver(EigenvalueSolver):
    """
    PRIMME eigenvalue problem solver.

    https://github.com/primme/primme

    Installation: pip install primme
    """
    name = 'eig.primme'

    _parameters = [
        ('which',
         "{'LM', 'SM', 'LA', 'SA', 'CLT', 'CGT'}", 'LM', False,
         'Which eigenvectors and eigenvalues to find.'),
        ('sigma', 'float', None, False,
         'Find eigenvalues near sigma.'),
        ('maxiter', 'int', None, False, 'Maximum number of iterations.'),
        ('tol', 'float', 0, False,
         'Tolerance for eigenpairs (stopping criterion).'),
        ('*', '*', None, False,
         'Additional parameters supported by eigsh().'),
    ]

    def __init__(self, conf, comm=None, context=None, **kwargs):
        import primme

        EigenvalueSolver.__init__(self, conf, primme=primme, context=context,
                                  **kwargs)

    @standard_call
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None, eigenvectors=None,
                 status=None, conf=None, comm=None, context=None):
        if n_eigs is None:
            n_eigs = mtx_a.shape[0]

        solver_kwargs = self.build_solver_kwargs(conf)

        out = self.primme.eigsh(mtx_a, n_eigs, M=mtx_b,
                                which=conf.which.upper(),
                                tol=conf.tol,
                                maxiter=conf.maxiter,
                                sigma=conf.sigma,
                                return_eigenvectors=eigenvectors,
                                **solver_kwargs)

        return out
