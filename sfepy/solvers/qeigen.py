"""
Quadratic eigenvalue problem solvers.
"""
from __future__ import absolute_import
import time

import numpy as nm
import scipy.sparse as sps

from sfepy.base.base import get_default, structify
from sfepy.solvers.solvers import SolverMeta, Solver

def standard_call(call):
    """
    Decorator handling argument preparation and timing for quadratic
    eigensolvers.
    """
    def _standard_call(self, mtx_m, mtx_d, mtx_k, n_eigs=None,
                       eigenvectors=None, status=None, conf=None, **kwargs):
        tt = time.clock()

        conf = get_default(conf, self.conf)
        mtx_m = get_default(mtx_m, self.mtx_m)
        mtx_d = get_default(mtx_d, self.mtx_d)
        mtx_k = get_default(mtx_k, self.mtx_k)
        n_eigs = get_default(n_eigs, self.n_eigs)
        eigenvectors = get_default(eigenvectors, self.eigenvectors)
        status = get_default(status, self.status)

        result = call(self, mtx_m, mtx_d, mtx_k,
                      n_eigs, eigenvectors, status, conf,
                      **kwargs)

        ttt = time.clock() - tt
        if status is not None:
            status['time'] = ttt

        return result

    return _standard_call

class QuadraticEVPSolver(Solver):
    """
    Abstract quadratic eigenvalue problem solver class.
    """

    def __init__(self, conf, mtx_m=None, mtx_d=None, mtx_k=None, n_eigs=None,
                 eigenvectors=None, status=None, context=None, **kwargs):
        Solver.__init__(self, conf=conf, mtx_m=mtx_m, mtx_d=mtx_d, mtx_k=mtx_k,
                        n_eigs=n_eigs, eigenvectors=eigenvectors,
                        status=status, context=context)
        solver_conf = structify(conf.solver)
        self.solver = Solver.any_from_conf(solver_conf)

    def __call__(self, mtx_m, mtx_d, mtx_k, n_eigs=None,
                 eigenvectors=None, status=None, conf=None):
        raise ValueError('called an abstract QuadraticEVPSolver instance!')

class LQuadraticEVPSolver(QuadraticEVPSolver):
    """
    Quadratic eigenvalue problem solver based on the problem linearization.

    (w^2 M + w D + K) x = 0.
    """
    name = 'eig.qevp'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', "{'companion', 'cholesky'}", 'companion', False,
         'The linearization method.'),
        ('solver', 'dict', {'kind': 'eig.scipy', 'method': 'eig'}, False,
         """The configuration of an eigenvalue solver for
            the linearized problem (A - w B) x = 0."""),
        ('mode', "{'normal', 'inverted'}", 'normal', False,
         'Solve either A - w B (normal), or B - 1/w A (inverted).'),
    ]

    @standard_call
    def __call__(self, mtx_m, mtx_d, mtx_k, n_eigs=None,
                 eigenvectors=None, status=None, conf=None):

        if conf.method == 'companion':
            mtx_eye = -sps.eye(mtx_m.shape[0], dtype=mtx_m.dtype)

            mtx_a = sps.bmat([[mtx_d, mtx_k],
                              [mtx_eye, None]])
            mtx_b = sps.bmat([[-mtx_m, None],
                              [None, mtx_eye]])

        elif conf.method == 'cholesky':
            from sksparse.cholmod import cholesky

            factor = cholesky(mtx_m)
            perm = factor.P()
            ir = nm.arange(len(perm))
            mtx_p = sps.coo_matrix((nm.ones_like(perm), (ir, perm)))
            mtx_l = mtx_p.T * factor.L()

            mtx_eye = sps.eye(mtx_l.shape[0], dtype=nm.float64)

            mtx_a = sps.bmat([[-mtx_k, None],
                              [None, mtx_eye]])
            mtx_b = sps.bmat([[mtx_d, mtx_l],
                              [mtx_l.T, None]])

        else:
            raise ValueError('unknown method! (%s)' % conf.method)

        if conf.mode == 'normal':
            out = self.solver(mtx_a, mtx_b, n_eigs=n_eigs,
                              eigenvectors=eigenvectors, status=status)

            if eigenvectors:
                eigs, vecs = out
                out = (eigs, vecs[:mtx_m.shape[0], :])

        else:
            out = self.solver(mtx_b, mtx_a, n_eigs=n_eigs,
                              eigenvectors=eigenvectors, status=status)

            if eigenvectors:
                eigs, vecs = out
                out = (1.0 / eigs, vecs[:mtx_m.shape[0], :])

            else:
                out = 1.0 / out

        return out
