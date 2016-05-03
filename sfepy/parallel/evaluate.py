"""
PETSc-related parallel evaluation of problem equations.
"""
import numpy as nm

from sfepy.base.base import Struct
import sfepy.parallel.parallel as pp
from sfepy.discrete.evaluate import BasicEvaluator

class PETScParallelEvaluator(BasicEvaluator):
    """
    The parallel evaluator of the problem equations for
    :class:`PETScNonlinearSolver <sfepy.solvers.nls.PETScNonlinearSolver>`.

    Its methods can be used as the function and Jacobian callbacks of the PETSc
    SNES (Scalable Nonlinear Equations Solvers).

    Notes
    -----
    Assumes ``problem.active_only == False``.
    """

    def __init__(self, problem, pdofs, drange, is_overlap, psol,
                 comm, matrix_hook=None, verbose=False):
        BasicEvaluator.__init__(self, problem, matrix_hook=matrix_hook)
        Struct.__init__(self, pdofs=pdofs, drange=drange, is_overlap=is_overlap,
                        comm=comm, verbose=verbose)

        variables = problem.get_variables()
        di = variables.di

        ebc_rows = []
        for ii, variable in enumerate(variables.iter_state(ordered=True)):
            eq_map = variable.eq_map
            ebc_rows.append(eq_map.eq_ebc + di.ptr[ii])
        self.ebc_rows = nm.concatenate(ebc_rows)

        self.psol_i = pp.create_local_petsc_vector(pdofs)

        self.gather, self.scatter = pp.create_gather_scatter(pdofs, self.psol_i,
                                                             psol, comm=comm)

    def eval_residual(self, snes, psol, prhs):
        self.scatter(self.psol_i, psol)

        rhs_if = BasicEvaluator.eval_residual(self, self.psol_i[...],
                                              is_full=True)

        pp.assemble_rhs_to_petsc(prhs, rhs_if, self.pdofs, self.drange,
                                 self.is_overlap,
                                 self.comm, verbose=self.verbose)

    def eval_tangent_matrix(self, snes, psol, pmtx, ppmtx):
        self.scatter(self.psol_i, psol)

        mtx_if = BasicEvaluator.eval_tangent_matrix(self, self.psol_i[...],
                                                    is_full=True)

        pp.apply_ebc_to_matrix(mtx_if, self.ebc_rows)
        pp.assemble_mtx_to_petsc(pmtx, mtx_if, self.pdofs, self.drange,
                                 self.is_overlap,
                                 self.comm, verbose=self.verbose)
