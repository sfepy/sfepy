from __future__ import absolute_import
from sfepy.base.base import Struct
from sfepy.solvers.solvers import Solver, use_first_available


class AutoFallbackSolver(Solver):
    """
    Base class for virtual solvers with the automatic fallback.
    """
    _ls_solvers = []

    def __new__(cls, conf, **kwargs):
        """
        Choose a available solver from `self._ls_solvers`.

        Parameters
        ----------
        conf : dict
            The solver configuration.
        """
        ls_solvers = [(ls, Struct(**conf) + Struct(kind=ls))
                      for ls, conf in cls._ls_solvers]

        return use_first_available(ls_solvers)


class AutoDirect(AutoFallbackSolver):
    """The automatically selected linear direct solver.

    The first available solver from the following list is used:
    `ls.mumps <sfepy.solvers.ls.MUMPSSolver>`,
    `ls.scipy_umfpack <sfepy.solvers.ls.ScipyUmfpack>` and
    `ls.scipy_superlu <sfepy.solvers.ls.ScipySuperLU>`.
    """
    name = 'ls.auto_direct'

    _ls_solvers = [
        ('ls.mumps', {}),
        ('ls.scipy_umfpack', {}),
        ('ls.scipy_superlu', {})
    ]


class AutoIterative(AutoFallbackSolver):
    """The automatically selected linear iterative solver.

    The first available solver from the following list is used:
    `ls.petsc <sfepy.solvers.ls.PETScKrylovSolver>` and
    `ls.scipy_iterative <sfepy.solvers.ls.ScipyIterative>`
    """
    name = 'ls.auto_iterative'

    _ls_solvers = [
        ('ls.petsc', {'method': 'cg', 'precond': 'icc'}),
        ('ls.scipy_iterative', {'method': 'cg'}),
    ]
