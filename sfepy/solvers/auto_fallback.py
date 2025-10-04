from sfepy.base.base import Struct
from sfepy.solvers.solvers import Solver, use_first_available


class AutoFallbackSolver(Solver):
    """
    Base class for virtual solvers with the automatic fallback.
    """
    _ls_solvers = []

    def __new__(cls, conf, **kwargs):
        """
        Choose an available solver from `self._ls_solvers`.

        Parameters
        ----------
        conf : dict or Struct
            The solver configuration.
        **kwargs : keyword arguments
            Additional solver options, see the particular __init__() methods.
        """
        if isinstance(conf, Struct):
            dconf = conf.to_dict()

        else:
            dconf = conf

        dconf.pop('kind', None)
        ls_solvers = [(ls, Struct(**_conf) + Struct(kind=ls) + Struct(**dconf))
                      for ls, _conf in cls._ls_solvers]

        return use_first_available(ls_solvers, **kwargs)


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
