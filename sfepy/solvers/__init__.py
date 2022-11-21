from __future__ import absolute_import
import os
import sfepy
from sfepy.base.base import load_classes, insert_static_method
from .solvers import *
from .eigen import eig
from .auto_fallback import AutoFallbackSolver

solver_files = sfepy.get_paths('sfepy/solvers/*.py')
remove = ['setup.py', 'solvers.py', 'ls_mumps_parallel.py']
solver_files = [name for name in solver_files
                if os.path.basename(name) not in remove]
solver_table = load_classes(solver_files,
                            [AutoFallbackSolver,
                             LinearSolver, NonlinearSolver,
                             TimeStepController, TimeSteppingSolver,
                             EigenvalueSolver, QuadraticEVPSolver,
                             OptimizationSolver], package_name='sfepy.solvers')


def register_solver(cls):
    """
    Register a custom solver.
    """
    solver_table[cls.name] = cls

def any_from_conf(conf, **kwargs):
    """Create an instance of a solver class according to the configuration."""
    try:
        return solver_table[conf.kind](conf, **kwargs)
    except Exception:
        raise
insert_static_method(Solver, any_from_conf)
del any_from_conf
del sfepy
