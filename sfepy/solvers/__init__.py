import os
import sfepy
from solvers import *
from eigen import eig

solver_files = sfepy.get_paths('sfepy/solvers/*.py')
remove = ['setup.py', 'solvers.py', 'generic.py']
solver_files = [name for name in solver_files
                if os.path.basename(name) not in remove]
solver_table = load_classes(solver_files,
                            [LinearSolver, NonlinearSolver,
                             TimeSteppingSolver, EigenvalueSolver,
                             OptimizationSolver], package_name='sfepy.solvers')

def any_from_conf(conf, **kwargs):
    """Create an instance of a solver class according to the configuration."""
    return solver_table[conf.kind](conf, **kwargs)
insert_static_method(Solver, any_from_conf)
del any_from_conf
del sfepy
