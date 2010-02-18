from solvers import *

from ls import *
from nls import *
from oseen import *
from ts import *
from eigen import *
from optimize import *
from semismooth_newton import *

solver_table = find_subclasses(vars(),
                               [LinearSolver, NonlinearSolver,
                                TimeSteppingSolver, EigenvalueSolver,
                                OptimizationSolver])

def any_from_conf(conf, **kwargs):
    """Create an instance of a solver class according to the configuration."""
    return solver_table[conf.kind](conf, **kwargs)
insert_static_method(Solver, any_from_conf)
del any_from_conf
