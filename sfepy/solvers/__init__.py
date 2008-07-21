from solvers import *

from ls import *
from nls import *
from oseen import *
from ts import *
from eigen import *
try:
    from optimize import *
except:
    pass

##
# c: 16.10.2007, r: 03.03.2008
var_dict = vars().items()
solver_table = {}

for key, var in var_dict:
    try:
        if is_derived_class( var, LinearSolver ) or \
               is_derived_class( var, NonlinearSolver ) or \
               is_derived_class( var, TimeSteppingSolver ) or \
               is_derived_class( var, EigenvalueSolver ) or \
               is_derived_class( var, OptimizationSolver ):
            solver_table[var.name] = var
    except TypeError:
        pass
del var_dict

## print solver_table
## pause()

##
# 23.10.2007, c
def any_from_conf( conf, **kwargs ):
    return solver_table[conf.kind]( conf, **kwargs )
insert_static_method( Solver, any_from_conf )
del any_from_conf
