from solvers import *

from ls import *
from nls import *
from oseen import *
try:
    from optimize import *
except:
    pass

##
# 16.10.2007, c
# 17.10.2007
varDict = vars().items()
solverTable = {}

for key, var in varDict:
    try:
        if isDerivedClass( var, LinearSolver ) or \
               isDerivedClass( var, NonlinearSolver ) or \
               isDerivedClass( var, OptimizationSolver ):
            solverTable[var.name] = var
    except TypeError:
        pass
del varDict

## print solverTable
## pause()

##
# 23.10.2007, c
def anyFromConf( conf, **kwargs ):
    return solverTable[conf.kind]( conf, **kwargs )
insertStaticMethod( Solver, anyFromConf )
del anyFromConf
