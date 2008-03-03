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
varDict = vars().items()
solverTable = {}

for key, var in varDict:
    try:
        if isDerivedClass( var, LinearSolver ) or \
               isDerivedClass( var, NonlinearSolver ) or \
               isDerivedClass( var, TimeSteppingSolver ) or \
               isDerivedClass( var, EigenvalueSolver ) or \
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
