from sfe.base.base import *

##
# 10.10.2007, c
class Solver( Struct ):
    ##
    # 10.10.2007, c
    def __init__( self, conf, **kwargs ):
        Struct.__init__( self, conf = conf, **kwargs )

    ##
    # 10.10.2007, c
    def __call__( self, **kwargs ):
        print 'called an abstract Solver instance!'
        raise ValueError

##
# 10.10.2007, c
class LinearSolver( Solver ):
    ##
    # 10.10.2007, c
    def __init__( self, conf, mtx = None, status = None, **kwargs ):
        Solver.__init__( self, conf = conf, mtx = mtx, status = status,
                         **kwargs )
    
    ##
    # 10.10.2007, c
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        print 'called an abstract LinearSolver instance!'
        raise ValueError

##
# 10.10.2007, c
class NonlinearSolver( Solver ):
    ##
    # 10.10.2007, c
    def __init__( self, conf, evaluator = None, linSolver = None,
                  status = None, **kwargs ):
        Solver.__init__( self, conf = conf, evaluator = evaluator,
                         linSolver = linSolver, status = status,
                         **kwargs )
    
    ##
    # 10.10.2007, c
    def __call__( self, state0, conf = None, evaluator = None,
                  linSolver = None, status = None ):
        print 'called an abstract NonlinearSolver instance!'
        raise ValueError

##
# 17.10.2007, c
class OptimizationSolver( Solver ):
    ##
    # 17.10.2007, c
    def __init__( self, conf, objFun = None, objFunGrad = None, status = None,
                  objArgs = None,  **kwargs ):
        Solver.__init__( self, conf = conf, objFun = objFun,
                         objFunGrad = objFunGrad, status = status,
                         objArgs = objArgs, **kwargs )
    ##
    # 17.10.2007, c
    def __call__( self, state0, conf = None, objFun = None, objFunGrad = None,
                  status = None, objArgs = None ):
        print 'called an abstract OptimizationSolver instance!'
        raise ValueError
