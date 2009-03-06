from sfepy.base.base import *

##
# 10.10.2007, c
class Solver( Struct ):
    def process_conf( conf ):
        """
        Ensures conf contains 'name' and 'kind'.
        """
        get = conf.get_default_attr
        name = get( 'name', None, 'missing "name" in options!' )
        name = get( 'kind', None, 'missing "kind" in options!' )

        return Struct( **locals() )
    process_conf = staticmethod( process_conf )

    ##
    # 10.10.2007, c
    def __init__( self, conf, **kwargs ):
        conf = self.__class__.process_conf( conf )
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
    def __init__( self, conf, fun = None, fun_grad = None, lin_solver = None,
                  status = None, **kwargs ):
        Solver.__init__( self, conf = conf, fun = fun, fun_grad = fun_grad,
                         lin_solver = lin_solver, status = status,
                         **kwargs )
    
    ##
    # 10.10.2007, c
    def __call__( self, state0, conf = None, fun = None, fun_grad = None,
                  lin_solver = None, status = None ):
        print 'called an abstract NonlinearSolver instance!'
        raise ValueError

##
# c: 06.02.2008, r: 06.02.2008
class TimeSteppingSolver( Solver ):
    ##
    # c: 06.02.2008, r: 06.02.2008
    def __init__( self, conf, step_fun = None, step_args = None, **kwargs ):
        Solver.__init__( self, conf = conf,
                         step_fun = step_fun, step_args = step_args, **kwargs )

    ##
    # c: 06.02.2008, r: 06.02.2008
    def __call__( self, state0 = None, conf = None,
                  step_fun = None, step_args = None ):
        print 'called an abstract TimeSteppingSolver instance!'
        raise ValueError

##
# 17.10.2007, c
class OptimizationSolver( Solver ):
    ##
    # 17.10.2007, c
    def __init__( self, conf, obj_fun = None, obj_fun_grad = None, status = None,
                  obj_args = None,  **kwargs ):
        Solver.__init__( self, conf = conf, obj_fun = obj_fun,
                         obj_fun_grad = obj_fun_grad, status = status,
                         obj_args = obj_args, **kwargs )
    ##
    # 17.10.2007, c
    def __call__( self, state0, conf = None, obj_fun = None, obj_fun_grad = None,
                  status = None, obj_args = None ):
        print 'called an abstract OptimizationSolver instance!'
        raise ValueError

##
# c: 03.03.2008, r: 03.03.2008
class EigenvalueSolver( Solver ):
    ##
    # c: 03.03.2008, r: 03.03.2008
    def __init__( self, conf, mtx_a = None, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None ):
        Solver.__init__( self, conf = conf, mtx_a = mtx_a, mtx_b = mtx_b,
                         n_eigs = n_eigs, eigenvectors = eigenvectors,
                         status = status )
                         
    ##
    # c: 03.03.2008, r: 03.03.2008
    def __call__( self, mtx_a, mtx_b = None, n_eigs = None,
                  eigenvectors = None, status = None, conf = None ):
        print 'called an abstract EigenvalueSolver instance!'
        raise ValueError

    ##
    # c: 08.04.2008, r: 08.04.2008
    def _to_array( self, mtx_a, mtx_b = None ):
        if hasattr( mtx_a, 'toarray' ):
            mtx_a = mtx_a.toarray()
        if mtx_b is not None:
            if hasattr( mtx_b, 'toarray' ):
                mtx_b = mtx_b.toarray()
        return mtx_a, mtx_b
