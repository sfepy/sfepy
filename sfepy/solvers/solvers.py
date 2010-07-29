"""
Base (abstract) solver classes.
"""
from sfepy.base.base import *

class Solver(Struct):
    """
    Base class for all solver kinds. Takes care of processing of common
    configuration options.

    The factory method any_from_conf() can be used to create an instance of any
    subclass.

    The subclasses have to reimplement __init__() and __call__(). The
    subclasses that implement process_conf() have to call Solver.process_conf().
    """

    @staticmethod
    def process_conf(conf):
        """
        Ensures conf contains 'name' and 'kind'.
        """
        get = conf.get_default_attr
        name = get('name', None, 'missing "name" in options!')
        name = get('kind', None, 'missing "kind" in options!')

        return Struct(**locals())

    def __init__(self, conf, **kwargs):
        if isinstance(conf, dict):
            conf = Struct(**conf)

            if conf.get('name', None) is None:
                conf.name = 'auto_' + self.__class__.__name__

            if conf.get('kind', None) is None:
                if hasattr(self.__class__, 'name'):
                    conf.kind = self.__class__.name

        conf = self.__class__.process_conf(conf)
        Struct.__init__(self, conf=conf, **kwargs)

    def __call__(self, **kwargs):
        raise ValueError('called an abstract Solver instance!')

class LinearSolver(Solver):
    """
    Abstract linear solver class.
    """

    def __init__(self, conf, mtx=None, status=None, **kwargs):
        Solver.__init__(self, conf=conf, mtx=mtx, status=status, **kwargs)
    
    def __call__(self, rhs, x0=None, conf=None, mtx=None, status=None):
        raise ValueError('called an abstract LinearSolver instance!')

class NonlinearSolver(Solver):
    """
    Abstract nonlinear solver class.
    """

    def __init__(self, conf, fun=None, fun_grad=None, lin_solver=None,
                  status=None, **kwargs):
        Solver.__init__(self, conf=conf, fun=fun, fun_grad=fun_grad,
                        lin_solver=lin_solver, status=status, **kwargs)
    
    def __call__(self, state0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, status=None):
        raise ValueError('called an abstract NonlinearSolver instance!')

class TimeSteppingSolver(Solver):
    """
    Abstract time stepping solver class.
    """

    def __init__(self, conf, step_fun=None, step_args=None, **kwargs):
        Solver.__init__(self, conf=conf,
                        step_fun=step_fun, step_args=step_args, **kwargs)

    def __call__(self, state0=None, conf=None, step_fun=None, step_args=None ):
        raise ValueError('called an abstract TimeSteppingSolver instance!')

class OptimizationSolver(Solver):
    """
    Abstract optimization solver class.
    """

    def __init__(self, conf, obj_fun=None, obj_fun_grad=None, status=None,
                 obj_args=None,  **kwargs ):
        Solver.__init__(self, conf=conf, obj_fun=obj_fun,
                        obj_fun_grad=obj_fun_grad, status=status,
                        obj_args=obj_args, **kwargs)

    def __call__(self, state0, conf=None, obj_fun=None, obj_fun_grad=None,
                 status=None, obj_args=None):
        raise ValueError('called an abstract OptimizationSolver instance!')

class EigenvalueSolver(Solver):
    """
    Abstract eigenvalue solver class.
    """

    def __init__(self, conf, mtx_a=None, mtx_b=None, n_eigs=None,
                 eigenvectors=None, status=None):
        Solver.__init__(self, conf=conf, mtx_a=mtx_a, mtx_b=mtx_b,
                        n_eigs=n_eigs, eigenvectors=eigenvectors,
                        status=status)
                         
    def __call__(self, mtx_a, mtx_b=None, n_eigs=None,
                 eigenvectors=None, status=None, conf=None):
        raise ValueError('called an abstract EigenvalueSolver instance!')

    def _to_array(self, mtx_a, mtx_b=None):
        if hasattr(mtx_a, 'toarray'):
            mtx_a = mtx_a.toarray()
        if mtx_b is not None:
            if hasattr(mtx_b, 'toarray'):
                mtx_b = mtx_b.toarray()
        return mtx_a, mtx_b
