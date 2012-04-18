"""
Base (abstract) solver classes.
"""
from sfepy.base.base import get_default, Struct

def make_get_conf(conf, kwargs):
    def _get_conf_item(name, default=None, msg_if_none=None):
        return kwargs.get(name,
                          conf.get_default_attr(name, default=default,
                                                msg_if_none=msg_if_none))

    return _get_conf_item

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
    def process_conf(conf, kwargs=None):
        """
        Ensures conf contains 'name' and 'kind'.
        """
        get = conf.get_default_attr
        name = get('name', None, 'missing "name" in options!')
        kind = get('kind', None, 'missing "kind" in options!')
        verbose = get('verbose', False)

        return Struct(name=name, kind=kind, verbose=verbose)

    def __init__(self, conf=None, **kwargs):
        if conf is None:
            conf = Struct()

        elif isinstance(conf, dict):
            conf = Struct(**conf)

        if conf.get('name', None) is None:
            conf.name = 'auto_' + self.__class__.__name__

        if conf.get('kind', None) is None:
            if hasattr(self.__class__, 'name'):
                conf.kind = self.__class__.name

            else:
                raise ValueError('solver kind cannot be determined!')

        new_conf = self.process_conf(conf, kwargs)
        Struct.__init__(self, conf=new_conf, orig_conf=conf, **kwargs)

    def __call__(self, **kwargs):
        raise ValueError('called an abstract Solver instance!')

class LinearSolver(Solver):
    """
    Abstract linear solver class.
    """
    def __init__(self, conf, mtx=None, status=None, **kwargs):
        Solver.__init__(self, conf=conf, mtx=mtx, status=status, **kwargs)

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, **kwargs):
        raise ValueError('called an abstract LinearSolver instance!')

    def get_tolerance(self):
        """
        Return tuple `(eps_a, eps_r)` of absolute and relative tolerance
        settings. Either value can be `None`, meaning that the solver
        does not use that setting.
        """
        return self.conf.eps_a, self.conf.eps_r

class NonlinearSolver(Solver):
    """
    Abstract nonlinear solver class.
    """

    def __init__(self, conf, fun=None, fun_grad=None, lin_solver=None,
                 iter_hook=None, status=None, **kwargs):
        Solver.__init__(self, conf=conf, fun=fun, fun_grad=fun_grad,
                        lin_solver=lin_solver, iter_hook=iter_hook,
                        status=status, **kwargs)
    
    def __call__(self, state0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None):
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

    def set_step_fun(self, step_fun, step_args=None):
        """
        Set time step function and its optional arguments.
        """
        self.step_fun = step_fun
        self.step_args = get_default(step_args, ())

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
