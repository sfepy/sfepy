"""
Base (abstract) solver classes.
"""
import numpy as nm
from sfepy.base.base import Struct, output, structify

def make_get_conf(conf, kwargs):
    def _get_conf_item(name, default=None, msg_if_none=None):
        return kwargs.get(name, conf.get(name, default=default,
                                         msg_if_none=msg_if_none))

    return _get_conf_item

def format_next(text, new_text, pos, can_newline, width, ispaces):
    new_len = len(new_text)

    if (pos + new_len > width) and can_newline:
        text += '\n' + ispaces + new_text
        pos = new_len
        can_newline = False

    else:
        if pos > 0:
            text += ' ' + new_text
            pos += new_len + 1

        else:
            text += ispaces + new_text
            pos += len(ispaces) + new_len

        can_newline = True

    return text, pos, can_newline

def typeset_to_indent(txt, indent, width):
    if not len(txt): return txt

    txt_lines = txt.strip().split('\n')
    ispaces = ' ' * indent

    can_newline = False
    pos = 0
    text = ''
    for line in txt_lines:
        for word in line.split():
            text, pos, can_newline = format_next(text, word, pos, can_newline,
                                                 width, ispaces)

    return text

def make_option_docstring(name, kind, default, required, doc):
    if default is None:
        entry = '    %s : %s\n' % (name, kind)

    else:
        entry = '    %s : %s (default: %s)\n' % (name, kind, repr(default))

    entry += typeset_to_indent(doc, 8, 75)

    return entry


def use_first_available(solver_list, context=None, **kwargs):
    """
    Use the first available solver from `solver_list`.

    Parameters
    ----------
    solver_list : list of str or Struct
        The list of solver names or configuration objects.
    context : object, optional
        An optional solver context to pass to the solver.
    **kwargs : keyword arguments
        Additional solver options, see the particular __init__() methods.

    Returns
    -------
    out : Solver
        The first available solver.
    """
    prev_name = None
    for solver, conf in solver_list:
        name = solver if isinstance(solver, str) else solver.name
        if prev_name is not None:
            output("'%s' not available, trying '%s'" % (prev_name, name))

        try:
            _conf = conf.copy()
            _conf.__dict__.pop('fallback', None)
            if isinstance(solver, str):
                out = Solver.any_from_conf(_conf, context=context, **kwargs)
            else:
                out = solver(_conf, contex=context, **kwargs)

            output("using '%s' solver" % name)
            break

        except (KeyboardInterrupt, SystemExit):
            raise

        except Exception as exc:
            output('solver construction failed with: "{0!r}"'.format(exc))
            prev_name = name

    else:
        raise ValueError('no solver available from %s!'
                         % [ii[0] for ii in solver_list])

    return out

par_template = \
"""
    For common configuration parameters, see :class:`Solver
    <sfepy.solvers.solvers.Solver>`.

    Specific configuration parameters:

    Parameters
    ----------
"""

class SolverMeta(type):
    """
    Metaclass for solver classes that automatically adds configuration
    parameters to the solver class docstring from the ``_parameters`` class
    attribute.
    """

    def __new__(cls, name, bases, ndict):
        if '__doc__' in ndict:
            ndict['__doc__'] = ndict['__doc__'].rstrip()

            solver_kind = ndict.get('name')
            if solver_kind is not None:
                ndict['__doc__'] += '\n\n    Kind: %s\n' % repr(solver_kind)
                ndict['__doc__'] += par_template

            options = ndict.get('_parameters')
            if options is not None:
                docs = [make_option_docstring(*option) for option in options]
                ndict['__doc__'] = (ndict['__doc__'].rstrip()
                                    + '\n' + '\n'.join(docs))

        return super(SolverMeta, cls).__new__(cls, name, bases, ndict)

class Solver(Struct, metaclass=SolverMeta):
    """
    Base class for all solver kinds. Takes care of processing of common
    configuration options.

    The factory method any_from_conf() can be used to create an instance of any
    subclass.

    The subclasses have to reimplement __init__() and __call__().

    All solvers use the following configuration parameters:

    Parameters
    ----------
    """

    _parameters = [
        ('name', 'str', None, True,
         'The name referred to in problem description options.'),
        ('kind', 'str', None, True,
         """The solver kind, as given by the `name` class attribute of the
            Solver subclasses."""),
        ('report_status', 'bool', False, False,
         'If True, the solver can report the solution status.'),
        ('verbose', 'bool', False, False,
         """If True, the solver can print more information about the
            solution."""),
    ]

    @classmethod
    def process_conf(cls, conf, kwargs):
        """
        Process configuration parameters.
        """
        get = make_get_conf(conf, kwargs)

        if len(cls._parameters) and cls._parameters[0][0] != 'name':
            options = Solver._parameters + cls._parameters

        else:
            options = Solver._parameters

        opts = Struct()
        allow_extra = False
        for name, _, default, required, _ in options:
            if name == '*':
                allow_extra = True
                continue

            msg = ('missing "%s" in options!' % name) if required else None
            setattr(opts, name, get(name, default, msg))

        if allow_extra:
            all_keys = set(conf.to_dict().keys())
            other = all_keys.difference(list(opts.to_dict().keys()))
            for name in other:
                setattr(opts, name, get(name, None, None))

        return opts

    def __init__(self, conf=None, context=None, **kwargs):
        if conf is None:
            conf = Struct()

        elif isinstance(conf, dict):
            conf = Struct(**conf)

        else:
            conf = conf.copy()

        if conf.get('name', None) is None:
            conf.name = 'auto_' + self.__class__.__name__

        if conf.get('kind', None) is None:
            if hasattr(self.__class__, 'name'):
                conf.kind = self.__class__.name

            else:
                raise ValueError('solver kind cannot be determined!')

        new_conf = self.process_conf(conf, kwargs)
        Struct.__init__(self, conf=new_conf, orig_conf=conf, context=context,
                        **kwargs)

    def build_solver_kwargs(self, conf):
        """
        Build the `kwargs` dict for the underlying solver function using the
        extra options (marked by '*' in ``_parameters``) in `conf`. The
        declared parameters are omitted.
        """
        if len(self._parameters) and self._parameters[0][0] != 'name':
            options = Solver._parameters + self._parameters

        else:
            options = self._parameters

        std = set([ii[0] for ii in options if ii[0] != '*'])

        kwargs = {}
        for key, val in conf.to_dict().items():
            if key not in std:
                kwargs[key] = val

        return kwargs

    def __call__(self, **kwargs):
        raise ValueError('called an abstract Solver instance!')

    def set_field_split(self, field_ranges, **kwargs):
        pass

class LinearSolver(Solver):
    """
    Abstract linear solver class.
    """
    def __init__(self, conf, mtx=None, status=None, context=None, **kwargs):
        Solver.__init__(self, conf=conf, mtx=mtx, status=status,
                        context=context, **kwargs)
        self.mtx_digest = (0, '')

    def __call__(self, rhs, x0=None, conf=None, eps_a=None, eps_r=None,
                 i_max=None, mtx=None, status=None, context=None, **kwargs):
        raise ValueError('called an abstract LinearSolver instance!')

    def get_tolerance(self):
        """
        Return tuple `(eps_a, eps_r)` of absolute and relative tolerance
        settings. Either value can be `None`, meaning that the solver
        does not use that setting.
        """
        return self.conf.get('eps_a', None), self.conf.get('eps_r', None)

    def clear(self):
        pass

    def presolve(self, mtx):
        pass

class NonlinearSolver(Solver):
    """
    Abstract nonlinear solver class.
    """

    def __init__(self, conf, fun=None, fun_grad=None, lin_solver=None,
                 iter_hook=None, status=None, context=None, **kwargs):
        Solver.__init__(self, conf=conf, fun=fun, fun_grad=fun_grad,
                        lin_solver=lin_solver, iter_hook=iter_hook,
                        status=status, context=context, **kwargs)

    def __call__(self, state0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None):
        raise ValueError('called an abstract NonlinearSolver instance!')

class TimeStepController(Solver):
    """
    Abstract time step controller class.
    """

    def __init__(self, conf, **kwargs):
        Solver.__init__(self, conf=conf, **kwargs)

    def get_initial_dt(self, ts, vec, **kwargs):
        return ts.dt

    def __call__(self, ts, vec0, vec1, **kwargs):
        raise ValueError('called an abstract TimeStepController instance!')

class TimeSteppingSolver(Solver):
    """
    Abstract time stepping solver class.
    """

    def __init__(self, conf, nls=None, status=None, tsc=None, context=None,
                 **kwargs):
        Solver.__init__(self, conf=conf, nls=nls, tsc=tsc, status=status,
                        context=context, **kwargs)

    def set_dof_info(self, di):
        self.di = di

    def __call__(self, vec0=None, nls=None, tsc=None, init_fun=None,
                 prestep_fun=None, poststep_fun=None,
                 status=None, **kwargs):
        raise ValueError('called an abstract TimeSteppingSolver instance!')

class OptimizationSolver(Solver):
    """
    Abstract optimization solver class.
    """

    def __init__(self, conf, obj_fun=None, obj_fun_grad=None, status=None,
                 obj_args=None, context=None, **kwargs):
        Solver.__init__(self, conf=conf, obj_fun=obj_fun,
                        obj_fun_grad=obj_fun_grad, status=status,
                        obj_args=obj_args, context=context, **kwargs)

    def __call__(self, state0, conf=None, obj_fun=None, obj_fun_grad=None,
                 status=None, obj_args=None):
        raise ValueError('called an abstract OptimizationSolver instance!')

class EigenvalueSolver(Solver):
    """
    Abstract eigenvalue solver class.
    """

    def __init__(self, conf, mtx_a=None, mtx_b=None, n_eigs=None,
                 eigenvectors=None, status=None, context=None, **kwargs):
        Solver.__init__(self, conf=conf, mtx_a=mtx_a, mtx_b=mtx_b,
                        n_eigs=n_eigs, eigenvectors=eigenvectors,
                        status=status, context=context, **kwargs)

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

    def _ret_zero(self, mtx, eigenvectors=None):
        out = nm.zeros(0, dtype=nm.float64)
        if eigenvectors:
            out = (out,
                   nm.zeros((mtx.shape[0], 0), dtype=nm.float64))
        return out

class QuadraticEVPSolver(Solver):
    """
    Abstract quadratic eigenvalue problem solver class.
    """

    def __init__(self, conf, mtx_m=None, mtx_d=None, mtx_k=None, n_eigs=None,
                 eigenvectors=None, status=None, context=None, **kwargs):
        Solver.__init__(self, conf=conf, mtx_m=mtx_m, mtx_d=mtx_d, mtx_k=mtx_k,
                        n_eigs=n_eigs, eigenvectors=eigenvectors,
                        status=status, context=context)
        solver_conf = structify(self.conf.solver)
        self.solver = Solver.any_from_conf(solver_conf)

    def __call__(self, mtx_m, mtx_d, mtx_k, n_eigs=None,
                 eigenvectors=None, status=None, conf=None):
        raise ValueError('called an abstract QuadraticEVPSolver instance!')
