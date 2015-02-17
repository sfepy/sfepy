"""
Nonlinear solvers.
"""
import time

import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import output, get_default, debug, Struct
from sfepy.base.log import Log, get_logging_conf
from sfepy.solvers.solvers import SolverMeta, NonlinearSolver

def check_tangent_matrix(conf, vec_x0, fun, fun_grad):
    """
    Verify the correctness of the tangent matrix as computed by `fun_grad()` by
    comparing it with its finite difference approximation evaluated by
    repeatedly calling `fun()` with `vec_x0` items perturbed by a small delta.
    """
    vec_x = vec_x0.copy()
    delta = conf.delta

    vec_r = fun(vec_x) # Update state.
    mtx_a0 = fun_grad(vec_x)

    mtx_a = mtx_a0.tocsc()
    mtx_d = mtx_a.copy()
    mtx_d.data[:] = 0.0

    vec_dx = nm.zeros_like(vec_r)

    for ic in range(vec_dx.shape[0]):
        vec_dx[ic] = delta
        xx = vec_x.copy() - vec_dx
        vec_r1 = fun(xx)

        vec_dx[ic] = -delta
        xx = vec_x.copy() - vec_dx
        vec_r2 = fun(xx)

        vec_dx[ic] = 0.0;

        vec = 0.5 * (vec_r2 - vec_r1) / delta

        ir = mtx_a.indices[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]]
        mtx_d.data[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]] = vec[ir]

    vec_r = fun(vec_x) # Restore.

    tt = time.clock()
    output(mtx_a, '.. analytical')
    output(mtx_d, '.. difference')
    import sfepy.base.plotutils as plu
    plu.plot_matrix_diff(mtx_d, mtx_a, delta, ['difference', 'analytical'],
                         conf.check)

    return time.clock() - tt

def conv_test(conf, it, err, err0):
    """
    Nonlinear solver convergence test.

    Parameters
    ----------
    conf : Struct instance
        The nonlinear solver configuration.
    it : int
        The current iteration.
    err : float
        The current iteration error.
    err0 : float
        The initial error.

    Returns
    -------
    status : int
        The convergence status: -1 = no convergence (yet), 0 = solver converged
        - tolerances were met, 1 = max. number of iterations reached.
    """
    status = -1
    if (abs(err0) < conf.macheps):
        err_r = 0.0
    else:
        err_r = err / err0

    output('nls: iter: %d, residual: %e (rel: %e)' % (it, err, err_r))

    conv_a = err < conf.eps_a
    if it > 0:
        conv_r = err_r < conf.eps_r

        if conv_a and conv_r:
            status = 0

        elif (conf.get('eps_mode', '') == 'or') and (conv_a or conv_r):
            status = 0

    else:
        if conv_a:
            status = 0

    if (status == -1) and (it >= conf.i_max):
        status = 1

    return status

class Newton(NonlinearSolver):
    r"""
    Solves a nonlinear system :math:`f(x) = 0` using the Newton method with
    backtracking line-search, starting with an initial guess :math:`x^0`.
    """
    name = 'nls.newton'

    __metaclass__ = SolverMeta

    _parameters = [
        ('i_max', 'int', 1, False,
         'The maximum number of iterations.'),
        ('eps_a', 'float', 1e-10, False,
         'The absolute tolerance for the residual, i.e. :math:`||f(x^i)||`.'),
        ('eps_r', 'float', 1.0, False,
         """The relative tolerance for the residual, i.e. :math:`||f(x^i)|| /
            ||f(x^0)||`."""),
        ('eps_mode', "'and' or 'or'", 'and', False,
         """The logical operator to use for combining the absolute and relative
            tolerances."""),
        ('macheps', 'float', nm.finfo(nm.float64).eps, False,
         'The float considered to be machine "zero".'),
        ('lin_red', 'float', 1.0, False,
         """The linear system solution error should be smaller than (`eps_a` *
            `lin_red`), otherwise a warning is printed."""),
        ('lin_precision', 'float or None', None, False,
         """If not None, the linear system solution tolerances are set in each
            nonlinear iteration relative to the current residual norm by the
            `lin_precision` factor. Ignored for direct linear solvers."""),
        ('ls_on', 'float', 0.99999, False,
         """Start the backtracking line-search by reducing the step, if
            :math:`||f(x^i)|| / ||f(x^{i-1})||` is larger than `ls_on`."""),
        ('ls_red', '0.0 < float < 1.0', 0.1, False,
         'The step reduction factor in case of correct residual assembling.'),
        ('ls_red_warp', '0.0 < float < 1.0', 0.001, False,
         """The step reduction factor in case of failed residual assembling
            (e.g. the "warp violation" error caused by a negative volume
            element resulting from too large deformations)."""),
        ('ls_min', '0.0 < float < 1.0', 1e-5, False,
         'The minimum step reduction factor.'),
        ('give_up_warp', 'bool', False, False,
         'If True, abort on the "warp violation" error.'),
        ('check', '0, 1 or 2', 0, False,
         """If >= 1, check the tangent matrix using finite differences.  If 2,
            plot the resulting sparsity patterns."""),
        ('delta', 'float', 1e-6, False,
         r"""If `check >= 1`, the finite difference matrix is taken as
            :math:`A_{ij} = \frac{f_i(x_j + \delta) - f_i(x_j - \delta)}{2
            \delta}`."""),
        ('log', 'dict or None', None, False,
         """If not None, log the convergence according to the configuration in
            the following form: ``{'text' : 'log.txt', 'plot' : 'log.pdf'}``.
            Each of the dict items can be None."""),
        ('is_linear', 'bool', False, False,
         'If True, the problem is considered to be linear.'),
    ]

    def __init__(self, conf, **kwargs):
        NonlinearSolver.__init__(self, conf, **kwargs)

        conf = self.conf

        log = get_logging_conf(conf)
        conf.log = log = Struct(name='log_conf', **log)
        conf.is_any_log = (log.text is not None) or (log.plot is not None)

        if conf.is_any_log:
            self.log = Log([[r'$||r||$'], ['iteration']],
                           xlabels=['', 'all iterations'],
                           ylabels=[r'$||r||$', 'iteration'],
                           yscales=['log', 'linear'],
                           is_plot=conf.log.plot is not None,
                           log_filename=conf.log.text,
                           formats=[['%.8e'], ['%d']])

        else:
            self.log = None

    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None):
        """
        Nonlinear system solver call.

        Solves a nonlinear system :math:`f(x) = 0` using the Newton method with
        backtracking line-search, starting with an initial guess :math:`x^0`.

        Parameters
        ----------
        vec_x0 : array
            The initial guess vector :math:`x_0`.
        conf : Struct instance, optional
            The solver configuration parameters,
        fun : function, optional
            The function :math:`f(x)` whose zero is sought - the residual.
        fun_grad : function, optional
            The gradient of :math:`f(x)` - the tangent matrix.
        lin_solver : LinearSolver instance, optional
            The linear solver for each nonlinear iteration.
        iter_hook : function, optional
            User-supplied function to call before each iteration.
        status : dict-like, optional
            The user-supplied object to hold convergence statistics.

        Notes
        -----
        * The optional parameters except `iter_hook` and `status` need
          to be given either here or upon `Newton` construction.
        * Setting `conf.is_linear == True` means a pre-assembled and possibly
          pre-solved matrix. This is mostly useful for linear time-dependent
          problems.
        """
        conf = get_default(conf, self.conf)
        fun = get_default(fun, self.fun)
        fun_grad = get_default(fun_grad, self.fun_grad)
        lin_solver = get_default(lin_solver, self.lin_solver)
        iter_hook = get_default(iter_hook, self.iter_hook)
        status = get_default(status, self.status)

        ls_eps_a, ls_eps_r = lin_solver.get_tolerance()
        eps_a = get_default(ls_eps_a, 1.0)
        eps_r = get_default(ls_eps_r, 1.0)
        lin_red = conf.eps_a * conf.lin_red

        time_stats = {}

        vec_x = vec_x0.copy()
        vec_x_last = vec_x0.copy()
        vec_dx = None

        if self.log is not None:
            self.log.plot_vlines(color='r', linewidth=1.0)

        err = err0 = -1.0
        err_last = -1.0
        it = 0
        while 1:
            if iter_hook is not None:
                iter_hook(self, vec_x, it, err, err0)

            ls = 1.0
            vec_dx0 = vec_dx;
            while 1:
                tt = time.clock()

                try:
                    vec_r = fun(vec_x)

                except ValueError:
                    if (it == 0) or (ls < conf.ls_min):
                        output('giving up!')
                        raise

                    else:
                        ok = False

                else:
                    ok = True

                time_stats['rezidual'] = time.clock() - tt
                if ok:
                    try:
                        err = nla.norm(vec_r)
                    except:
                        output('infs or nans in the residual:', vec_r)
                        output(nm.isfinite(vec_r).all())
                        debug()

                    if self.log is not None:
                        self.log(err, it)

                    if it == 0:
                        err0 = err;
                        break
                    if err < (err_last * conf.ls_on): break
                    red = conf.ls_red;
                    output('linesearch: iter %d, (%.5e < %.5e) (new ls: %e)'
                           % (it, err, err_last * conf.ls_on, red * ls))
                else: # Failure.
                    if conf.give_up_warp:
                        output('giving up!')
                        break

                    red = conf.ls_red_warp;
                    output('rezidual computation failed for iter %d'
                           ' (new ls: %e)!' % (it, red * ls))

                if ls < conf.ls_min:
                    output('linesearch failed, continuing anyway')
                    break

                ls *= red;

                vec_dx = ls * vec_dx0;
                vec_x = vec_x_last.copy() - vec_dx
            # End residual loop.

            if self.log is not None:
                self.log.plot_vlines([1], color='g', linewidth=0.5)

            err_last = err;
            vec_x_last = vec_x.copy()

            condition = conv_test(conf, it, err, err0)
            if condition >= 0:
                break

            if (not ok) and conf.give_up_warp:
                condition = 2
                break

            tt = time.clock()
            if not conf.is_linear:
                mtx_a = fun_grad(vec_x)

            else:
                mtx_a = fun_grad('linear')

            time_stats['matrix'] = time.clock() - tt

            if conf.check:
                tt = time.clock()
                wt = check_tangent_matrix(conf, vec_x, fun, fun_grad)
                time_stats['check'] = time.clock() - tt - wt

            if conf.lin_precision is not None:
                if ls_eps_a is not None:
                    eps_a = max(err * conf.lin_precision, ls_eps_a)

                elif ls_eps_r is not None:
                    eps_r = max(conf.lin_precision, ls_eps_r)

                lin_red = max(eps_a, err * eps_r)

            if conf.verbose:
                output('solving linear system...')

            tt = time.clock()
            vec_dx = lin_solver(vec_r, x0=vec_x,
                                eps_a=eps_a, eps_r=eps_r, mtx=mtx_a)
            time_stats['solve'] = time.clock() - tt

            if conf.verbose:
                output('...done')

            for kv in time_stats.iteritems():
                output('%10s: %7.2f [s]' % kv)

            vec_e = mtx_a * vec_dx - vec_r
            lerr = nla.norm(vec_e)
            if lerr > lin_red:
                output('warning: linear system solution precision is lower')
                output('then the value set in solver options! (err = %e < %e)'
                       % (lerr, lin_red))

            vec_x -= vec_dx
            it += 1

        if status is not None:
            status['time_stats'] = time_stats
            status['err0'] = err0
            status['err'] = err
            status['n_iter'] = it
            status['condition'] = condition

        if conf.log.plot is not None:
            if self.log is not None:
                self.log(save_figure=conf.log.plot)

        return vec_x

class ScipyBroyden(NonlinearSolver):
    """
    Interface to Broyden and Anderson solvers from ``scipy.optimize``.
    """
    name = 'nls.scipy_broyden_like'

    __metaclass__ = SolverMeta

    _parameters = [
        ('method', 'str', 'anderson', False,
         'The name of the solver in ``scipy.optimize``.'),
        ('i_max', 'int', 10, False,
         'The maximum number of iterations.'),
        ('alpha', 'float', 0.9, False,
         'See ``scipy.optimize``.'),
        ('M', 'float', 5, False,
         'See ``scipy.optimize``.'),
        ('f_tol', 'float', 1e-6, False,
         'See ``scipy.optimize``.'),
        ('w0', 'float', 0.1, False,
         'See ``scipy.optimize``.'),
    ]

    def __init__(self, conf, **kwargs):
        NonlinearSolver.__init__(self, conf, **kwargs)
        self.set_method(self.conf)

    def set_method(self, conf):
        import scipy.optimize as so

        try:
            solver = getattr(so, conf.method)
        except AttributeError:
            output('scipy solver %s does not exist!' % conf.method)
            output('using broyden3 instead')
            solver = so.broyden3
        self.solver = solver

    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None):
        if conf is not None:
            self.set_method(conf)
        else:
            conf = self.conf
        fun = get_default(fun, self.fun)
        status = get_default(status, self.status)

        tt = time.clock()

        kwargs = {'iter' : conf.i_max,
                  'alpha' : conf.alpha,
                  'verbose' : conf.verbose}

        if conf.method == 'broyden_generalized':
            kwargs.update({'M' : conf.M})

        elif conf.method in ['anderson', 'anderson2']:
            kwargs.update({'M' : conf.M, 'w0' : conf.w0})

        if conf.method in ['anderson', 'anderson2',
                           'broyden', 'broyden2' , 'newton_krylov']:
            kwargs.update({'f_tol' : conf.f_tol })

        vec_x = self.solver(fun, vec_x0, **kwargs)
        vec_x = nm.asarray(vec_x)

        if status is not None:
            status['time_stats'] = time.clock() - tt

        return vec_x
