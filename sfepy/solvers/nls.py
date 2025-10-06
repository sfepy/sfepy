"""
Nonlinear solvers.
"""
from functools import partial

import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import output, get_default, Struct
from sfepy.base.log import Log, get_logging_conf
from sfepy.base.timing import Timer, Timers
from sfepy.solvers.solvers import NonlinearSolver

def standard_nls_call(call):
    """
    Decorator handling argument preparation and timing for nonlinear solvers.
    """

    def _standard_nls_call(self, vec_x0, conf=None, fun=None, fun_grad=None,
                           lin_solver=None, iter_hook=None, status=None,
                           **kwargs):
        timer = Timer(start=True)

        status = get_default(status, self.status)
        fun = get_default(fun, self.fun, 'function has to be specified!')
        result = call(self, vec_x0, conf=conf, fun=fun, fun_grad=fun_grad,
                      lin_solver=lin_solver, iter_hook=iter_hook,
                      status=status, **kwargs)

        elapsed = timer.stop()
        if status is not None:
            status['time'] = elapsed

        return result

    return _standard_nls_call

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

    timer = Timer(start=True)
    output(mtx_a, '.. analytical')
    output(mtx_d, '.. difference')
    import sfepy.base.plotutils as plu
    plu.plot_matrix_diff(mtx_d, mtx_a, delta, ['difference', 'analytical'],
                         conf.check)

    return timer.stop()

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

    conv_a = err <= conf.eps_a
    if it > 0:
        conv_r = err_r <= conf.eps_r

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

def apply_line_search_bt(vec_x0, vec_r0, vec_dx0, it, err_last, conf, fun,
                         apply_lin_solver, timers, log=None, context=None):
    """
    Apply a backtracking line-search on iteration residuals or error estimates
    depending on `conf.ls_mode` parameter. The error mode corresponds to the
    damped Newton method [1]_.

    .. [1] P. Deuflhard, A modified Newton method for the solution of
           ill-conditioned systems of nonlinear equations with application to
           multiple shooting, Numer. Math. 22, 289 (1974).
    """
    ls = 1.0
    vec_dx = vec_dx0

    while 1:
        vec_x = vec_x0 - vec_dx

        timers.residual.start()
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

        timers.residual.stop()

        if ok:
            err = nla.norm(vec_r)
            if not nm.isfinite(err):
                output('residual:', vec_r)
                output(nm.isfinite(vec_r).all())
                raise ValueError('infs or nans in the residual')

            if log is not None:
                log(err, it)

            if (it == 0):
                break

            if conf.ls_mode == 'residual':
                if (err < (err_last * conf.ls_on)):
                    break

            elif conf.ls_mode == 'error':
                vec_dx = apply_lin_solver(vec_r, resolve=True)

                edx0 = nla.norm(vec_dx0)
                edx = nla.norm(vec_dx)
                if edx < (edx0 * conf.ls_on):
                    break

            else:
                raise ValueError("ls_mode can be 'residual' or 'error',"
                                 f" not '{conf.ls_mode}'!")

            red = conf.ls_red
            output('linesearch: iter %d, (%.5e < %.5e) (new ls: %e)'
                   % (it, err, err_last * conf.ls_on, red * ls))

        else: # Failure.
            if conf.give_up_warp:
                output('giving up!')
                break

            red = conf.ls_red_warp
            output('residual computation failed for iter %d'
                   ' (new ls: %e)!' % (it, red * ls))

        if ls < conf.ls_min:
            output('linesearch failed, continuing anyway')
            break

        ls *= red
        vec_dx = ls * vec_dx0

    return vec_x, vec_r, err, ok

def apply_line_search_ainc(vec_x0, vec_r0, vec_dx0, it, err_last, conf, fun,
                           apply_lin_solver, timers, log=None, context=None):
    r"""
    Apply an initial step size given by Affine-Invariant Cubic Newton [1]_,
    then continue with :func:`apply_line_search_bt()`.

    .. [1] S. Hanzely, D. Kamzolov, D. Pasechnyuk, A. Gasnikov, P. Richtarik,
           and M. Takac, A Damped Newton Method Achieves Global :math:`\mathcal
           O \left(\frac{1}{k^2}\right)` and Local Quadratic Convergence Rate,
           Advances in Neural Information Processing Systems 35, 25320 (2022).
    """
    if it > 0:
        lest = conf.ls_red
        gnorm = lest * nm.sqrt(vec_r0 @ vec_dx0, dtype=nm.complex128)
        if gnorm.imag != 0.0:
            output('warning: matrix is not s.p.d!'
                   f' (r^T inv(A) r: {vec_r0 @ vec_dx0})')
            alpha = 1.0

        else:
            gnorm = gnorm.real

            if gnorm == 0.0:
                alpha = 0.0

            else:
                alpha = (-1.0 + nm.sqrt(1.0 + 2.0 * gnorm)) / gnorm

        vec_dx = alpha * vec_dx0

        output(f'alpha: {alpha:.16f}')

    else:
        vec_dx = vec_dx0

    out = apply_line_search_bt(vec_x0, vec_r0, vec_dx, it, err_last, conf, fun,
                               apply_lin_solver, timers, log=log,
                               context=context)
    return out

def find_contact_ipc_term(equations):
    for eq in equations:
        for term in eq.terms:
            if term.name == 'dw_contact_ipc':
                break

        else:
            continue

        break

    else:
        raise ValueError('no dw_contact_ipc in equations!')

    return term

def apply_line_search_ipc(vec_x0, vec_r0, vec_dx0, it, err_last, conf, fun,
                          apply_lin_solver, timers, log=None, context=None,
                          clog=None):
    """
    Apply a backtracking line-search with continuous collision detection from
    IPC toolkit.
    """
    pb = context
    term = find_contact_ipc_term(pb.equations)

    select_other = lambda term: term.name != 'dw_contact_ipc'

    ls = 1.0
    vec_dx = vec_dx0

    variables = pb.get_variables()

    ci = term.get_contact_info(variables['u'])
    ci.it = it

    ls_it = 0
    while 1:
        vec_x = vec_x0 - vec_dx
        if it > 0:
            # Determine max. step size using continuous collision detection.
            u0 = variables.make_full_vec(vec_x0).reshape((-1, ci.dim))
            u1 = variables.make_full_vec(vec_x).reshape((-1, ci.dim))
            x0 = ci.smesh.coors + u0[ci.nods]
            x1 = ci.smesh.coors + u1[ci.nods]

            max_step_size = term.ipc.compute_collision_free_stepsize(
                ci.collision_mesh, x0, x1,
            )
            if max_step_size < 1.0:
                vec_x = vec_x0 - max_step_size * vec_dx

            ls = min(ls, max_step_size)

        timers.residual.start()

        ci.ls_it = ls_it

        try:
            vec_r1 = fun(vec_x, select_term=select_other)

            vec_r1_full = pb.equations.make_full_vec(vec_r1, force_value=0.0)
            ci.e_grad_full = vec_r1_full

            val, iels = term.evaluate(mode='weak', dw_mode='vector',
                                      standalone=False, ci=ci)

            vec_r2 = nm.zeros_like(vec_r1)
            term.assemble_to(vec_r2, val, iels, mode='vector')

            vec_r = vec_r1 + vec_r2

        except ValueError:
            if (it == 0) or (ls < conf.ls_min):
                output('giving up!')
                raise

            else:
                ok = False

        else:
            ok = True

        timers.residual.stop()

        if ok:
            err = nm.linalg.norm(vec_r)
            if not nm.isfinite(err):
                output('residual:', vec_r)
                output(nm.isfinite(vec_r).all())
                raise ValueError('infs or nans in the residual')

            if log is not None:
                log(err, it)

            if clog is not None:
                if it > 0:
                    clog(ci.min_distance, ci.barrier_stiffness,
                         ci.bp_grad_norm, ci.e_grad_norm)

                else:
                    clog(nm.nan, nm.nan, nm.nan, nm.nan)

            if (it == 0) or (err < (err_last * conf.ls_on)):
                break

            red = conf.ls_red
            output('linesearch: iter %d, (%.5e < %.5e) (new ls: %e)'
                   % (it, err, err_last * conf.ls_on, red * ls))

        else: # Failure.
            if conf.give_up_warp:
                output('giving up!')
                break

            red = conf.ls_red_warp
            output('residual computation failed for iter %d'
                   ' (new ls: %e)!' % (it, red * ls))

        if ls < conf.ls_min:
            output('linesearch failed, continuing anyway')
            break

        ls *= red
        vec_dx = ls * vec_dx0
        ls_it += 1

    return vec_x, vec_r, err, ok

class Newton(NonlinearSolver):
    r"""
    Solves a nonlinear system :math:`f(x) = 0` using the Newton method.

    The solver uses a backtracking line-search on divergence.
    """
    name = 'nls.newton'

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
        ('is_new_jacobian_fun',
         'function(it, it_jac, rhs1, rhs0, x1, x0, context)', None, False,
         """If given, the Jacobian (tangent) matrix is updated only when it
            returns True. `it` is the current iteration, `it_jac` the last
            iteration when the Jacobian was updated, `rhs1`, `x1` are the
            current iteration residual, solution and `rhs0`, `x0` are from the
            previous iteration.
            Can be used to implement the modified Newton method."""),
        ('scale_system_fun', 'function(mtx, rhs, x0, context)', None, False,
         """User-defined function for scaling the linear system and initial
            guess in each iteration."""),
        ('scale_solution_fun', 'function(x, context)', None, False,
         'User-defined function for scaling the solution in each iteration.'),
        ('scaled_error', 'bool', False, False,
         """If True, the error of the linear solver is calculated using
            the scaled values."""),
        ('lin_red', 'float or None', 1.0, False,
         """The linear system solution error should be smaller than (`eps_a` *
            `lin_red`), otherwise a warning is printed. If None, the check is
            skipped."""),
        ('lin_precision', 'float or None', None, False,
         """If not None, the linear system solution tolerances are set in each
            nonlinear iteration relative to the current residual norm by the
            `lin_precision` factor. Ignored for direct linear solvers."""),
        ('step_red', '0.0 < float <= 1.0', 1.0, False,
         """Step reduction factor. Equivalent to the mixing parameter :math:`a`:
            :math:`(1 - a) x + a (x + dx) = x + a dx`"""),
        ('line_search_fun',
         'function(it, vec_x0, vec_r0, vec_dx0, err_last, conf, fun,'
         ' apply_lin_solver, timers, log=None)',
         apply_line_search_bt, False, 'The line search function.'),
        ('ls_mode', "'residual' or 'error'", 'residual', False,
         """The line search mode: when it is 'residual', the solver tries to
            make the iteration residuals decreasing while for 'error'
            the solution error estimates should decrease."""),
        ('ls_on', 'float', 0.99999, False,
         r"""Start the backtracking line-search by reducing the step, if
            :math:`||d(x^i)|| / ||d(x^{i-1})||` is larger than `ls_on`,
            where :math:`d` is either :math:`f` or :math:`\Delta x`
            depending on `ls_mode`."""),
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
        ('log_vlines', "'iteration' or 'solve' or both",
         ('solve',), False,
         """Put log vertical lines after each iteration and/or before the
            solve."""),
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

    def _apply_lin_solver(self, vec_r, mtx_a, vec_x, lin_solver, ls_status,
                          err, conf, timers, resolve=False):
        scale_system = get_default(conf.scale_system_fun,
                                   lambda mtx, rhs, x0, context: (mtx, rhs, x0))
        scale_solution = get_default(conf.scale_solution_fun,
                                   lambda x, context: x)
        scaled_error = get_default(conf.scaled_error, False)
        if conf.lin_red is not None:
            lin_red = conf.eps_a * conf.lin_red

        else:
            lin_red = None

        ls_eps_a, ls_eps_r = lin_solver.get_tolerance()
        eps_a = get_default(ls_eps_a, 1.0)
        eps_r = get_default(ls_eps_r, 1.0)

        if conf.lin_precision is not None:
            if ls_eps_a is not None:
                eps_a = max(err * conf.lin_precision, ls_eps_a)

            elif ls_eps_r is not None:
                eps_r = max(conf.lin_precision, ls_eps_r)

            if lin_red is not None:
                lin_red = max(eps_a, err * eps_r)

        if conf.verbose:
            aux = 'modified ' if resolve else ''
            output(f'solving {aux}linear system...')

        timers.solve.start()
        if not resolve:
            smtx_a, svec_r, svec_x = scale_system(mtx_a, vec_r, vec_x,
                                                  self.context)
            svec_dx = lin_solver(svec_r, x0=svec_x,
                                 eps_a=eps_a, eps_r=eps_r, mtx=smtx_a,
                                 status=ls_status)
            lin_solver.mtx = smtx_a

        else:
            _, svec_r, svec_x = scale_system(None, vec_r, vec_x,
                                             self.context)
            svec_dx = lin_solver(svec_r, x0=svec_x,
                                 eps_a=eps_a, eps_r=eps_r,
                                 status=ls_status)

        vec_dx = scale_solution(svec_dx, self.context)
        timers.solve.stop()

        if conf.verbose:
            output('...done')

        if lin_red is not None:
            if scaled_error:
                vec_e = smtx_a @ svec_dx - svec_r
            else:
                vec_e = mtx_a @ vec_dx - vec_r
            lerr = nla.norm(vec_e)
            if lerr > lin_red:
                output('warning: linear system solution precision is lower'
                       ' then the value set in solver options!'
                       ' (err = %e < %e)' % (lerr, lin_red))

        return vec_dx

    @standard_nls_call
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
        is_new_jacobian = get_default(conf.is_new_jacobian_fun,
                                      lambda *args, **kwargs: True)
        apply_line_search = get_default(conf.line_search_fun,
                                        apply_line_search_bt)

        timers = Timers(['residual', 'matrix', 'solve'])
        if conf.check:
            timers.create('check')

        vec_x = vec_x0.copy()
        vec_x_last = vec_x0.copy()
        vec_r_last = 0.0
        vec_dx = 0.0
        mtx_a = None

        if (self.log is not None) and ('solve' in conf.log_vlines):
            self.log.plot_vlines(color='r', linewidth=1.0)

        err = err0 = -1.0
        err_last = -1.0
        it = 0
        it_jac = 0
        ls_status = {}
        ls_n_iter = 0
        while 1:
            if iter_hook is not None:
                iter_hook(self.context, self, vec_x, it, err, err0)

            apply_lin_solver = partial(self._apply_lin_solver,
                                       mtx_a=mtx_a, vec_x=vec_x,
                                       lin_solver=lin_solver,
                                       ls_status=ls_status,
                                       err=err, conf=conf, timers=timers)
            vec_x, vec_r, err, ok = apply_line_search(
                vec_x_last, vec_r_last, vec_dx, it, err_last, conf, fun,
                apply_lin_solver, timers, log=self.log, context=self.context,
            )
            if it == 0:
                err0 = err

            if (self.log is not None) and ('iteration' in conf.log_vlines):
                self.log.plot_vlines([1], color='g', linewidth=0.5)

            condition = conv_test(conf, it, err, err0)
            if condition >= 0:
                break

            if (not ok) and conf.give_up_warp:
                condition = 2
                break

            timers.matrix.start()
            if not conf.is_linear:
                if (is_new_jacobian(it, it_jac, vec_r, vec_r_last, vec_x,
                                    vec_x_last, context=self.context)
                    or (mtx_a is None)):
                    mtx_a = fun_grad(vec_x)
                    it_jac = it

            else:
                mtx_a = fun_grad('linear')

            timers.matrix.stop()

            if conf.check:
                timers.check.start()
                wt = check_tangent_matrix(conf, vec_x, fun, fun_grad)
                timers.check.stop()
                timers.check.add(-wt)

            vec_dx = self._apply_lin_solver(
                vec_r, mtx_a, vec_x, lin_solver, ls_status, err, conf, timers,
            )
            ls_n_iter += ls_status['n_iter']

            for key, val in timers.get_dts().items():
                output('%10s: %7.2f [s]' % (key, val))

            vec_dx *= conf.step_red
            it += 1

            err_last = err
            vec_x_last = vec_x.copy()
            vec_r_last = vec_r.copy()

        time_stats = timers.get_totals()
        ls_n_iter = ls_n_iter if ls_n_iter >= 0 else -1
        if status is not None:
            status['time_stats'] = time_stats
            status['err0'] = err0
            status['err'] = err
            status['n_iter'] = it
            status['ls_n_iter'] = ls_n_iter
            status['condition'] = condition

        if conf.report_status:
            output(f'cond: {condition}, iter: {it}, ls_iter: {ls_n_iter},'
                   f' err0: {err0:.8e}, err: {err:.8e}')
            for key, val in time_stats.items():
                output('%8s: %.8f [s]' % (key, val))
            output('     sum: %.8f [s]' % sum(time_stats.values()))

        if conf.log.plot is not None:
            if self.log is not None:
                self.log(save_figure=conf.log.plot)

        return vec_x

class ScipyRoot(NonlinearSolver):
    """
    Interface to ``scipy.optimize.root()``.
    """
    name = 'nls.scipy_root'

    methods = ['hybr', 'lm', 'broyden1', 'broyden2', 'anderson', 'linearmixing',
               'diagbroyden', 'excitingmixing', 'krylov', 'df-sane']

    _parameters = [
        ('method', 'str', 'anderson', False,
         f"""Type of solver supported in ``scipy.optimize.root()``, one of:
             {methods}"""),
        ('use_jacobian', 'bool', False, False,
         'If True, use the exact Jacobian.'),
        ('tol', 'float', None, False,
         """Tolerance for termination. For detailed control, use solver-specific
            options."""),
        ('callback', 'callback(x, f)', None, False,
         """Optional callback function. It is called on every iteration as
            with x the current solution and f the corresponding residual.
            For all methods but ‘hybr’ and ‘lm’."""),
        ('options', 'dict', None, False,
         """A dictionary of solver options. E.g., `xtol` or `maxiter`, see
            ``scipy.optimize.show_options('root')`` for details."""),
    ]

    def __init__(self, conf, **kwargs):
        from scipy.optimize import root

        NonlinearSolver.__init__(self, conf, root=root, **kwargs)
        if self.conf.method not in self.methods:
            raise ValueError(f"'method' option must be on of {self.methods}!")

    @standard_nls_call
    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None):
        conf = get_default(conf, self.conf)
        fun = get_default(fun, self.fun)
        fun_grad = get_default(fun_grad, self.fun_grad)
        status = get_default(status, self.status)

        if self.conf.method not in self.methods:
            raise ValueError(f"'method' option must be on of {self.methods}!")

        timer = Timer(start=True)

        options = conf.options.copy() if conf.options is not None else {}
        options['disp'] = conf.verbose

        sol = self.root(fun, vec_x0,
                        method=conf.method,
                        jac=fun_grad if conf.use_jacobian else None,
                        tol=conf.tol,
                        callback=conf.callback,
                        options=options)
        err = nla.norm(sol.fun)

        nit = sol.get('nit', -1)
        if status is not None:
            status['time_stats'] = {'solver' : timer.stop()}
            status['err'] = err
            status['n_iter'] = nit
            status['n_fev'] = sol.nfev
            status['condition'] = sol.status

        if conf.report_status:
            output(sol.message)
            output(f'status: {sol.status}, iter: {nit},'
                   f' nfev: {sol.nfev}, err: {err:.8e}')
            output('solver: %.8f [s]' % status['time_stats']['solver'])

        return sol.x

class PETScNonlinearSolver(NonlinearSolver):
    """
    Interface to PETSc SNES (Scalable Nonlinear Equations Solvers).

    The solver supports parallel use with a given MPI communicator (see `comm`
    argument of :func:`PETScNonlinearSolver.__init__()`). Returns a (global)
    PETSc solution vector instead of a (local) numpy array, when given a PETSc
    initial guess vector.

    For parallel use, the `fun` and `fun_grad` callbacks should be provided by
    :class:`PETScParallelEvaluator
    <sfepy.parallel.evaluate.PETScParallelEvaluator>`.
    """
    name = 'nls.petsc'

    _parameters = [
        ('method', 'str', 'newtonls', False,
         'The SNES type.'),
        ('i_max', 'int', 10, False,
         'The maximum number of iterations.'),
        ('if_max', 'int', 100, False,
         'The maximum number of function evaluations.'),
        ('eps_a', 'float', 1e-10, False,
         'The absolute tolerance for the residual, i.e. :math:`||f(x^i)||`.'),
        ('eps_r', 'float', 1.0, False,
         """The relative tolerance for the residual, i.e. :math:`||f(x^i)|| /
            ||f(x^0)||`."""),
        ('eps_s', 'float', 0.0, False,
         r"""The convergence tolerance in terms of the norm of the change in
            the solution between steps,
            i.e. $||delta x|| < \epsilon_s ||x||$"""),
    ]

    def __init__(self, conf, pmtx=None, prhs=None, comm=None, **kwargs):
        if comm is None:
            try:
                import petsc4py
                petsc4py.init([])
            except ImportError:
                msg = 'cannot import petsc4py!'
                raise ImportError(msg)

        from petsc4py import PETSc as petsc

        converged_reasons = {}
        for key, val in petsc.SNES.ConvergedReason.__dict__.items():
            if isinstance(val, int):
                converged_reasons[val] = key

        ksp_converged_reasons = {}
        for key, val in petsc.KSP.ConvergedReason.__dict__.items():
            if isinstance(val, int):
                ksp_converged_reasons[val] = key

        NonlinearSolver.__init__(self, conf, petsc=petsc,
                                 pmtx=pmtx, prhs=prhs, comm=comm,
                                 converged_reasons=converged_reasons,
                                 ksp_converged_reasons=ksp_converged_reasons,
                                 **kwargs)

    @standard_nls_call
    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None,
                 pmtx=None, prhs=None, comm=None):
        conf = self.conf
        fun = get_default(fun, self.fun)
        fun_grad = get_default(fun_grad, self.fun_grad)
        lin_solver = get_default(lin_solver, self.lin_solver)
        iter_hook = get_default(iter_hook, self.iter_hook)
        status = get_default(status, self.status)
        pmtx = get_default(pmtx, self.pmtx)
        prhs = get_default(prhs, self.prhs)
        comm = get_default(comm, self.comm)

        timer = Timer(start=True)

        if isinstance(vec_x0, self.petsc.Vec):
            psol = vec_x0

        else:
            psol = pmtx.getVecLeft()
            psol[...] = vec_x0

        snes = self.petsc.SNES()
        snes.create(comm)
        snes.setType(conf.method)

        ksp = lin_solver.create_ksp()
        snes.setKSP(ksp)
        ls_conf = lin_solver.conf
        ksp.setTolerances(atol=ls_conf.eps_a, rtol=ls_conf.eps_r,
                          divtol=ls_conf.eps_d, max_it=ls_conf.i_max)

        snes.setFunction(fun, prhs)
        snes.setJacobian(fun_grad, pmtx)

        snes.setTolerances(atol=conf.eps_a, rtol=conf.eps_r,
                           stol=conf.eps_s, max_it=conf.i_max)
        snes.setMaxFunctionEvaluations(conf.if_max)
        snes.setFromOptions()

        fun(snes, psol, prhs)
        err0 = prhs.norm()

        snes.solve(prhs.duplicate(), psol)

        if status is not None:
            status['time_stats'] = {'solver' : timer.stop()}

        if snes.reason in self.converged_reasons:
            reason = 'snes: %s' % self.converged_reasons[snes.reason]

        else:
            reason = 'ksp: %s' % self.ksp_converged_reasons[snes.reason]

        output('%s(%s): %d iterations in the last step'
               % (ksp.getType(), ksp.getPC().getType(),
                  ksp.getIterationNumber()),
               verbose=conf.verbose)

        output('%s convergence: %s (%s, %d iterations, %d function evaluations)'
               % (snes.getType(), snes.reason, reason,
                  snes.getIterationNumber(), snes.getFunctionEvaluations()),
               verbose=conf.verbose)

        converged = snes.reason >= 0
        condition = 0 if converged else -1
        n_iter = snes.getLinearSolveIterations()
        ls_n_iter = snes.getLinearSolveIterations()

        if not converged:
            # PETSc does not update the solution if KSP have not converged.
            dpsol = snes.getSolutionUpdate()
            psol -= dpsol

            fun(snes, psol, prhs)
            err = prhs.norm()

        else:
            try:
                err = snes.getFunctionNorm()

            except AttributeError:
                fun(snes, psol, prhs)
                err = prhs.norm()

        if status is not None:
            status['err0'] = err0
            status['err'] = err
            status['n_iter'] = n_iter
            status['ls_n_iter'] = ls_n_iter
            status['condition'] = condition

        if conf.report_status:
            output(f'cond: {condition}, iter: {n_iter}, ls_iter: {ls_n_iter},'
                   f' err0: {err0:.8e}, err: {err:.8e}')
            output('solver: %.8f [s]' % status['time_stats']['solver'])

        if isinstance(vec_x0, self.petsc.Vec):
            sol = psol

        else:
            sol = psol[...].copy()

        return sol
