from __future__ import absolute_import

import numpy as nm
import numpy.linalg as nla
import scipy.sparse as sp

from sfepy.base.base import output, get_default, debug
from sfepy.base.timing import Timer
from sfepy.solvers.nls import Newton, conv_test
from sfepy.linalg import compose_sparse
import six
from six.moves import range

class SemismoothNewton(Newton):
    r"""
    The semi-smooth Newton method.

    This method is suitable for solving problems of the following structure:

    .. math::
        \begin{split}
          & F(y) = 0 \\
          & A(y) \ge 0 \;,\ B(y) \ge 0 \;,\ \langle A(y), B(y) \rangle = 0
        \end{split}

    The function :math:`F(y)` represents the smooth part of the problem.

    Regular step: :math:`y \leftarrow y - J(y)^{-1} \Phi(y)`

    Steepest descent step: :math:`y \leftarrow y - \beta J(y) \Phi(y)`

    Although ``fun_smooth_grad()`` computes the gradient of the smooth part
    only, it should return the global matrix, where the non-smooth part is
    uninitialized, but pre-allocated.
    """
    name = 'nls.semismooth_newton'

    _parameters = [
        ('semismooth', 'bool', True, False,
         """If True, use the semi-smooth algorithm. Otherwise a non-smooth
            equation is assumed (use a brute force)."""),
        ('i_max', 'int', 1, False,
         'The maximum number of iterations.'),
        ('eps_a', 'float', 1e-10, False,
         'The absolute tolerance for the residual, i.e. :math:`||f(x^i)||`.'),
        ('eps_r', 'float', 1.0, False,
         """The relative tolerance for the residual, i.e. :math:`||f(x^i)|| /
            ||f(x^0)||`."""),
        ('macheps', 'float', nm.finfo(nm.float64).eps, False,
         'The float considered to be machine "zero".'),
        ('lin_red', 'float', 1.0, False,
         """The linear system solution error should be smaller than (`eps_a` *
            `lin_red`), otherwise a warning is printed."""),
        ('ls_on', 'float', 0.99999, False,
         """Start the backtracking line-search by reducing the step, if
            :math:`||f(x^i)|| / ||f(x^{i-1})||` is larger than `ls_on`."""),
        ('ls_red', 'dict', {'regular' : 0.1, 'steepest_descent' : 0.01}, False,
         """The step reduction factor in case of correct residual assembling
            for regular and steepest descent modes."""),
        ('ls_red_warp', '0.0 < float < 1.0', 0.001, False,
         """The step reduction factor in case of failed residual assembling
            (e.g. the "warp violation" error caused by a negative volume
            element resulting from too large deformations)."""),
        ('ls_min', '0.0 < float < 1.0', 1e-5, False,
         'The minimum step reduction factor.'),
    ]

    _colors = {'regular' : 'g', 'steepest_descent' : 'k'}

    def __call__(self, vec_x0, conf=None, fun_smooth=None, fun_smooth_grad=None,
                 fun_a=None, fun_a_grad=None, fun_b=None, fun_b_grad=None,
                 lin_solver=None, status=None):

        conf = get_default(conf, self.conf)

        fun_smooth = get_default(fun_smooth, self.fun_smooth)
        fun_smooth_grad = get_default(fun_smooth_grad, self.fun_smooth_grad)
        fun_a = get_default(fun_a, self.fun_a)
        fun_a_grad = get_default(fun_a_grad, self.fun_a_grad)
        fun_b = get_default(fun_b, self.fun_b)
        fun_b_grad = get_default(fun_b_grad, self.fun_b_grad)

        lin_solver = get_default(lin_solver, self.lin_solver)
        status = get_default(status, self.status)

        timer = Timer()
        time_stats = {}

        vec_x = vec_x0.copy()
        vec_x_last = vec_x0.copy()
        vec_dx = None

        if self.log is not None:
            self.log.plot_vlines(color='r', linewidth=1.0)

        err0 = -1.0
        err_last = -1.0
        it = 0
        step_mode = 'regular'
        r_last = None
        reuse_matrix = False
        while 1:

            ls = 1.0
            vec_dx0 = vec_dx;
            i_ls = 0
            while 1:
                timer.start()

                try:
                    vec_smooth_r = fun_smooth(vec_x)
                    vec_a_r = fun_a(vec_x)
                    vec_b_r = fun_b(vec_x)

                except ValueError:
                    vec_smooth_r = vec_semismooth_r = None
                    if (it == 0) or (ls < conf.ls_min):
                        output('giving up!')
                        raise

                    else:
                        ok = False

                else:
                    if conf.semismooth:
                        # Semi-smooth equation.
                        vec_semismooth_r = (nm.sqrt(vec_a_r**2.0 + vec_b_r**2.0)
                                            - (vec_a_r + vec_b_r))

                    else:
                        # Non-smooth equation (brute force).
                        vec_semismooth_r = nm.where(vec_a_r < vec_b_r,
                                                    vec_a_r, vec_b_r)

                    r_last = (vec_smooth_r, vec_a_r, vec_b_r, vec_semismooth_r)

                    ok = True

                time_stats['residual'] = timer.stop()

                if ok:
                    vec_r = nm.r_[vec_smooth_r, vec_semismooth_r]

                    try:
                        err = nla.norm(vec_r)
                    except:
                        output('infs or nans in the residual:',
                               vec_semismooth_r)
                        output(nm.isfinite(vec_semismooth_r).all())
                        debug()

                    if self.log is not None:
                        self.log(err, it)

                    if it == 0:
                        err0 = err;
                        break

                    if err < (err_last * conf.ls_on):
                        step_mode = 'regular'
                        break

                    else:
                        output('%s step line search' % step_mode)

                        red = conf.ls_red[step_mode];
                        output('iter %d, (%.5e < %.5e) (new ls: %e)'\
                               % (it, err, err_last * conf.ls_on, red * ls))

                else: # Failed to compute residual.
                    red = conf.ls_red_warp;
                    output('residual computation failed for iter %d'
                           ' (new ls: %e)!' % (it, red * ls))

                if ls < conf.ls_min:
                    if step_mode == 'regular':
                        output('restore previous state')
                        vec_x = vec_x_last.copy()
                        (vec_smooth_r, vec_a_r, vec_b_r,
                         vec_semismooth_r) = r_last
                        err = err_last
                        reuse_matrix = True

                        step_mode = 'steepest_descent'

                    else:
                        output('linesearch failed, continuing anyway')

                    break

                ls *= red;

                vec_dx = ls * vec_dx0;
                vec_x = vec_x_last.copy() - vec_dx

                i_ls += 1

            # End residual loop.

            output('%s step' % step_mode)

            if self.log is not None:
                self.log.plot_vlines([1],
                                     color=self._colors[step_mode],
                                     linewidth=0.5)

            err_last = err;
            vec_x_last = vec_x.copy()

            condition = conv_test(conf, it, err, err0)
            if condition >= 0:
                break

            timer.start()

            if not reuse_matrix:
                mtx_jac = self.compute_jacobian(vec_x, fun_smooth_grad,
                                                fun_a_grad, fun_b_grad,
                                                vec_smooth_r,
                                                vec_a_r, vec_b_r)

            else:
                reuse_matrix = False

            time_stats['matrix'] = timer.stop()

            timer.start()

            if step_mode == 'regular':
                vec_dx = lin_solver(vec_r, mtx=mtx_jac)

                vec_e = mtx_jac * vec_dx - vec_r
                lerr = nla.norm(vec_e)
                if lerr > (conf.eps_a * conf.lin_red):
                    output('linear system not solved! (err = %e)' % lerr)

                    output('switching to steepest descent step')
                    step_mode = 'steepest_descent'
                    vec_dx = mtx_jac.T * vec_r

            else:
                vec_dx = mtx_jac.T * vec_r

            time_stats['solve'] = timer.stop()

            for kv in six.iteritems(time_stats):
                output('%10s: %7.2f [s]' % kv)

            vec_x -= vec_dx
            it += 1

        if status is not None:
            status['time_stats'] = time_stats
            status['err0'] = err0
            status['err'] = err
            status['condition'] = condition

        if conf.log.plot is not None:
            if self.log is not None:
                self.log(save_figure=conf.log.plot)

        return vec_x

    def compute_jacobian(self, vec_x, fun_smooth_grad, fun_a_grad, fun_b_grad,
                         vec_smooth_r, vec_a_r, vec_b_r):
        conf = self.conf

        mtx_s = fun_smooth_grad(vec_x)
        mtx_a = fun_a_grad(vec_x)
        mtx_b = fun_b_grad(vec_x)

        n_s = vec_smooth_r.shape[0]
        n_ns = vec_a_r.shape[0]

        if conf.semismooth:
            aa = nm.abs(vec_a_r)
            ab = nm.abs(vec_b_r)
            iz = nm.where((aa < (conf.macheps * max(aa.max(), 1.0)))
                          & (ab < (conf.macheps * max(ab.max(), 1.0))))[0]
            inz = nm.setdiff1d(nm.arange(n_ns), iz)

            output('non_active/active: %d/%d' % (len(inz), len(iz)))

            mul_a = nm.empty_like(vec_a_r)
            mul_b = nm.empty_like(mul_a)

            # Non-active part of the jacobian.
            if len(inz) > 0:
                a_r_nz = vec_a_r[inz]
                b_r_nz = vec_b_r[inz]

                sqrt_ab = nm.sqrt(a_r_nz**2.0 + b_r_nz**2.0)
                mul_a[inz] = (a_r_nz / sqrt_ab) - 1.0
                mul_b[inz] = (b_r_nz / sqrt_ab) - 1.0

            # Active part of the jacobian.
            if len(iz) > 0:
                vec_z = nm.zeros_like(vec_x)
                vec_z[n_s+iz] = 1.0

                mtx_a_z = mtx_a[iz]
                mtx_b_z = mtx_b[iz]

                sqrt_ab = nm.empty((iz.shape[0],), dtype=vec_a_r.dtype)
                for ir in range(len(iz)):
                    row_a_z = mtx_a_z[ir]
                    row_b_z = mtx_b_z[ir]
                    sqrt_ab[ir] = nm.sqrt((row_a_z * row_a_z.T).todense()
                                          + (row_b_z * row_b_z.T).todense())
                mul_a[iz] = ((mtx_a_z * vec_z) / sqrt_ab) - 1.0
                mul_b[iz] = ((mtx_b_z * vec_z) / sqrt_ab) - 1.0

        else:
            iz = nm.where(vec_a_r > vec_b_r)[0]
            mul_a = nm.zeros_like(vec_a_r)
            mul_b = nm.ones_like(mul_a)

            mul_a[iz] = 1.0
            mul_b[iz] = 0.0

        mtx_ns = sp.spdiags(mul_a, 0, n_ns, n_ns) * mtx_a \
                 + sp.spdiags(mul_b, 0, n_ns, n_ns) * mtx_b

        mtx_jac = compose_sparse([[mtx_s], [mtx_ns]]).tocsr()
        mtx_jac.sort_indices()

        return mtx_jac
