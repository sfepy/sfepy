
import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import output, get_default, pause, Struct
from sfepy.base.log import Log, get_logging_conf
from sfepy.base.timing import Timer
from sfepy.solvers.solvers import OptimizationSolver

import scipy.optimize as sopt
import scipy.optimize.linesearch as linesearch

def conv_test(conf, it, of, of0, ofg_norm=None):
    """
    Returns
    -------
    flag : int
        * -1 ... continue
        *  0 ... small OF -> stop
        *  1 ... i_max reached -> stop
        *  2 ... small OFG -> stop
        *  3 ... small relative decrase of OF
     """

    status = -1
    output('opt: iter: %d, of: %e (||ofg||: %e)' % (it, of, ofg_norm))

    if (abs(of) < conf.eps_of):
        status = 0
    elif ofg_norm and (ofg_norm < conf.eps_ofg):
        status = 2
    elif (it > 0) and (abs(of0 - of) < (conf.eps_rd * abs(of0))):
        status = 3

    if (status == -1) and (it >= conf.i_max):
        status = 1

    return status

def wrap_function(function, args):
    ncalls = [0]
    times = []
    timer = Timer()
    def function_wrapper(x):
        ncalls[0] += 1
        timer.start()
        out = function(x, *args)
        times.append(timer.stop())
        return out
    return ncalls, times, function_wrapper

def check_gradient(xit, aofg, fn_of, delta, check):

    dofg = nm.zeros_like(aofg)
    xd = xit.copy()
    for ii in range(xit.shape[0]):
        xd[ii] = xit[ii] + delta
        ofp = fn_of(xd)

        xd[ii] = xit[ii] - delta
        ofm = fn_of(xd)

        xd[ii] = xit[ii]

        dofg[ii] = 0.5 * (ofp - ofm) / delta

        output('**********', ii, aofg[ii], dofg[ii])

    diff = abs(aofg - dofg)
    aux = nm.concatenate((aofg[:,nm.newaxis], dofg[:,nm.newaxis],
                          diff[:,nm.newaxis]), 1)
    output(aux)
    output(nla.norm(diff, nm.inf))
    aofg.tofile('aofg.txt', ' ')
    dofg.tofile('dofg.txt', ' ')
    diff.tofile('diff.txt', ' ')
    if check == 2:
        import pylab
        pylab.plot(aofg)
        pylab.plot(dofg)
        pylab.legend(('analytical', 'finite difference'))
        pylab.show()
    pause('gradient checking done')

class FMinSteepestDescent(OptimizationSolver):
    """
    Steepest descent optimization solver.
    """
    name = 'opt.fmin_sd'

    _parameters = [
        ('i_max', 'int', 10, False,
         'The maximum number of iterations.'),
        ('eps_rd', 'float', 1e-5, False,
         'The relative delta of the objective function.'),
        ('eps_of', 'float', 1e-4, False,
         'The tolerance for the objective function.'),
        ('eps_ofg', 'float', 1e-8, False,
         'The tolerance for the objective function gradient.'),
        ('norm', 'numpy norm', nm.inf, False,
         'The norm to be used.'),
        ('ls', 'bool', True, False,
         'If True, use a line-search.'),
        ('ls_method', "{'backtracking', 'full'}", 'backtracking', False,
         'The line-search method.'),
        ('ls_on', 'float', 0.99999, False,
         """Start the backtracking line-search by reducing the step, if
            :math:`||f(x^i)|| / ||f(x^{i-1})||` is larger than `ls_on`."""),
        ('ls0', '0.0 < float < 1.0', 1.0, False,
         'The initial step.'),
        ('ls_red', '0.0 < float < 1.0', 0.5, False,
         'The step reduction factor in case of correct residual assembling.'),
        ('ls_red_warp', '0.0 < float < 1.0', 0.1, False,
         """The step reduction factor in case of failed residual assembling
            (e.g. the "warp violation" error caused by a negative volume
            element resulting from too large deformations)."""),
        ('ls_min', '0.0 < float < 1.0', 1e-5, False,
         'The minimum step reduction factor.'),
        ('check', '0, 1 or 2', 0, False,
         """If >= 1, check the tangent matrix using finite differences.  If 2,
            plot the resulting sparsity patterns."""),
        ('delta', 'float', 1e-6, False,
         r"""If `check >= 1`, the finite difference matrix is taken as
            :math:`A_{ij} = \frac{f_i(x_j + \delta) - f_i(x_j - \delta)}{2
            \delta}`."""),
        ('output', 'function', None, False,
         """If given, use it instead of :func:`output()
            <sfepy.base.base.output()>` function."""),
        ('yscales', 'list of str', ['linear', 'log', 'log', 'linear'], False,
         'The list of four convergence log subplot scales.'),
        ('log', 'dict or None', None, False,
         """If not None, log the convergence according to the configuration in
            the following form: ``{'text' : 'log.txt', 'plot' : 'log.pdf'}``.
            Each of the dict items can be None."""),
    ]

    def __init__(self, conf, **kwargs):
        OptimizationSolver.__init__(self, conf, **kwargs)

        conf = self.conf

        log = get_logging_conf(conf)
        conf.log = log = Struct(name='log_conf', **log)
        conf.is_any_log = (log.text is not None) or (log.plot is not None)

        if conf.is_any_log:
            self.log = Log([[r'$||\Psi||$'], [r'$||\nabla \Psi||$'],
                            [r'$\alpha$'], ['iteration']],
                           xlabels=['', '', 'all iterations', 'all iterations'],
                           yscales=conf.yscales,
                           is_plot=conf.log.plot is not None,
                           log_filename=conf.log.text,
                           formats=[['%.8e'], ['%.3e'], ['%.3e'], ['%d']])
        else:
            self.log = None

    def __call__(self, x0, conf=None, obj_fun=None, obj_fun_grad=None,
                 status=None, obj_args=None):
        conf = get_default(conf, self.conf)
        obj_fun = get_default(obj_fun, self.obj_fun)
        obj_fun_grad = get_default(obj_fun_grad, self.obj_fun_grad)
        status = get_default(status, self.status)
        obj_args = get_default(obj_args, self.obj_args)

        if conf.output:
            globals()['output'] = conf.output

        output('entering optimization loop...')

        nc_of, tt_of, fn_of = wrap_function(obj_fun, obj_args)
        nc_ofg, tt_ofg, fn_ofg = wrap_function(obj_fun_grad, obj_args)

        timer = Timer()
        time_stats = {'of' : tt_of, 'ofg': tt_ofg, 'check' : []}

        ofg = None

        it = 0
        xit = x0.copy()
        while 1:

            of = fn_of(xit)

            if it == 0:
                of0 = ofit0 = of_prev = of
                of_prev_prev = of + 5000.0

            if ofg is None:
                ofg = fn_ofg(xit)

            if conf.check:
                timer.start()
                check_gradient(xit, ofg, fn_of, conf.delta, conf.check)
                time_stats['check'].append(timer.stop())

            ofg_norm = nla.norm(ofg, conf.norm)

            ret = conv_test(conf, it, of, ofit0, ofg_norm)
            if ret >= 0:
                break
            ofit0 = of

            ##
            # Backtrack (on errors).
            alpha = conf.ls0
            can_ls = True
            while 1:
                xit2 = xit - alpha * ofg
                aux = fn_of(xit2)

                if self.log is not None:
                    self.log(of, ofg_norm, alpha, it)

                if aux is None:
                    alpha *= conf.ls_red_warp
                    can_ls = False
                    output('warp: reducing step (%f)' % alpha)
                elif conf.ls and conf.ls_method == 'backtracking':
                    if aux < of * conf.ls_on: break
                    alpha *= conf.ls_red
                    output('backtracking: reducing step (%f)' % alpha)
                else:
                    of_prev_prev = of_prev
                    of_prev = aux
                    break

                if alpha < conf.ls_min:
                    if aux is None:
                        raise RuntimeError('giving up...')
                    output('linesearch failed, continuing anyway')
                    break

            # These values are modified by the line search, even if it fails
            of_prev_bak = of_prev
            of_prev_prev_bak = of_prev_prev

            if conf.ls and can_ls and conf.ls_method == 'full':
                output('full linesearch...')
                alpha, fc, gc, of_prev, of_prev_prev, ofg1 = \
                    linesearch.line_search(fn_of,fn_ofg,xit,
                                           -ofg,ofg,of_prev,of_prev_prev,
                                           c2=0.4)
                if alpha is None:  # line search failed -- use different one.
                    alpha, fc, gc, of_prev, of_prev_prev, ofg1 = \
                        sopt.line_search(fn_of,fn_ofg,xit,
                                         -ofg,ofg,of_prev_bak,
                                         of_prev_prev_bak)
                    if alpha is None or alpha == 0:
                        # This line search also failed to find a better
                        # solution.
                        ret = 3
                        break
                output(' -> alpha: %.8e' % alpha)
            else:
                if conf.ls_method == 'full':
                    output('full linesearch off (%s and %s)'
                           % (conf.ls, can_ls))
                ofg1 = None

            if self.log is not None:
                self.log.plot_vlines(color='g', linewidth=0.5)

            xit = xit - alpha * ofg
            if ofg1 is None:
                ofg = None
            else:
                ofg = ofg1.copy()

            for key, val in time_stats.items():
                if len(val):
                    output('%10s: %7.2f [s]' % (key, val[-1]))

            it = it + 1

        output('status:               %d' % ret)
        output('initial value:        %.8e' % of0)
        output('current value:        %.8e' % of)
        output('iterations:           %d' % it)
        output('function evaluations: %d in %.2f [s]'
               % (nc_of[0], nm.sum(time_stats['of'])))
        output('gradient evaluations: %d in %.2f [s]'
               % (nc_ofg[0], nm.sum(time_stats['ofg'])))

        if self.log is not None:
            self.log(of, ofg_norm, alpha, it)

            if conf.log.plot is not None:
                self.log(save_figure=conf.log.plot, finished=True)

            else:
                self.log(finished=True)

        if status is not None:
            status['log'] = self.log
            status['status'] = status
            status['of0'] = of0
            status['of'] = of
            status['it'] = it
            status['nc_of'] = nc_of[0]
            status['nc_ofg'] = nc_ofg[0]
            status['time_stats'] = time_stats

        return xit

class ScipyFMinSolver(OptimizationSolver):
    """
    Interface to SciPy optimization solvers scipy.optimize.fmin_*.
    """
    name = 'nls.scipy_fmin_like'

    _i_max_name  = {
        'fmin' : 'maxiter',
        'fmin_bfgs' : 'maxiter',
        'fmin_cg' : 'maxiter',
        'fmin_cobyla' : 'maxfun',
        'fmin_l_bfgs_b' : 'maxfun',
        'fmin_ncg' : 'maxiter',
        'fmin_powell' : 'maxiter',
        'fmin_slsqp' : 'iter',
        'fmin_tnc' : 'maxfun',
    }
    _has_grad = ('fmin_bfgs', 'fmin_cg', 'fmin_l_bfgs_b', 'fmin_ncg',
                 'fmin_slsqp', 'fmin_tnc')

    _parameters = [
        ('method',
         '{%s}' % ', '.join(sorted(repr(ii) for ii in _i_max_name.keys())),
         'fmin', False,
         'The actual optimization method to use.'),
        ('i_max', 'int', 10, False,
         'The maximum number of iterations.'),
        ('*', '*', None, False,
         'Additional parameters supported by the method.'),
    ]

    def __init__(self, conf, **kwargs):
        OptimizationSolver.__init__(self, conf, **kwargs)
        self.set_method(self.conf)

    def set_method(self, conf):
        import scipy.optimize as so

        try:
            solver = getattr(so, conf.method)
        except AttributeError:
            raise ValueError('scipy solver %s does not exist!' % conf.method)

        self.solver = solver

    def __call__(self, x0, conf=None, obj_fun=None, obj_fun_grad=None,
                 status=None, obj_args=None):
        import inspect

        if conf is not None:
            self.set_method(conf)

        else:
            conf = self.conf

        obj_fun = get_default(obj_fun, self.obj_fun)
        obj_fun_grad = get_default(obj_fun_grad, self.obj_fun_grad)
        status = get_default(status, self.status)
        obj_args = get_default(obj_args, self.obj_args)

        timer = Timer(start=True)

        kwargs = {self._i_max_name[conf.method] : conf.i_max,
                  'args' : obj_args}

        if conf.method in self._has_grad:
            kwargs['fprime'] = obj_fun_grad

        if 'disp' in inspect.getargspec(self.solver)[0]:
            kwargs['disp'] = conf.verbose

        kwargs.update(self.build_solver_kwargs(conf))

        out = self.solver(obj_fun, x0, **kwargs)

        if status is not None:
            status['time_stats'] = timer.stop()

        return out
