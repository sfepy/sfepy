import time

import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import output, get_default, pause, debug, Struct
from sfepy.base.log import Log, get_logging_conf
from sfepy.solvers.solvers import make_get_conf, NonlinearSolver

def check_tangent_matrix( conf, vec_x0, fun, fun_grad ):
    """Verify the correctness of the tangent matrix as computed by fun_grad()
    by comparing it with its finite difference approximation evaluated by
    repeatedly calling fun() with vec_x items perturbed by a small delta."""
    vec_x = vec_x0.copy()
    delta = conf.delta

    vec_r = fun( vec_x ) # Update state.
    mtx_a0 = fun_grad( vec_x )

    mtx_a = mtx_a0.tocsc()
    mtx_d = mtx_a.copy()
    mtx_d.data[:] = 0.0

    vec_dx = nm.zeros_like( vec_r )

    for ic in range( vec_dx.shape[0] ):
        vec_dx[ic] = delta
        xx = vec_x.copy() - vec_dx
        vec_r1 = fun( xx )

        vec_dx[ic] = -delta
        xx = vec_x.copy() - vec_dx
        vec_r2 = fun( xx )

        vec_dx[ic] = 0.0;

        vec = 0.5 * (vec_r2 - vec_r1) / delta

##         ir = mtx_a.indices[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]]
##         for ii in ir:
##             mtx_d[ii,ic] = vec[ii]

        ir = mtx_a.indices[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]]
        mtx_d.data[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]] = vec[ir]


    vec_r = fun( vec_x ) # Restore.

    tt = time.clock()
    print mtx_a, '.. analytical'
    print mtx_d, '.. difference'
    import sfepy.base.plotutils as plu
    plu.plot_matrix_diff( mtx_d, mtx_a, delta, ['difference', 'analytical'],
                        conf.check )

    return time.clock() - tt

##
# c: 02.12.2005, r: 02.04.2008
def conv_test( conf, it, err, err0 ):

    status = -1
    if (abs( err0 ) < conf.macheps):
        err_r = 0.0
    else:
        err_r = err / err0

    output( 'nls: iter: %d, residual: %e (rel: %e)' % (it, err, err_r) )
    if it > 0:
        if (err < conf.eps_a) and (err_r < conf.eps_r):
            status = 0
    else:
        if err < conf.eps_a:
            status = 0

    if (status == -1) and (it >= conf.i_max):
        status = 1

    return status

##
# 10.10.2007, c
class Newton( NonlinearSolver ):
    name = 'nls.newton'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are set to default values for a linear problem.

        Example configuration, all items::

            solver_1 = {
                'name' : 'newton',
                'kind' : 'nls.newton',

                'i_max' : 2,
                'eps_a' : 1e-8,
                'eps_r' : 1e-2,
                'macheps' : 1e-16,
                'lin_red' : 1e-2, # Linear system error < (eps_a * lin_red).
                'lin_precision' : None,
                'ls_red' : 0.1,
                'ls_red_warp' : 0.001,
                'ls_on' : 0.99999,
                'ls_min' : 1e-5,
                'give_up_warp' : False,
                'check' : 0,
                'delta' : 1e-6,
                'is_plot' : False,
                'log' : None, # 'nonlinear' or 'linear' (ignore i_max)
                'problem' : 'nonlinear',
            }
        """
        get = make_get_conf(conf, kwargs)
        common = NonlinearSolver.process_conf(conf)

        log = get_logging_conf(conf)
        log = Struct(name='log_conf', **log)
        is_any_log = (log.text is not None) or (log.plot is not None)

        return Struct(i_max=get('i_max', 1),
                      eps_a=get('eps_a', 1e-10),
                      eps_r=get('eps_r', 1.0),
                      macheps=get('macheps', nm.finfo(nm.float64).eps),
                      lin_red=get('lin_red', 1.0),
                      lin_precision=get('lin_precision', None),
                      ls_red=get('ls_red', 0.1),
                      ls_red_warp=get('ls_red_warp', 0.001),
                      ls_on=get('ls_on', 0.99999),
                      ls_min=get('ls_min', 1e-5),
                      give_up_warp=get('give_up_warp', False),
                      check=get('check', 0),
                      delta=get('delta', 1e-6),
                      is_plot=get('is_plot', False),
                      problem=get('problem', 'nonlinear'),
                      log=log,
                      is_any_log=is_any_log) + common

    def __init__(self, conf, **kwargs):
        NonlinearSolver.__init__( self, conf, **kwargs )

        conf = self.conf
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

        Solves :math:`f(x) = 0` by the Newton method with backtracking
        linesearch.

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
        * Setting `conf.problem == 'linear'` means 1 iteration and no
          rezidual check!
        """
        import sfepy.base.plotutils as plu
        conf = get_default( conf, self.conf )
        fun = get_default( fun, self.fun )
        fun_grad = get_default( fun_grad, self.fun_grad )
        lin_solver = get_default( lin_solver, self.lin_solver )
        iter_hook = get_default(iter_hook, self.iter_hook)
        status = get_default( status, self.status )

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
                    vec_r = fun( vec_x )

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
                        err = nla.norm( vec_r )
                    except:
                        output( 'infs or nans in the residual:', vec_r )
                        output( nm.isfinite( vec_r ).all() )
                        debug()

                    if self.log is not None:
                        self.log(err, it)

                    if it == 0:
                        err0 = err;
                        break
                    if err < (err_last * conf.ls_on): break
                    red = conf.ls_red;
                    output( 'linesearch: iter %d, (%.5e < %.5e) (new ls: %e)'\
                            % (it, err, err_last * conf.ls_on, red * ls) )
                else: # Failure.
                    if conf.give_up_warp:
                        output('giving up!')
                        break

                    red = conf.ls_red_warp;
                    output(  'rezidual computation failed for iter %d'
                             ' (new ls: %e)!' % (it, red * ls) )

                if ls < conf.ls_min:
                    output( 'linesearch failed, continuing anyway' )
                    break

                ls *= red;

                vec_dx = ls * vec_dx0;
                vec_x = vec_x_last.copy() - vec_dx
            # End residual loop.

            if self.log is not None:
                self.log.plot_vlines([1], color='g', linewidth=0.5)

            err_last = err;
            vec_x_last = vec_x.copy()

            condition = conv_test( conf, it, err, err0 )
            if condition >= 0:
                break

            if (not ok) and conf.give_up_warp:
                condition = 2
                break

            tt = time.clock()
            if conf.problem == 'nonlinear':
                mtx_a = fun_grad(vec_x)

            else:
                mtx_a = fun_grad( 'linear' )

            time_stats['matrix'] = time.clock() - tt

            if conf.check:
                tt = time.clock()
                wt = check_tangent_matrix( conf, vec_x, fun, fun_grad )
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
                output( '%10s: %7.2f [s]' % kv )

            vec_e = mtx_a * vec_dx - vec_r
            lerr = nla.norm( vec_e )
            if lerr > lin_red:
                output('linear system not solved! (err = %e < %e)'
                       % (lerr, lin_red))

            vec_x -= vec_dx

            if conf.is_plot:
                plu.plt.ion()
                plu.plt.gcf().clear()
                plu.plt.subplot( 2, 2, 1 )
                plu.plt.plot( vec_x_last )
                plu.plt.ylabel( r'$x_{i-1}$' )
                plu.plt.subplot( 2, 2, 2 )
                plu.plt.plot( vec_r )
                plu.plt.ylabel( r'$r$' )
                plu.plt.subplot( 2, 2, 4 )
                plu.plt.plot( vec_dx )
                plu.plt.ylabel( r'$\_delta x$' )
                plu.plt.subplot( 2, 2, 3 )
                plu.plt.plot( vec_x )
                plu.plt.ylabel( r'$x_i$' )
                plu.plt.draw()
                plu.plt.ioff()
                pause()

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

class ScipyBroyden( NonlinearSolver ):
    """Interface to Broyden and Anderson solvers from scipy.optimize."""

    name = 'nls.scipy_broyden_like'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Missing items are left to scipy defaults. Unused options are ignored.

        Example configuration, all items::

            solver_1 = {
                'name' : 'broyden',
                'kind' : 'nls.scipy_broyden_like',

                'method'  : 'broyden3',
                'i_max'   : 10,
                'alpha'   : 0.9,
                'M'       : 5,
                'w0'      : 0.1,
                'verbose' : True,
            }
        """
        get = make_get_conf(conf, kwargs)
        common = NonlinearSolver.process_conf(conf)

        return Struct(method=get('method', 'broyden3'),
                      i_max=get('i_max', 10),
                      alpha=get('alpha', 0.9),
                      M=get('M', 5),
                      w0=get('w0', 0.1),
                      verbose=get('verbose', False)) + common

    def __init__( self, conf, **kwargs ):
        NonlinearSolver.__init__( self, conf, **kwargs )
        self.set_method( self.conf )

    def set_method( self, conf ):
        import scipy.optimize as so

        try:
            solver = getattr( so, conf.method )
        except AttributeError:
            output( 'scipy solver %s does not exist!' % conf.method )
            output( 'using broyden3 instead' )
            solver = so.broyden3
        self.solver = solver

    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None):
        if conf is not None:
            self.set_method( conf )
        else:
            conf = self.conf
        fun = get_default( fun, self.fun )
        status = get_default( status, self.status )

        tt = time.clock()

        kwargs = {'iter' : conf.i_max,
                  'alpha' : conf.alpha,
                  'verbose' : conf.verbose}

        if conf.method == 'broyden_generalized':
            kwargs.update( {'M' : conf.M} )

        elif conf.method in ['anderson', 'anderson2']:
            kwargs.update( {'M' : conf.M, 'w0' : conf.w0} )

        vec_x = self.solver( fun, vec_x0, **kwargs )
        
        if status is not None:
            status['time_stats'] = time.clock() - tt

        return vec_x
