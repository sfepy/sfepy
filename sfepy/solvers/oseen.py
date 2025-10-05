
import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import output, get_default, Struct
from sfepy.base.log import Log, get_logging_conf
from sfepy.base.timing import Timer
from sfepy.solvers.solvers import NonlinearSolver
from .nls import conv_test

class StabilizationFunction(Struct):
    """
    Definition of stabilization material function for the Oseen solver.

    Notes
    -----
    - tau_red <= 1.0; if tau is None: tau = tau_red * delta
    - diameter mode: 'edge': longest edge 'volume': volume-based, 'max': max. of
      previous
    """

    def __init__(self, name_map, gamma=None, delta=None, tau=None, tau_red=1.0,
                 tau_mul=1.0, delta_mul=1.0, gamma_mul=1.0,
                 diameter_mode='max'):
        Struct.__init__(self, name_map=name_map,
                        gamma=gamma, delta=delta, tau=tau,
                        tau_red=tau_red, tau_mul=tau_mul, delta_mul=delta_mul,
                        gamma_mul=gamma_mul, diameter_mode=diameter_mode)

    def setup(self, problem):
        """
        Setup common problem-dependent data.
        """
        variables = problem.get_variables()

        ns = self.name_map

        # Indices to the state vector.
        ii = {}
        ii['u'] = variables.get_indx(ns['u'])
        ii['us'] = variables.get_indx(ns['u'], reduced=True)
        ii['ps'] = variables.get_indx(ns['p'], reduced=True)
        self.indices = ii

        materials = problem.get_materials()

        # The viscosity.
        fluid_mat = materials[ns['fluid']]
        self.viscosity = fluid_mat.function()[ns['viscosity']]

        # The Friedrich's constant.
        self.c_friedrichs = problem.domain.get_diameter()
        self.sigma = 1e-12 # 1 / dt.

        self.b_norm = 1.0

    def get_maps(self):
        """
        Get the maps of names and indices of variables in state vector.
        """
        return self.name_map, self.indices

    def __call__(self, ts, coor, mode=None, term=None, problem=None,
                 b_norm=None, **kwargs):
        """
        The actual material function.
        """
        if mode != 'qp': return

        if not hasattr(self, 'viscosity'):
            self.setup(problem)

        ns = self.name_map

        # Update stored b_norm.
        self.b_norm = get_default(b_norm, self.b_norm)

        output('|b|_max (mat_fun):', self.b_norm)
        gamma = self.viscosity + self.b_norm * self.c_friedrichs

        data = {}
        if self.gamma is None:
            _gamma = self.gamma_mul * gamma

        else:
            _gamma = nm.asarray(self.gamma_mul * self.gamma, dtype=nm.float64)
        _gamma = nm.tile(_gamma, (coor.shape[0], 1, 1))

        if self.delta is None:
            # Element diameter modes.
            dm = {'edge': 0, 'volume': 1, 'max': 2}[self.diameter_mode]

            field = problem.fields[ns['velocity']]
            region = term.region
            vg, _ = field.get_mapping(region, term.integral, 'volume')
            cells = region.get_cells()
            d2 = problem.domain.get_element_diameters(cells, vg.volume, dm)

            self.diameters2 = d2

            val1 = min(1.0, 1.0 / self.sigma)
            val2 = self.sigma * self.c_friedrichs**2
            val3 = ((self.b_norm**2)
                    * min((self.c_friedrichs**2) / self.viscosity,
                          1.0 / self.sigma))

            n_qp = coor.shape[0] / self.diameters2.shape[0]
            diameters2 = nm.repeat(self.diameters2, n_qp)
            diameters2.shape = diameters2.shape + (1, 1)

            _delta = self.delta_mul * val1 * diameters2 / (_gamma + val2 + val3)

        else:
            val = nm.asarray(self.delta_mul * self.delta, dtype=nm.float64)
            _delta = nm.tile(val, (coor.shape[0], 1, 1))

        if self.tau is None:
            _tau = self.tau_red * _delta

        else:
            _tau = nm.asarray(self.tau_mul * self.tau, dtype=nm.float64)
            _tau = nm.tile(_tau, (coor.shape[0], 1, 1))

        data[ns['gamma']] = _gamma
        data[ns['delta']] = _delta
        data[ns['tau']] = _tau

        return data

def are_close(a, b, rtol=0.2, atol=1e-8):
    return False
#    return abs(a - b) <= max(atol, rtol * abs(b))

def scale_matrix(mtx, indx, factor):
    ptr0 = mtx.indptr[indx.start]
    ptr1 = mtx.indptr[indx.stop]
    mtx.data[ptr0:ptr1] *= factor

class Oseen(NonlinearSolver):
    """
    The Oseen solver for Navier-Stokes equations.
    """
    name = 'nls.oseen'

    _parameters = [
        ('stabil_mat', 'str', None, True,
         'The name of stabilization material.'),
        ('adimensionalize', 'bool', False, False,
         'If True, adimensionalize the problem (not implemented!).'),
        ('check_navier_stokes_residual', 'bool', False, False,
         'If True, check the Navier-Stokes residual after the nonlinear loop.'),
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
        ('lin_precision', 'float or None', None, False,
         """If not None, the linear system solution tolerances are set in each
            nonlinear iteration relative to the current residual norm by the
            `lin_precision` factor. Ignored for direct linear solvers."""),
    ]

    def __init__(self, conf, context=None, **kwargs):
        NonlinearSolver.__init__(self, conf, context=context, **kwargs)

        conf = self.conf

        log = get_logging_conf(conf)
        conf.log = log = Struct(name='log_conf', **log)
        conf.is_any_log = (log.text is not None) or (log.plot is not None)

        conf.problem = context

        conf = self.conf
        if conf.is_any_log:
            self.log = Log([[r'$||r||$'], ['iteration'],
                            [r'$\gamma$', r'$\max(\delta)$', r'$\max(\tau)$']],
                           xlabels=['', '', 'all iterations'],
                           ylabels=[r'$||r||$', 'iteration', 'stabilization'],
                           yscales=['log', 'linear', 'log'],
                           is_plot=conf.log.plot is not None,
                           log_filename=conf.log.text,
                           formats=[['%.8e'], ['%d'],
                                    ['%.8e', '%.8e', '%.8e']])

        else:
            self.log = None

    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, status=None, problem=None):
        """
        Oseen solver is problem-specific - it requires a Problem instance.
        """
        conf = get_default(conf, self.conf)
        fun = get_default(fun, self.fun)
        fun_grad = get_default(fun_grad, self.fun_grad)
        lin_solver = get_default(lin_solver, self.lin_solver)
        status = get_default(status, self.status)
        problem = get_default(problem, conf.problem,
                              '`problem` parameter needs to be set!')

        timer = Timer()
        time_stats = {}

        stabil = problem.get_materials()[conf.stabil_mat]
        ns, ii = stabil.function.function.get_maps()

        variables = problem.get_variables()
        make_full_vec = variables.make_full_vec

        output('problem size:')
        output('    velocity: %s' % ii['us'])
        output('    pressure: %s' % ii['ps'])

        vec_x = vec_x0.copy()
        vec_x_prev = vec_x0.copy()
        vec_dx = None

        if self.log is not None:
            self.log.plot_vlines(color='r', linewidth=1.0)

        err0 = -1.0
        it = 0
        ls_n_iter = 0
        ls_status = {}
        while 1:
            vec_x_prev_f = make_full_vec(vec_x_prev)
            variables[ns['b']].set_data(vec_x_prev_f, variables[ns['u']].indx)

            vec_b = vec_x_prev_f[ii['u']]
            b_norm = nla.norm(vec_b, nm.inf)
            output('|b|_max: %.12e' % b_norm)

            vec_x_f = make_full_vec(vec_x)
            vec_u = vec_x_f[ii['u']]
            u_norm = nla.norm(vec_u, nm.inf)
            output('|u|_max: %.2e' % u_norm)

            stabil.function.set_extra_args(b_norm=b_norm)
            stabil.time_update(None, problem.equations, mode='force',
                               problem=problem)
            max_pars = stabil.reduce_on_datas(lambda a, b: max(a, b.max()))
            output('stabilization parameters:')
            output('                   gamma: %.12e' % max_pars[ns['gamma']])
            output('            max(delta): %.12e' % max_pars[ns['delta']])
            output('              max(tau): %.12e' % max_pars[ns['tau']])

            if (not are_close(b_norm, 1.0)) and conf.adimensionalize:
                adimensionalize = True
            else:
                adimensionalize = False

            timer.start()
            try:
                vec_r = fun(vec_x)
            except ValueError:
                ok = False
            else:
                ok = True
            time_stats['residual'] = timer.stop()
            if ok:
                err = nla.norm(vec_r)
                if it == 0:
                    err0 = err;
                else:
                    err += nla.norm(vec_dx)
            else: # Failure.
                output('residual computation failed for iter %d!' % it)
                raise RuntimeError('giving up...')

            if self.log is not None:
                self.log(err, it,
                         max_pars[ns['gamma']], max_pars[ns['delta']],
                         max_pars[ns['tau']])

            condition = conv_test(conf, it, err, err0)
            if condition >= 0:
                break

            if adimensionalize:
                output('adimensionalizing')
                ## mat.viscosity = viscosity / b_norm
                ## vec_r[indx_us] /= b_norm

            timer.start()
            try:
                mtx_a = fun_grad(vec_x)
            except ValueError:
                ok = False
            else:
                ok = True
            time_stats['matrix'] = timer.stop()
            if not ok:
                raise RuntimeError('giving up...')

            timer.start()
            vec_dx = lin_solver(vec_r, x0=vec_x, mtx=mtx_a, status=ls_status)
            time_stats['solve'] = timer.stop()
            ls_n_iter += ls_status['n_iter']

            vec_e = mtx_a * vec_dx - vec_r
            lerr = nla.norm(vec_e)
            if lerr > (conf.eps_a * conf.lin_red):
                output('linear system not solved! (err = %e)' % lerr)

            if adimensionalize:
                output('restoring pressure...')
                ## vec_dx[indx_ps] *= b_norm

            dx_norm = nla.norm(vec_dx)
            output('||dx||: %.2e' % dx_norm)

            for kv in time_stats.items():
                output('%10s: %7.2f [s]' % kv)

            vec_x_prev = vec_x.copy()
            vec_x -= vec_dx
            it += 1

        if conf.check_navier_stokes_residual:

            t1 = '+ dw_div_grad.%s.%s(%s.viscosity, %s, %s)' \
                 % (ns['i2'], ns['omega'], ns['fluid'], ns['v'], ns['u'])
            # t2 = '+ dw_lin_convect.%s(%s, %s, %s)' % (ns['omega'],
            #                                           ns['v'], b_name, ns['u'])
            t2 = '+ dw_convect.%s.%s(%s, %s)' % (ns['i2'], ns['omega'],
                                                 ns['v'], ns['u'])
            t3 = '- dw_stokes.%s.%s(%s, %s)' % (ns['i1'], ns['omega'],
                                                ns['v'], ns['p'])
            t4 = 'dw_stokes.%s.%s(%s, %s)' % (ns['i1'], ns['omega'],
                                              ns['u'], ns['q'])
            equations = {
                'balance' : ' '.join((t1, t2, t3)),
                'incompressibility' : t4,
            }
            problem.set_equations(equations)
            try:
                vec_rns0 = fun(vec_x0)
                vec_rns = fun(vec_x)
            except ValueError:
                ok = False
            else:
                ok = True
            if not ok:
                output('Navier-Stokes residual computation failed!')
                err_ns = err_ns0 = None
            else:
                err_ns0 = nla.norm(vec_rns0)
                err_ns = nla.norm(vec_rns)
            output('Navier-Stokes residual0: %.8e' % err_ns0)
            output('Navier-Stokes residual : %.8e' % err_ns)
            output('b - u: %.8e' % nla.norm(vec_b - vec_u))
            output(condition)

        else:
            err_ns = None

        if status is not None:
            status['time_stats'] = time_stats
            status['err0'] = err0
            status['err'] = err
            status['err_ns'] = err_ns
            status['condition'] = condition
            status['n_iter'] = it
            status['ls_n_iter'] = -1
            status['time'] = sum(v for v in time_stats.values())

        if conf.log.plot is not None:
            if self.log is not None:
                self.log(save_figure=conf.log.plot)

        return vec_x
