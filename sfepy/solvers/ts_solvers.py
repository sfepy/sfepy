"""
Time stepping solvers.
"""
from inspect import signature
from functools import partial
import os.path as osp
import numpy as nm

from sfepy.base.base import (get_default, output, assert_, Struct)
from sfepy.base.timing import Timer
from sfepy.linalg.utils import output_array_stats
from sfepy.solvers.solvers import TimeSteppingSolver, NonlinearSolver
from sfepy.solvers.ls import RMMSolver
from sfepy.solvers.ts_controllers import FixedTSC
from sfepy.solvers.ts import TimeStepper, VariableTimeStepper

def standard_ts_call(call):
    """
    Decorator handling argument preparation and timing for time-stepping
    solvers.
    """

    def _standard_ts_call(self, vec0=None, nls=None,
                         init_fun=None, prestep_fun=None, poststep_fun=None,
                         status=None, log_nls_status=False, **kwargs):
        timer = Timer(start=True)

        nls = get_default(nls, self.nls,
                          'nonlinear solver has to be specified!')

        init_fun = get_default(init_fun, lambda ts, vec0: vec0)
        prestep_fun = get_default(prestep_fun, lambda ts, vec: None)
        poststep_fun = get_default(poststep_fun, lambda ts, vec: None)

        _nls = nls
        if status is not None:
            nls_status = status.get('nls_status')
            if nls_status is not None:
                if log_nls_status:
                    pb = self.context
                    filename_log = osp.join(pb.output_dir,
                                            pb.ofn_trunk + '_log.csv')
                    status.log_file = open(filename_log, 'wt',
                                           encoding="utf-8")

                class _TimingNLS(type(nls)):
                    def __call__(self, *args, **kwargs):
                        # Call the original nls...
                        out = super().__call__(*args, **kwargs)

                        log_stats = {}
                        # ...and collect its time stats.
                        time_stats = nls_status.get('time_stats')
                        if time_stats is not None:
                            all_stats = status.setdefault('time_stats', {})
                            for key, val in time_stats.items():
                                all_stats.setdefault(key, 0.0)
                                all_stats[key] += val

                            log_stats.update(time_stats)

                            all_stats.setdefault('time', 0.0)
                            all_stats['time'] += nls_status.get('time', 0.0)

                        # ...and collect its step stats.
                        all_stats = status.setdefault('step_stats', [])
                        _nls_status = nls_status.copy()
                        _nls_status.step = self.context.ts.step
                        _nls_status.step_time = self.context.ts.time
                        all_stats.append(_nls_status)
                        log_stats.update(_nls_status.to_dict())

                        if getattr(status, 'log_file', None) is not None:
                            if 'time_stats' in log_stats:
                                del log_stats['time_stats']

                            keys = list(log_stats.keys())
                            keys.sort()

                            if status.log_file.tell() == 0:
                                status.log_file.write(','.join(keys) + '\n')

                            log_vals = [f'{log_stats[key]}' for key in keys]
                            status.log_file.write(','.join(log_vals) + '\n')
                            status.log_file.flush()

                        return out

                _nls = nls.copy()
                _nls.__class__ = _TimingNLS

        result = call(self, vec0=vec0, nls=_nls, init_fun=init_fun,
                      prestep_fun=prestep_fun, poststep_fun=poststep_fun,
                      status=status, **kwargs)

        elapsed = timer.stop()
        if status is not None:
            status['time'] = elapsed
            status['n_step'] = self.ts.n_step
            if getattr(status, 'log_file', None) is not None:
                status.log_file.close()

        return result

    return _standard_ts_call

#
# General solvers.
#
class StationarySolver(TimeSteppingSolver):
    """
    Solver for stationary problems without time stepping.

    This class is provided to have a unified interface of the time stepping
    solvers also for stationary problems.
    """
    name = 'ts.stationary'

    def __init__(self, conf, nls=None, context=None, **kwargs):
        TimeSteppingSolver.__init__(self, conf, nls=nls, context=context,
                                    **kwargs)

        self.ts = TimeStepper(0.0, 1.0, n_step=1, is_quasistatic=True)

    @standard_ts_call
    def __call__(self, vec0=None, nls=None, init_fun=None, prestep_fun=None,
                 poststep_fun=None, status=None, **kwargs):
        ts = self.ts
        nls = get_default(nls, self.nls)

        vec0 = init_fun(ts, vec0)

        vec0 = prestep_fun(ts, vec0)

        vec = nls(vec0)

        vec = poststep_fun(ts, vec)

        return vec

class SimpleTimeSteppingSolver(TimeSteppingSolver):
    """
    Implicit time stepping solver with a fixed time step.
    """
    name = 'ts.simple'

    _parameters = [
        ('t0', 'float', 0.0, False,
         'The initial time.'),
        ('t1', 'float', 1.0, False,
         'The final time.'),
        ('dt', 'float', None, False,
         'The time step. Used if `n_step` is not given.'),
        ('n_step', 'int', 10, False,
         'The number of time steps. Has precedence over `dt`.'),
        ('quasistatic', 'bool', False, False,
         """If True, assume a quasistatic time-stepping. Then the non-linear
            solver is invoked also for the initial time."""),
    ]

    def __init__(self, conf, nls=None, context=None, **kwargs):
        TimeSteppingSolver.__init__(self, conf, nls=nls, context=context,
                                    **kwargs)
        self.ts = TimeStepper.from_conf(self.conf)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format
        self.verbose = self.conf.verbose

    def solve_step0(self, nls, vec0):
        if self.conf.quasistatic:
            vec = nls(vec0)

        else:
            res = nls.fun(vec0)
            err = nm.linalg.norm(res)
            output('initial residual: %e' % err, verbose=self.verbose)
            vec = vec0.copy()

        return vec

    def solve_step(self, ts, nls, vec, prestep_fun=None):
        return nls(vec)

    def output_step_info(self, ts):
        output(self.format % (ts.time, ts.step + 1, ts.n_step),
               verbose=self.verbose)

    @standard_ts_call
    def __call__(self, vec0=None, nls=None, init_fun=None, prestep_fun=None,
                 poststep_fun=None, status=None, **kwargs):
        """
        Solve the time-dependent problem.
        """
        ts = self.ts
        nls = get_default(nls, self.nls)

        vec0 = init_fun(ts, vec0)

        self.output_step_info(ts)
        if ts.step == 0:
            vec0 = prestep_fun(ts, vec0)

            vec = self.solve_step0(nls, vec0)

            vec = poststep_fun(ts, vec)
            ts.advance()

        else:
            vec = vec0

        for step, time in ts.iter_from(ts.step):
            self.output_step_info(ts)

            vec = prestep_fun(ts, vec)

            vect = self.solve_step(ts, nls, vec, prestep_fun)

            vect = poststep_fun(ts, vect)

            vec = vect

        return vec

def get_min_dt(adt):
    red = adt.red
    while red >= adt.red_max:
        red *= adt.red_factor

    dt = adt.dt0 * red

    return dt

def adapt_time_step(ts, status, adt, context=None, verbose=False):
    """
    Adapt the time step of `ts` according to the exit status of the
    nonlinear solver.

    The time step dt is reduced, if the nonlinear solver did not converge. If it
    converged in less then a specified number of iterations for several time
    steps, the time step is increased. This is governed by the following
    parameters:

    - `red_factor` : time step reduction factor
    - `red_max` : maximum time step reduction factor
    - `inc_factor` : time step increase factor
    - `inc_on_iter` : increase time step if the nonlinear solver converged in
      less than this amount of iterations...
    - `inc_wait` : ...for this number of consecutive time steps

    Parameters
    ----------
    ts : VariableTimeStepper instance
        The time stepper.
    status : IndexedStruct instance
        The nonlinear solver exit status.
    adt : Struct instance
        The object with the adaptivity parameters of the time-stepping solver
        such as `red_factor` (see above) as attributes.
    context : object, optional
        The context can be used in user-defined adaptivity functions. Not used
        here.

    Returns
    -------
    is_break : bool
        If True, the adaptivity loop should stop.
    """
    is_break = False

    if status.condition == 0:
        if status.n_iter <= adt.inc_on_iter:
            adt.wait += 1

            if adt.wait > adt.inc_wait:
                if adt.red < 1.0:
                    adt.red = adt.red * adt.inc_factor
                    ts.set_time_step(adt.dt0 * adt.red)
                    output('+++++ new time step: %e +++++' % ts.dt,
                           verbose=verbose)
                adt.wait = 0

        else:
            adt.wait = 0

        is_break = True

    else:
        adt.red = adt.red * adt.red_factor
        if adt.red < adt.red_max:
            is_break = True

        else:
            ts.set_time_step(adt.dt0 * adt.red, update_time=True)
            output('----- new time step: %e -----' % ts.dt,
                   verbose=verbose)
            adt.wait = 0

    return is_break

class AdaptiveTimeSteppingSolver(SimpleTimeSteppingSolver):
    """
    Implicit time stepping solver with an adaptive time step.

    Either the built-in or user supplied function can be used to adapt the time
    step.
    """
    name = 'ts.adaptive'

    _parameters = SimpleTimeSteppingSolver._parameters + [
        ('adapt_fun', 'callable(ts, status, adt, context, verbose)',
         None, False,
         """If given, use this function to set the time step in `ts`. The
            function return value is a bool - if True, the adaptivity loop
            should stop. The other parameters below are collected in `adt`,
            `status` is the nonlinear solver status, `context` is
            a user-defined context and `verbose` is a verbosity flag.
            Solvers created by
            :class:`Problem <sfepy.discrete.problem.Problem>` use the
            Problem instance as the context."""),
        ('dt_red_factor', 'float', 0.2, False,
         'The time step reduction factor.'),
        ('dt_red_max', 'float', 1e-3, False,
         'The maximum time step reduction factor.'),
        ('dt_inc_factor', 'float', 1.25, False,
         'The time step increase factor.'),
        ('dt_inc_on_iter', 'int', 4, False,
         """Increase the time step if the nonlinear solver converged in less
            than this amount of iterations for `dt_inc_wait` consecutive time
            steps."""),
        ('dt_inc_wait', 'int', 5, False,
         'The number of consecutive time steps, see `dt_inc_on_iter`.'),
    ]

    def __init__(self, conf, nls=None, context=None, **kwargs):
        TimeSteppingSolver.__init__(self, conf, nls=nls, context=context,
                                    **kwargs)

        self.ts = VariableTimeStepper.from_conf(self.conf)

        get = self.conf.get
        adt = Struct(red_factor=get('dt_red_factor', 0.2),
                     red_max=get('dt_red_max', 1e-3),
                     inc_factor=get('dt_inc_factor', 1.25),
                     inc_on_iter=get('dt_inc_on_iter', 4),
                     inc_wait=get('dt_inc_wait', 5),
                     red=1.0, wait=0, dt0=0.0)
        self.adt = adt

        adt.dt0 = self.ts.get_default_time_step()
        self.ts.set_n_digit_from_min_dt(get_min_dt(adt))

        self.format = '====== time %e (dt %e, wait %d, step %d of %d) ====='
        self.verbose = self.conf.verbose

        self.adapt_time_step = self.conf.adapt_fun
        if self.adapt_time_step is None:
            self.adapt_time_step = adapt_time_step

    def solve_step(self, ts, nls, vec, prestep_fun):
        """
        Solve a single time step.
        """
        while 1:
            vect = nls(vec, status=nls.status)

            is_break = self.adapt_time_step(ts, nls.status, self.adt,
                                            self.context, verbose=self.verbose)

            if is_break:
                break

            vec = prestep_fun(ts, vec)

        return vect

    def output_step_info(self, ts):
        output(self.format % (ts.time, ts.dt, self.adt.wait,
                              ts.step + 1, ts.n_step),
               verbose=self.verbose)

#
# Elastodynamics solvers.
#
def transform_equations_ed(equations, materials):
    """
    Transform equations and variables for :class:`ElastodynamicsBaseTS`-based
    time stepping solvers. The displacement variable name is automatically
    detected by seeking the second time derivative, i.e. the 'dd' prefix in
    variable names.
    """
    from sfepy.terms.terms import Terms, Term
    from sfepy.discrete.variables import FieldVariable
    from sfepy.discrete.equations import Equations, Equation

    eq_mass = []
    eq_damping = []
    eq_other = []
    new_var_primary = set()
    variables = equations.variables
    for eq in equations:
        for term in eq.terms:
            for name in term.names.state:
                der = term.arg_derivatives[name]
                if ((der is not None) and
                    isinstance(der, int) and
                    (der > 0)):
                    new_var_primary.add(name)
                    if der == 2:
                        eq_mass.append(term)

                    else:
                        eq_damping.append(term)

                else:
                    eq_other.append(term)

                continue

    assert len(new_var_primary) == 1
    uname = new_var_primary.pop()
    vname = variables[uname].dual_var_name
    duname = 'd' + uname
    dduname = 'dd' + uname
    dvname = 'd' + vname
    ddvname = 'dd' + vname
    new_variables = variables.copy()
    if new_variables[uname]._order is None:
        raise ValueError('state variable orders have to be specified when using'
                         ' auto_transform_equations!')
    num = len(new_variables.state)
    new_variables.extend([
        FieldVariable(duname, 'unknown', variables[uname].field, order=num),
        FieldVariable(dduname, 'unknown', variables[uname].field, order=num + 1),
        FieldVariable(dvname, 'test', variables[uname].field,
                      primary_var_name=duname),
        FieldVariable(ddvname, 'test', variables[uname].field,
                      primary_var_name=dduname),
    ])
    for it, term in enumerate(eq_mass.copy()):
        aux = ','.join([ii.strip() for ii in term.arg_str.split(',')])
        arg_str = aux.replace(f',{vname},', f',{ddvname},')
        new_term = term.__class__(term.name, arg_str, term.integral, term.region)
        new_term.setup(allow_derivatives=False)
        new_term.assign_args(new_variables, materials, user=None)


        eq_mass[it] = new_term

    for it, term in enumerate(eq_damping.copy()):
        aux = ','.join([ii.strip() for ii in term.arg_str.split(',')])
        arg_str = aux.replace(f',{vname},', f',{dvname},')
        new_term = term.__class__(term.name, arg_str, term.integral, term.region)
        new_term.setup(allow_derivatives=False)
        new_term.assign_args(new_variables, materials, user=None)


        eq_damping[it] = new_term

    if not len(eq_damping):
        mterm = eq_mass[0]
        # Dummy damping to introduce du, could use a single cell region?
        dterm = Term.new(
            f'dw_zero({dvname}, {duname})', mterm.integral, mterm.region,
            **{dvname : new_variables[dvname],
               duname : new_variables[duname]},
        )
        dterm.setup(allow_derivatives=False)
        dterm.assign_args(new_variables, materials, user=None)
        eq_damping = [dterm]

    new_equations = Equations([Equation('M', Terms(eq_mass), setup=False),
                               Equation('C', Terms(eq_damping), setup=False),
                               Equation('K', Terms(eq_other), setup=False),])

    var_names = {
        'u' : uname, 'du' : duname, 'ddu' : dduname,
        'extra' : set(equations.variables.di.var_names)
        - set([uname, duname, dduname])
    }
    return new_equations, var_names

def gen_multi_vec_packing(di, names, extra_variables=False):
    """
    Return DOF vector (un)packing functions for nlst. For multiphysical
    problems (non-empty `ie` slice for extra variables) the `unpack()` function
    accepts an additional argument `mode` that can be set to 'full' or 'nls'.

    The following DOF ordering must be obeyed:

    - The full DOF vector:

      | ``---iue---|-iv-|-ia-``
      | ``-iu-|-ie-|``
    """
    iu = di.indx[names['u']]
    iv = di.indx[names['du']]
    ia = di.indx[names['ddu']]
    ie = slice(iu.stop, iv.start)
    iue = slice(0, iv.start)
    assert_(iu.stop == di.n_dof[names['u']])
    assert_(iv.start == ie.stop)
    assert_(ia.start == iv.stop)

    if not extra_variables:
        assert_(ie.stop == ie.start)
        n_arg = 3

        def unpack(vec, mode=None):
            return vec[iu], vec[iv], vec[ia]

    else:
        n_arg = 4

        def unpack(vec, mode='full'):
            if mode == 'nls':
                return vec[iu], vec[ie]

            else:
                return vec[iu], vec[ie], vec[iv], vec[ia]

    def pack(*args):
        return nm.concatenate(args)

    indices = dict(indices=(iue, iu, ie, iv, ia),
                   n_dof=di.n_dof_total, n_uedof=ie.stop, n_arg=n_arg)
    unpack.__dict__.update(indices)
    pack.__dict__.update(indices)
    return unpack, pack

def _cache(obj, attr, dep):
    def decorate(fun):
        def new_fun(*args, **kwargs):
            if dep:
                val = getattr(obj, attr)
                if val is None:
                    val = fun(*args, **kwargs)
                    setattr(obj, attr, val)

            else:
                val = fun(*args, **kwargs)

            return val
        return new_fun
    return decorate

class ElastodynamicsBaseTS(TimeSteppingSolver):
    """
    Base class for elastodynamics solvers.

    Assumes block-diagonal matrix in `u`, `v`, `a`.
    """

    _common_parameters = [
        ('t0', 'float', 0.0, False,
         'The initial time.'),
        ('t1', 'float', 1.0, False,
         'The final time.'),
        ('dt', 'float', None, False,
         'The time step. Used if `n_step` is not given.'),
        ('n_step', 'int', 10, False,
         'The number of time steps. Has precedence over `dt`.'),
        ('is_linear', 'bool', False, False,
         'If True, the problem is considered to be linear.'),
        ('has_time_derivatives', 'bool', False, False,
         """If True, the problem equations contain time derivatives of other
            variables besides displacements. In that case the cached constant
            matrices must be cleared on time step changes."""),
        ('var_names', 'dict', None, False,
         """The mapping of variables with keys 'u', 'du', 'ddu' and 'extra',
            and values corresponding to the names of the actual variables.
            See `var_names` returned from :func:`transform_equations_ed()`"""),
    ]

    def __init__(self, conf, nls=None, tsc=None, context=None, **kwargs):
        TimeSteppingSolver.__init__(self, conf, nls=nls, tsc=tsc,
                                    context=context, **kwargs)
        self.conf.quasistatic = False

        if self.tsc is None:
            self.tsc = FixedTSC({})

        if isinstance(self.tsc, FixedTSC):
            # Using TimeStepper instead of VariableTimeStepper ensures the
            # final time is reached "exactly".
            self.ts = TimeStepper.from_conf(self.conf)

        else:
            self.ts = VariableTimeStepper.from_conf(self.conf)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format
        self.verbose = self.conf.verbose
        self.constant_matrices = None
        self.matrix = None

    def get_matrices(self, nls, vec, unpack=None):
        if self.conf.is_linear and self.constant_matrices is not None:
            out = self.constant_matrices

        else:
            aux = nls.fun_grad(vec)

            if unpack is not None:
                iue, iu, ie, iv, ia = unpack.indices
                aux = nls.fun_grad(vec)

                M = aux[ia, ia]
                C = aux[iv, iv]
                K = aux[iue, iue]

            else:
                assert_((len(vec) % 3) == 0)
                i3 = len(vec) // 3

                K = aux[:i3, :i3]
                C = aux[i3:2*i3, i3:2*i3]
                M = aux[2*i3:, 2*i3:]

            out = (M, C, K)

            if self.conf.is_linear:
                M.eliminate_zeros()
                C.eliminate_zeros()
                K.eliminate_zeros()
                self.constant_matrices = (M, C, K)

        return out

    def get_a0(self, nls, u0, e0, v0, unpack):
        iue, iu, ie, iv, ia = unpack.indices

        vec = nm.r_[u0, e0, v0, nm.zeros_like(v0)]
        aux = nls.fun(vec)
        r = aux[iu] + aux[iv] + aux[ia]

        M = self.get_matrices(nls, vec, unpack)[0][iu, iu]
        a0 = nls.lin_solver(-r, mtx=M)
        nls.lin_solver.clear()
        output_array_stats(a0, 'initial acceleration', verbose=self.verbose)
        return a0

    def get_initial_vec(self, nls, vec0, init_fun, prestep_fun, poststep_fun):
        if not set(self.conf.var_names.keys()).issuperset(['ddu', 'du', 'u']):
            raise ValueError(
                "var_names have to contain 'u', 'du', 'ddu' keys!"
            )

        unpack, pack = gen_multi_vec_packing(
            self.di, self.conf.var_names,
            extra_variables=self.get('extra_variables', False),
        )
        self.unpack = unpack
        self.pack = pack

        ts = self.ts
        vec0 = init_fun(ts, vec0)

        output(self.format % (ts.time, ts.step + 1, ts.n_step),
               verbose=self.verbose)
        if ts.step == 0:
            vec0 = prestep_fun(ts, vec0)
            if unpack.n_arg == 4:
                u0, e0, v0, _ = unpack(vec0)

            else:
                u0, v0, _ = unpack(vec0)
                e0 = nm.empty(0, dtype=u0.dtype)

            a0 = self.get_a0(nls, u0, e0, v0, unpack)

            if unpack.n_arg == 4:
                vec = pack(u0, e0, v0, a0)

            else:
                vec = pack(u0, v0, a0)

            vec = poststep_fun(ts, vec)
            ts.advance()

        else:
            vec = vec0

        return vec, unpack, pack

    def clear_lin_solver(self, clear_constant_matrices=True):
        self.nls.lin_solver.clear()
        self.matrix = None
        if clear_constant_matrices:
            self.constant_matrices = None

    @standard_ts_call
    def __call__(self, vec0=None, nls=None, init_fun=None, prestep_fun=None,
                 poststep_fun=None, status=None, **kwargs):
        sig = signature(init_fun)
        if len(sig.parameters) == 3:
            init_fun = partial(init_fun, self)

        sig = signature(poststep_fun)
        if len(sig.parameters) == 3:
            poststep_fun = partial(poststep_fun, self)

        vec, unpack, pack = self.get_initial_vec(
            nls, vec0, init_fun, prestep_fun, poststep_fun)

        ts = self.ts
        dt0 = self.tsc.get_initial_dt(ts, vec, unpack=unpack)
        if not isinstance(self.tsc, FixedTSC):
            ts.set_time_step(dt0, update_time=True)
        while 1:
            output(self.format % (ts.time, ts.step + 1, ts.n_step),
                   verbose=self.verbose)
            # step, time = time step to compute = n+1
            # step-1, time-ts.dt = previous known step data = n
            # adaptivity modifies dt and time.
            while 1:
                # Previous step state q(t_n).
                # Both prestep_fun() and poststep_fun() apply EBCs to vec/vect
                # in case active_only is False.
                # TODO: Generalized alpha: EBCs for current time t_{n+1}. but
                # loads should be applied in the mid-step time t_{n+1-a}.
                vec = prestep_fun(ts, vec)
                vect = self.step(ts, vec, nls, pack, unpack,
                                 prestep_fun=prestep_fun)

                if isinstance(self.tsc, FixedTSC):
                    new_dt = ts.dt
                    break

                new_dt, status = self.tsc(ts, vec, vect, unpack=unpack)
                output('dt:', ts.dt, 'new dt:', new_dt, 'status:', status,
                       verbose=self.verbose)
                if new_dt != ts.dt:
                    self.clear_lin_solver(
                        clear_constant_matrices=self.conf.has_time_derivatives,
                    )

                if status.result == 'accept':
                    break

                ts.set_time_step(new_dt, update_time=True)

            # Current step state q(t_{n+1}).
            vect = poststep_fun(ts, vect)

            if ts.nt >= 1:
                break

            if new_dt != ts.dt:
                ts.set_time_step(new_dt, update_time=False)
            ts.advance()

            vec = vect

        return vec

class VelocityVerletTS(ElastodynamicsBaseTS):
    """
    Solve elastodynamics problems by the explicit velocity-Verlet method.

    The algorithm can be found in [1]_.

    It is mathematically equivalent to the :class:`CentralDifferenceTS` method.
    The current implementation code is essentially the same, as the mid-time
    velocities are not used for anything other than computing the new time
    velocities.

    .. [1] Swope, William C.; H. C. Andersen; P. H. Berens; K. R. Wilson (1
           January 1982). "A computer simulation method for the calculation of
           equilibrium constants for the formation of physical clusters of
           molecules: Application to small water clusters". The Journal of
           Chemical Physics. 76 (1): 648 (Appendix). doi:10.1063/1.442716
    """
    name = 'ts.velocity_verlet'

    _parameters = [
    ] + ElastodynamicsBaseTS._common_parameters

    def create_nlst(self, nls, dt, u0, v0, a0):
        vm = v0 + 0.5 * dt * a0
        u1 = u0 + dt * vm
        def v1(a):
            return vm + 0.5 * dt * a

        nlst = nls.copy()

        def fun(at):
            vec = nm.r_[u1, vm, at]

            aux = nls.fun(vec)

            i3 = len(at)
            rt = aux[:i3] + aux[i3:2*i3] + aux[2*i3:]
            return rt

        @_cache(self, 'matrix', self.conf.is_linear)
        def fun_grad(at):
            vec = None if self.conf.is_linear else nm.r_[u1, vm, at]
            M = self.get_matrices(nls, vec)[0]
            return M

        nlst.fun = fun
        nlst.fun_grad = fun_grad
        nlst.v1 = v1
        nlst.u1 = u1

        return nlst

    def step(self, ts, vec, nls, pack, unpack, **kwargs):
        """
        Solve a single time step.
        """
        dt = ts.dt
        ut, vt, at = unpack(vec)

        nlst = self.create_nlst(nls, dt, ut, vt, at)
        atp = nlst(at)
        vtp = nlst.v1(atp)
        utp = nlst.u1

        vect = pack(utp, vtp, atp)

        return vect

class CentralDifferenceTS(ElastodynamicsBaseTS):
    r"""
    Solve elastodynamics problems by the explicit central difference method.

    It is the same method as obtained by using :class:`NewmarkTS` with
    :math:`\beta = 0`, :math:`\gamma = 1/2`, but uses a simpler code.

    It is also mathematically equivalent to the :class:`VelocityVerletTS`
    method. The current implementation code is essentially the same.

    This solver supports, when used with :class:`RMMSolver
    <sfepy.solvers.ls.RMMSolver>`, the reciprocal mass matrix algorithm, see
    :class:`MassTerm <sfepy.terms.terms_mass.MassTerm>`.
    """
    name = 'ts.central_difference'

    _parameters = [
    ] + ElastodynamicsBaseTS._common_parameters

    def _create_nlst_a(self, nls, dt, ufun, vfun, cc, cache_name, is_rmm=False):
        nlst = nls.copy()

        if is_rmm:
            def fun(at):
                ut = ufun()
                zz = nm.zeros_like(ut)
                vec = nm.r_[ut, zz, zz]

                aux = nls.fun(vec)

                i3 = len(at)
                rt = aux[:i3]
                return rt

            @_cache(self, cache_name, self.conf.is_linear)
            def fun_grad(at):
                M = self.get_matrices(nls, None)[0]

                return M

        else:
            def fun(at):
                vec = nm.r_[ufun(), vfun(at), at]

                aux = nls.fun(vec)

                i3 = len(at)
                rt = aux[:i3] + aux[i3:2*i3] + aux[2*i3:]
                return rt

            @_cache(self, cache_name, self.conf.is_linear)
            def fun_grad(at):
                vec = (None if self.conf.is_linear
                       else nm.r_[ufun(), vfun(at), at])
                M, C = self.get_matrices(nls, vec)[:2]

                Kt = M + cc * C
                return Kt

        nlst.fun = fun
        nlst.fun_grad = fun_grad
        nlst.u = ufun
        nlst.v = vfun

        return nlst

    def create_nlst(self, nls, dt, u0, v0, a0):
        dt2 = dt**2

        def v(a):
            return v0 + dt * 0.5 * (a0 + a)

        def u():
            return u0 + dt * v0 + dt2 * 0.5 * a0

        if isinstance(nls.lin_solver, RMMSolver):
            import scipy.sparse as sps
            class NoNLS(NonlinearSolver):

                def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                             lin_solver=None, iter_hook=None, status=None):
                    vec_r = self.fun(vec_x0)
                    # Dummy all-zero matrix to make standard_call() happy.
                    mtx_a = sps.csr_matrix((vec_r.shape[0], vec_r.shape[0]))
                    return self.lin_solver(-vec_r, mtx=mtx_a)

            nlst = self._create_nlst_a(nls, dt, u, v, 0.5 * dt, 'matrix',
                                       is_rmm=True)
            nlst = NoNLS(Struct(name='nonls', kind='nls.nonls'),
                         fun=nlst.fun, fun_grad=nlst.fun_grad,
                         lin_solver=nlst.lin_solver, iter_hook=nlst.iter_hook,
                         status=nlst.status, context=nlst.context,
                         u=nlst.u, v=nlst.v)
            # nlst.lin_solver.a0 = a0

        else:
            nlst = self._create_nlst_a(nls, dt, u, v, 0.5 * dt, 'matrix')

        return nlst

    def step(self, ts, vec, nls, pack, unpack, **kwargs):
        """
        Solve a single time step.
        """
        dt = ts.dt
        ut, vt, at = unpack(vec)

        nlst = self.create_nlst(nls, dt, ut, vt, at)
        atp = nlst(at)
        vtp = nlst.v(atp)
        utp = nlst.u()

        vect = pack(utp, vtp, atp)

        return vect

class NewmarkTS(ElastodynamicsBaseTS):
    r"""
    Solve elastodynamics problems by the Newmark method.

    The method was introduced in [1]. Common settings [2]:

    ==================== ======== ==== ===== ==========
    name                 kind     beta gamma Omega_crit
    ==================== ======== ==== ===== ==========
    trapezoidal rule:    implicit 1/4  1/2   unconditional
    linear acceleration: implicit 1/6  1/2   :math:`2\sqrt{3}`
    Fox-Goodwin:         implicit 1/12 1/2   :math:`\sqrt{6}`
    central difference:  explicit 0    1/2   2
    ==================== ======== ==== ===== ==========

    All of these methods are 2-order of accuracy.

    [1] Newmark, N. M. (1959) A method of computation for structural dynamics.
    Journal of Engineering Mechanics, ASCE, 85 (EM3) 67-94.

    [2] Arnaud Delaplace, David Ryckelynck: Solvers for Computational Mechanics
    """
    name = 'ts.newmark'
    extra_variables = True

    _parameters = [
        ('beta', 'float', 0.25, False, 'The Newmark method parameter beta.'),
        ('gamma', 'float', 0.5, False, 'The Newmark method parameter gamma.'),
    ] + ElastodynamicsBaseTS._common_parameters

    def create_nlst(self, nls, dt, gamma, beta, u0, e0, v0, a0, pack, unpack):
        dt2 = dt**2
        iue, iu, ie, iv, ia = pack.indices

        cc0 = (1.0 - gamma) * dt
        cc = gamma * dt
        ck0 = (0.5 - beta) * dt2
        ck = beta * dt2

        def v(a):
            return v0 + cc0 * a0 + cc * a

        def u(a):
            return u0 + dt * v0 + ck0 * a0 + ck * a

        if iue == iu:
            def fun(at):
                vec = nm.r_[u(at), v(at), at]

                aux = nls.fun(vec)

                rt = aux[iu] + aux[iv] + aux[ia]
                return rt

            @_cache(self, 'matrix', self.conf.is_linear)
            def fun_grad(at):
                vec = None if self.conf.is_linear else nm.r_[u(at), v(at), at]
                M, C, K = self.get_matrices(nls, vec, unpack)

                Kt = M + cc * C + ck * K
                return Kt

        else: # Extra variables present.
            def fun(aet):
                at = aet[iu]
                vec = nm.r_[u(at), aet[ie], v(at), at]

                aux = nls.fun(vec)

                rt = nm.empty(pack.n_uedof, aux.dtype)
                rt[iu] = aux[iu] + aux[iv] + aux[ia]
                rt[ie] = aux[ie]
                return rt

            @_cache(self, 'matrix', self.conf.is_linear)
            def fun_grad(aet):
                if self.conf.is_linear:
                    vec = None

                else:
                    at = aet[iu]
                    vec = nm.r_[u(at), aet[ie], v(at), at]

                M, C, K = self.get_matrices(nls, vec, unpack)

                Kt = K.copy()
                Kt[:, iu] *= ck
                Kt[iu, iu] += M + cc * C
                return Kt

        nlst = nls.copy()
        nlst.fun = fun
        nlst.fun_grad = fun_grad
        nlst.u = u
        nlst.v = v

        return nlst

    def step(self, ts, vec, nls, pack, unpack, **kwargs):
        """
        Solve a single time step.
        """
        dt = ts.dt
        conf = self.conf
        ut, et, vt, at = unpack(vec)
        nlst = self.create_nlst(nls, dt, conf.gamma, conf.beta, ut, et, vt, at,
                                pack, unpack)

        aetp = nlst(pack(at, et))
        atp, etp = unpack(aetp, mode='nls')
        vtp = nlst.v(atp)
        utp = nlst.u(atp)

        vect = pack(utp, etp, vtp, atp)
        return vect

class GeneralizedAlphaTS(ElastodynamicsBaseTS):
    r"""
    Solve elastodynamics problems by the generalized :math:`\alpha` method.

    - The method was introduced in [1].
    - The method is unconditionally stable provided :math:`\alpha_m \leq
      \alpha_f \leq \frac{1}{2}`, :math:`\beta >= \frac{1}{4} +
      \frac{1}{2}(\alpha_f - \alpha_m)`.
    - The method is second-order accurate provided :math:`\gamma = \frac{1}{2} -
      \alpha_m + \alpha_f`. This is used when `gamma` is ``None``.
    - High frequency dissipation is maximized for :math:`\beta = \frac{1}{4}(1
      - \alpha_m + \alpha_f)^2`. This is used when `beta` is ``None``.
    - The default values of :math:`\alpha_m`, :math:`\alpha_f` (if `alpha_m` or
      `alpha_f`  are ``None``) are based on the user specified high-frequency
      dissipation parameter `rho_inf`.

    Special settings:

    - :math:`\alpha_m = 0` corresponds to the HHT-:math:`\alpha` method.
    - :math:`\alpha_f = 0` corresponds to the WBZ-:math:`\alpha` method.
    - :math:`\alpha_m = 0`, :math:`\alpha_f = 0` produces the Newmark method.

    [1] J. Chung, G.M.Hubert. "A Time Integration Algorithm for Structural
    Dynamics with Improved Numerical Dissipation: The
    Generalized-:math:`\alpha` Method" ASME Journal of Applied Mechanics, 60,
    371:375, 1993.
    """
    name = 'ts.generalized_alpha'

    _parameters = [
        ('rho_inf', 'float', 0.5, False,
         """The spectral radius in the high frequency limit (user specified
            high-frequency dissipation) in [0, 1]:
            1 = no dissipation, 0 = asymptotic annihilation."""),
        ('alpha_m', 'float', None, False,
         r'The parameter :math:`\alpha_m`.'),
        ('alpha_f', 'float', None, False,
         r'The parameter :math:`\alpha_f`.'),
        ('beta', 'float', None, False,
         r'The Newmark-like parameter :math:`\beta`.'),
        ('gamma', 'float', None, False,
         r'The Newmark-like parameter :math:`\gamma`.'),
    ] + ElastodynamicsBaseTS._common_parameters

    def __init__(self, conf, nls=None, tsc=None, context=None, **kwargs):
        ElastodynamicsBaseTS.__init__(self, conf, nls=nls, tsc=tsc,
                                      context=context, **kwargs)
        conf = self.conf

        rho_inf = conf.rho_inf
        alpha_m = get_default(conf.alpha_m,
                              (2.0 * rho_inf - 1.0) / (rho_inf + 1.0))
        alpha_f = get_default(conf.alpha_f, rho_inf / (rho_inf + 1.0))
        beta = get_default(conf.beta, 0.25 * (1.0 - alpha_m + alpha_f)**2)
        gamma = get_default(conf.gamma, 0.5 - alpha_m + alpha_f)
        self.pars = (alpha_m, alpha_f, gamma, beta)

        output('parameters rho_inf, alpha_m, alpha_f, beta, gamma:',
               verbose=self.verbose)
        output(rho_inf, alpha_m, alpha_f, beta, gamma,
               verbose=self.verbose)

    def _create_nlst_a(self, nls, dt, ufun, vfun, afun, cm, cc, ck, cache_name):
        nlst = nls.copy()

        def fun(at):
            vec = nm.r_[ufun(at), vfun(at), afun(at)]

            aux = nls.fun(vec)

            i3 = len(at)
            rt = aux[:i3] + aux[i3:2*i3] + aux[2*i3:]
            return rt

        @_cache(self, cache_name, self.conf.is_linear)
        def fun_grad(at):
            vec = (None if self.conf.is_linear
                   else nm.r_[ufun(at), vfun(at), afun(at)])
            M, C, K = self.get_matrices(nls, vec)

            Kt = cm * M + cc * C + ck * K
            return Kt

        nlst.fun = fun
        nlst.fun_grad = fun_grad

        return nlst

    def create_nlst(self, nls, dt, alpha_m, alpha_f, gamma, beta, u0, v0, a0):
        dt2 = dt**2

        def u1(a):
            return u0 + dt * v0 + dt2 * ((0.5 - beta) * a0 + beta * a)

        def v1(a):
            return v0 + dt * ((1.0 - gamma) * a0 + gamma * a)

        def um(a):
            return (1.0 - alpha_f) * u1(a) + alpha_f * u0

        def vm(a):
            return (1.0 - alpha_f) * v1(a) + alpha_f * v0

        def am(a):
            return (1.0 - alpha_m) * a + alpha_m * a0

        nlst = self._create_nlst_a(nls, dt, um, vm, am,
                                   (1.0 - alpha_m),
                                   (1.0 - alpha_f) * gamma * dt,
                                   (1.0 - alpha_f) * beta * dt2,
                                   'matrix')
        nlst.u = u1
        nlst.v = v1

        return nlst

    def step(self, ts, vec, nls, pack, unpack, **kwargs):
        """
        Solve a single time step.
        """
        dt = ts.dt
        alpha_m, alpha_f, gamma, beta = self.pars
        # Previous step state q(t_n).
        ut, vt, at = unpack(vec)

        nlst = self.create_nlst(nls, dt, alpha_m, alpha_f, gamma, beta,
                                ut, vt, at)

        # Notation: a = \alpha_f, t = t_{n+1}, t - dt = t_n.
        # Set time to t_{n+1-a} = (1 - a) t + a (t - dt) = t - a dt
        ts.set_substep_time(- alpha_f * dt)
        atp = nlst(at)
        # Restore t_{n+1}.
        ts.restore_step_time()

        vtp = nlst.v(atp)
        utp = nlst.u(atp)
        # Current step state q(t_{n+1}).
        vect = pack(utp, vtp, atp)

        return vect

class BatheTS(ElastodynamicsBaseTS):
    """
    Solve elastodynamics problems by the Bathe method.

    The method was introduced in [1].

    [1] Klaus-Juergen Bathe, Conserving energy and momentum in nonlinear
    dynamics: A simple implicit time integration scheme, Computers &
    Structures, Volume 85, Issues 7-8, 2007, Pages 437-445, ISSN 0045-7949,
    https://doi.org/10.1016/j.compstruc.2006.09.004.
    """
    name = 'ts.bathe'

    _parameters = [
    ] + ElastodynamicsBaseTS._common_parameters

    def __init__(self, conf, nls=None, context=None, **kwargs):
        ElastodynamicsBaseTS.__init__(self, conf, nls=nls, context=context,
                                      **kwargs)
        self.matrix1 = None
        self.ls1 = None
        self.ls2 = None

    def _create_nlst_u(self, nls, dt, vfun, afun, cm, cc, cache_name):
        nlst = nls.copy()

        def fun(ut):
            vt = vfun(ut)
            at = afun(vt)
            vec = nm.r_[ut, vt, at]

            aux = nls.fun(vec)

            i3 = len(at)
            rt = aux[:i3] + aux[i3:2*i3] + aux[2*i3:]
            return rt

        @_cache(self, cache_name, self.conf.is_linear)
        def fun_grad(ut):
            if self.conf.is_linear:
                vec = None
            else:
                vt = vfun(ut)
                at = afun(vt)
                vec = nm.r_[ut, vt, at]

            M, C, K = self.get_matrices(nls, vec)

            Kt = cm * M + cc * C + K
            return Kt

        nlst.fun = fun
        nlst.fun_grad = fun_grad
        nlst.v = vfun
        nlst.a = afun

        return nlst

    def create_nlst1(self, nls, dt, u0, v0, a0):
        """
        The first sub-step: the trapezoidal rule.
        """
        dt4 = 4.0 / dt

        def v(u):
            return dt4 * (u - u0) - v0

        def a(v):
            return dt4 * (v - v0) - a0

        nlst = self._create_nlst_u(nls, dt, v, a, dt4 * dt4, dt4, 'matrix1')
        if self.ls1 is None:
            self.ls1 = nls.lin_solver.copy()

        nlst.lin_solver = self.ls1
        return nlst

    def create_nlst2(self, nls, dt, u0, u1, v0, v1):
        """
        The second sub-step: the three-point Euler backward method.
        """
        dt1 = 1.0 / dt
        dt4 = 4.0 * dt1
        dt3 = 3.0 * dt1

        def v(u):
            return dt1 * u0 - dt4 * u1 + dt3 * u

        def a(v):
            return dt1 * v0 - dt4 * v1 + dt3 * v

        nlst = self._create_nlst_u(nls, dt, v, a, dt3 * dt3, dt3, 'matrix')
        if self.ls2 is None:
            self.ls2 = nls.lin_solver.copy()

        nlst.lin_solver = self.ls2
        return nlst

    def clear_lin_solver(self, clear_constant_matrices=True):
        ElastodynamicsBaseTS.clear_lin_solver(
            self, clear_constant_matrices=clear_constant_matrices,
        )
        self.ls1 = self.ls2 = None
        self.matrix1 = None

    def step(self, ts, vec, nls, pack, unpack, prestep_fun):
        """
        Solve a single time step.
        """
        dt = ts.dt
        ut, vt, at = unpack(vec)
        nlst1 = self.create_nlst1(nls, dt, ut, vt, at)
        ut1 = nlst1(ut)
        vt1 = nlst1.v(ut1)
        at1 = nlst1.a(vt1)

        ts.set_substep_time(0.5 * dt)

        vec1 = pack(ut1, vt1, at1)
        vec1 = prestep_fun(ts, vec1)

        nlst2 = self.create_nlst2(nls, dt, ut, ut1, vt, vt1)
        ut2 = nlst2(ut1)
        vt2 = nlst2.v(ut2)
        at2 = nlst2.a(vt2)

        ts.restore_step_time()

        vec2 = pack(ut2, vt2, at2)

        return vec2
