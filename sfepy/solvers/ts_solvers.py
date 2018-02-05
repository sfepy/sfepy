"""
Time stepping solvers.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import (get_default, output, assert_,
                             Struct, IndexedStruct, basestr)
from sfepy.linalg.utils import output_array_stats
from sfepy.solvers.solvers import SolverMeta, TimeSteppingSolver
from sfepy.solvers.ts import TimeStepper, VariableTimeStepper
import six

class NewmarkState(Struct):

    def __init__(self, uname, vname, aname):
        Struct.__init__(self, uname=uname, vname=vname, aname=aname)

    def unpack(self, state):
        parts = state.get_parts()
        vvs = state.variables
        u = vvs[self.uname].get_reduced(parts[self.uname])
        v = vvs[self.vname].get_reduced(parts[self.vname])
        a = vvs[self.aname].get_reduced(parts[self.aname])

        return u, v, a

    def pack(self, state, u, v, a):
        vec = nm.r_[u, v, a]

        state.set_reduced(vec)

def gen_multi_vec_packing(size, num):
    assert_((size % num) == 0)
    ii = size // num

    def unpack(vec):
        return [vec[ir:ir+ii] for ir in range(0, size, ii)]

    def pack(*args):
        return nm.concatenate(args)

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

def standard_ts_call(call):
    """
    Decorator handling argument preparation and timing for time-stepping
    solvers.
    """
    import time

    def _standard_ts_call(self, vec0=None, nls=None,
                         init_fun=None, prestep_fun=None, poststep_fun=None,
                         status=None, **kwargs):
        tt = time.clock()

        nls = get_default(nls, self.nls,
                          'nonlinear solver has to be specified!')

        init_fun = get_default(init_fun, lambda ts, vec0: vec0)
        prestep_fun = get_default(prestep_fun, lambda ts, vec: None)
        poststep_fun = get_default(poststep_fun, lambda ts, vec: None)

        result = call(self, vec0=vec0, nls=nls, init_fun=init_fun,
                      prestep_fun=prestep_fun, poststep_fun=poststep_fun,
                      status=status, **kwargs)

        ttt = time.clock() - tt
        if status is not None:
            status['time'] = ttt
            status['n_step'] = self.ts.n_step

        return result

    return _standard_ts_call

class NewmarkTS(TimeSteppingSolver):
    """
    Assumes block-diagonal matrix in `u`, `v`, `a`.

    Common settings [1]:

                                 beta gamma Omega_crit
    trapezoidal rule:    implicit 1/4  1/2  unconditional
    linear acceleration: implicit 1/6  1/2  2\sqrt{3}
    Fox-Goodwin:         implicit 1/12 1/2  \sqrt{6}
    central difference:  explicit 0    1/2  2

    All of these methods are 2-order of accuracy.

    [1] Arnaud Delaplace, David Ryckelynck: Solvers for Computational Mechanics
    """
    name = 'ts.newmark'

    __metaclass__ = SolverMeta

    _parameters = [
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
        ('beta', 'float', 0.25, False, 'The Newmark method parameter beta.'),
        ('gamma', 'float', 0.5, False, 'The Newmark method parameter gamma.'),
    ]

    def __init__(self, conf, nls=None, context=None, **kwargs):
        TimeSteppingSolver.__init__(self, conf, nls=nls, context=context,
                                    **kwargs)
        self.conf.quasistatic = False
        self.ts = TimeStepper.from_conf(self.conf)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format
        self.verbose = self.conf.verbose
        self.constant_matrices = None
        self.matrix = None

    def get_matrices(self, nls, vec):
        if self.conf.is_linear and self.constant_matrices is not None:
            out = self.constant_matrices

        else:
            aux = nls.fun_grad(vec)

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

    def get_a0(self, nls, u0, v0):
        vec = nm.r_[u0, v0, nm.zeros_like(u0)]

        aux = nls.fun(vec)
        i3 = len(u0)
        r = aux[:i3] + aux[i3:2*i3] + aux[2*i3:]

        M = self.get_matrices(nls, vec)[0]
        a0 = nls.lin_solver(-r, mtx=M)
        output_array_stats(a0, 'initial acceleration', verbose=self.verbose)
        return a0

    def create_nlst(self, nls, dt, gamma, beta, u0, v0, a0):
        nlst = nls.copy()

        dt2 = dt**2

        def v(a):
            return v0 + dt * ((1.0 - gamma) * a0 + gamma * a)

        def u(a):
            return u0 + dt * v0 + dt2 * ((0.5 - beta) * a0 + beta * a)

        def fun(at):
            vec = nm.r_[u(at), v(at), at]

            aux = nls.fun(vec)

            i3 = len(at)
            rt = aux[:i3] + aux[i3:2*i3] + aux[2*i3:]
            return rt

        @_cache(self, 'matrix', self.conf.is_linear)
        def fun_grad(at):
            vec = None if self.conf.is_linear else nm.r_[u(at), v(at), at]
            M, C, K = self.get_matrices(nls, vec)

            Kt = M + gamma * dt * C + beta * dt2 * K
            return Kt

        nlst.fun = fun
        nlst.fun_grad = fun_grad
        nlst.u = u
        nlst.v = v

        return nlst

    @standard_ts_call
    def __call__(self, vec0=None, nls=None, init_fun=None, prestep_fun=None,
                 poststep_fun=None, status=None, **kwargs):
        """
        Solve elastodynamics problems by the Newmark method.
        """
        conf = self.conf
        nls = get_default(nls, self.nls)

        ts = self.ts

        vec0 = init_fun(ts, vec0)

        unpack, pack = gen_multi_vec_packing(len(vec0), 3)

        gamma = conf.gamma
        beta = conf.beta

        output(self.format % (ts.time, ts.step + 1, ts.n_step),
               verbose=self.verbose)
        if ts.step == 0:
            prestep_fun(ts, vec0)
            u0, v0, _ = unpack(vec0)

            ut = u0
            vt = v0
            at = self.get_a0(nls, u0, v0)

            vec = pack(ut, vt, at)
            poststep_fun(ts, vec)
            ts.advance()

        else:
            vec = vec0

        for step, time in ts.iter_from(ts.step):
            output(self.format % (time, step + 1, ts.n_step),
                   verbose=self.verbose)
            dt = ts.dt

            prestep_fun(ts, vec)
            ut, vt, at = unpack(vec)

            nlst = self.create_nlst(nls, dt, gamma, beta, ut, vt, at)
            atp = nlst(at)
            vtp = nlst.v(atp)
            utp = nlst.u(atp)

            vect = pack(utp, vtp, atp)
            poststep_fun(ts, vect)

            vec = vect

        return vec

class StationarySolver(TimeSteppingSolver):
    """
    Solver for stationary problems without time stepping.

    This class is provided to have a unified interface of the time stepping
    solvers also for stationary problems.
    """
    name = 'ts.stationary'

    __metaclass__ = SolverMeta

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

        prestep_fun(ts, vec0)

        vec = nls(vec0)

        poststep_fun(ts, vec)

        return vec

def replace_virtuals(deps, pairs):
    out = {}
    for key, val in six.iteritems(deps):
        out[pairs[key]] = val

    return out

class EquationSequenceSolver(TimeSteppingSolver):
    """
    Solver for stationary problems with an equation sequence.
    """
    name = 'ts.equation_sequence'

    __metaclass__ = SolverMeta

    def __init__(self, conf, **kwargs):
        TimeSteppingSolver.__init__(self, conf, ts=None, **kwargs)

    def __call__(self, state0=None, save_results=True, step_hook=None,
                 post_process_hook=None, nls_status=None):
        from sfepy.base.base import invert_dict, get_subdict
        from sfepy.base.resolve_deps import resolve

        problem = self.problem

        if state0 is None:
            state0 = problem.create_state()

        variables = problem.get_variables()
        vtos = variables.get_dual_names()
        vdeps = problem.equations.get_variable_dependencies()
        sdeps = replace_virtuals(vdeps, vtos)

        sorder = resolve(sdeps)

        stov = invert_dict(vtos)
        vorder = [[stov[ii] for ii in block] for block in sorder]

        parts0 = state0.get_parts()
        state = state0.copy()
        solved = []
        for ib, block in enumerate(vorder):
            output('solving for %s...' % sorder[ib])

            subpb = problem.create_subproblem(block, solved)

            subpb.equations.print_terms()

            subpb.time_update()
            substate0 = subpb.create_state()

            vals = get_subdict(parts0, block)
            substate0.set_parts(vals)

            substate = subpb.solve(state0=substate0, nls_status=nls_status)

            state.set_parts(substate.get_parts())

            solved.extend(sorder[ib])
            output('...done')

        if step_hook is not None:
            step_hook(problem, None, state)

        if save_results:
            problem.save_state(problem.get_output_name(), state,
                               post_process_hook=post_process_hook,
                               file_per_var=None)

        yield 0, 0.0, state

    def init_time(self, nls_status=None):
        self.problem.init_solvers(nls_status=nls_status)

def get_initial_state(problem):
    """
    Create a zero state vector and apply initial conditions.
    """
    state = problem.create_state()

    problem.setup_ics()
    state.apply_ic()

    # Initialize variables with history.
    state.init_history()

    return state

def prepare_save_data(ts, conf):
    """
    Given a time stepper configuration, return a list of time steps when the
    state should be saved.
    """
    try:
        save_steps = conf.options.save_steps
    except:
        save_steps = -1

    if save_steps == -1:
        save_steps = ts.n_step

    is_save = nm.linspace(0, ts.n_step - 1, save_steps).astype(nm.int32)
    is_save = nm.unique(is_save)

    return ts.suffix, is_save

def prepare_matrix(problem, state):
    """
    Pre-assemble tangent system matrix.
    """
    problem.update_materials()

    ev = problem.get_evaluator()
    try:
        mtx = ev.eval_tangent_matrix(state(), is_full=True)

    except ValueError:
        output('matrix evaluation failed, giving up...')
        raise

    return mtx

def make_implicit_step(ts, state0, problem, nls_status=None):
    """
    Make a step of an implicit time stepping solver.
    """
    if ts.step == 0:
        state0.apply_ebc()
        state = state0.copy(deep=True)

        if not ts.is_quasistatic:
            ev = problem.get_evaluator()
            try:
                vec_r = ev.eval_residual(state(), is_full=True)
            except ValueError:
                output('initial residual evaluation failed, giving up...')
                raise
            else:
                err = nm.linalg.norm(vec_r)
                output('initial residual: %e' % err)

        if problem.is_linear():
            mtx = prepare_matrix(problem, state)
            problem.try_presolve(mtx)

        if ts.is_quasistatic:
            # Ordinary solve.
            state = problem.solve(state0=state0, nls_status=nls_status)

    else:
        problem.time_update(ts)
        state = problem.solve(state0=state0, nls_status=nls_status)

    return state

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

class SimpleTimeSteppingSolver(TimeSteppingSolver):
    """
    Implicit time stepping solver with a fixed time step.
    """
    name = 'ts.simple'

    __metaclass__ = SolverMeta

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
            prestep_fun(ts, vec0)

            vec = self.solve_step0(nls, vec0)

            poststep_fun(ts, vec)
            ts.advance()

        else:
            vec = vec0

        for step, time in ts.iter_from(ts.step):
            self.output_step_info(ts)

            prestep_fun(ts, vec)

            vect = self.solve_step(ts, nls, vec, prestep_fun)

            poststep_fun(ts, vect)

            vec = vect

        return vec

class AdaptiveTimeSteppingSolver(SimpleTimeSteppingSolver):
    """
    Implicit time stepping solver with an adaptive time step.

    Either the built-in or user supplied function can be used to adapt the time
    step.
    """
    name = 'ts.adaptive'

    __metaclass__ = SolverMeta

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
        status = IndexedStruct(n_iter=0, condition=0)
        while 1:
            vect = nls(vec, status=status)

            is_break = self.adapt_time_step(ts, status, self.adt, self.context,
                                            verbose=self.verbose)

            if is_break:
                break

            prestep_fun(ts, vec)

        return vect

    def output_step_info(self, ts):
        output(self.format % (ts.time, ts.dt, self.adt.wait,
                              ts.step + 1, ts.n_step),
               verbose=self.verbose)
