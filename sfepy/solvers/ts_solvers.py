"""
Time stepping solvers.
"""
import numpy as nm

from sfepy.base.base import output, Struct, IndexedStruct, basestr
from sfepy.solvers.solvers import make_get_conf, TimeSteppingSolver
from sfepy.discrete.mass_operator import MassOperator
from sfepy.solvers.ts import TimeStepper, VariableTimeStepper

class StationarySolver(TimeSteppingSolver):
    """
    Solver for stationary problems without time stepping.

    This class is provided to have a unified interface of the time stepping
    solvers also for stationary problems.
    """
    name = 'ts.stationary'

    def __init__(self, conf, **kwargs):
        TimeSteppingSolver.__init__(self, conf, ts=None, **kwargs)

    def __call__(self, state0=None, save_results=True, step_hook=None,
                 post_process_hook=None, nls_status=None):
        problem = self.problem

        problem.time_update()

        state = problem.solve(state0=state0, nls_status=nls_status)

        if step_hook is not None:
            step_hook(problem, None, state)

        if save_results:
            problem.save_state(problem.get_output_name(), state,
                               post_process_hook=post_process_hook,
                               file_per_var=None)

        return state

def replace_virtuals(deps, pairs):
    out = {}
    for key, val in deps.iteritems():
        out[pairs[key]] = val

    return out

class EquationSequenceSolver(TimeSteppingSolver):
    """
    Solver for stationary problems with an equation sequence.
    """
    name = 'ts.equation_sequence'

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

        return state

def get_initial_state(problem):
    """
    Create a zero state vector and apply initial conditions.
    """
    state = problem.create_state()

    problem.setup_ic()
    state.apply_ic()

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
    problem.time_update(ts)

    if ts.step == 0:
        state0.apply_ebc()
        state = state0.copy(deep=True)

        if not ts.is_quasistatic:
            problem.init_time(ts)

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

        else:
            mtx = None

        # Initialize solvers (and possibly presolve the matrix).
        presolve = mtx is not None
        problem.init_solvers(nls_status=nls_status, mtx=mtx, presolve=presolve)

        # Initialize variables with history.
        state0.init_history()
        if ts.is_quasistatic:
            # Ordinary solve.
            state = problem.solve(state0=state0, nls_status=nls_status)

    else:
        if (ts.step == 1) and ts.is_quasistatic and problem.is_linear():
            mtx = prepare_matrix(problem, state0)
            problem.init_solvers(nls_status=nls_status, mtx=mtx)

        state = problem.solve(state0=state0, nls_status=nls_status)

    return state

def make_explicit_step(ts, state0, problem, mass, nls_status=None):
    """
    Make a step of an explicit time stepping solver.
    """
    problem.time_update(ts)

    if ts.step == 0:
        state0.apply_ebc()
        state = state0.copy(deep=True)

        problem.init_time(ts)

        # Initialize variables with history.
        state0.init_history()

    ev = problem.get_evaluator()
    try:
        vec_r = ev.eval_residual(state0(), is_full=True)
    except ValueError:
        output('residual evaluation failed, giving up...')
        raise
    else:
        err = nm.linalg.norm(vec_r)
        output('residual: %e' % err)

    if ts.step > 0:
        variables = problem.get_variables()
        vec_rf = variables.make_full_vec(vec_r, force_value=0.0)

        rhs = -ts.dt * vec_rf + mass.action(state0())

        vec = mass.inverse_action(rhs)

        state = state0.copy(preserve_caches=True)
        state.set_full(vec)
        state.apply_ebc()

    return state

def get_min_dt(adt):
    red = adt.red
    while red >= adt.red_max:
        red *= adt.red_factor

    dt = adt.dt0 * red

    return dt

def adapt_time_step(ts, status, adt, problem=None):
    """
    Adapt the time step of `ts` according to the exit status of the
    nonlinear solver.

    The time step dt is reduced, if the nonlinear solver did not converge. If it
    converged in less then a specified number of iterations for several time
    steps, the time step is increased. This is governed by the following
    parameters:

    - red_factor : time step reduction factor
    - red_max : maximum time step reduction factor
    - inc_factor : time step increase factor
    - inc_on_iter : increase time step if the nonlinear solver converged in
      less than this amount of iterations...
    - inc_wait : ...for this number of consecutive time steps

    Parameters
    ----------
    ts : VariableTimeStepper instance
        The time stepper.
    status : IndexedStruct instance
        The nonlinear solver exit status.
    adt : Struct instance
        The adaptivity parameters of the time solver:
    problem : ProblemDefinition instance, optional
        This canbe used in user-defined adaptivity functions. Not used here.

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
                    output('+++++ new time step: %e +++++' % ts.dt)
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
            output('----- new time step: %e -----' % ts.dt)
            adt.wait = 0

    return is_break

class SimpleTimeSteppingSolver(TimeSteppingSolver):
    """
    Implicit time stepping solver with a fixed time step.
    """
    name = 'ts.simple'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Process configuration options.
        """
        get = make_get_conf(conf, kwargs)
        common = TimeSteppingSolver.process_conf(conf)

        return Struct(t0=get('t0', 0.0),
                      t1=get('t1', 1.0),
                      dt=get('dt', None),
                      n_step=get('n_step', 10),
                      quasistatic=get('quasistatic', False)) + common

    def __init__(self, conf, **kwargs):
        TimeSteppingSolver.__init__(self, conf, **kwargs)

        self.ts = TimeStepper.from_conf(self.conf)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format

    def __call__(self, state0=None, save_results=True, step_hook=None,
                 post_process_hook=None, nls_status=None):
        """
        Solve the time-dependent problem.
        """
        problem = self.problem
        ts = self.ts

        suffix, is_save = prepare_save_data(ts, problem.conf)

        if state0 is None:
            state0 = get_initial_state(problem)

        ii = 0
        for step, time in ts:
            output(self.format % (time, step + 1, ts.n_step))

            state = self.solve_step(ts, state0, nls_status=nls_status)
            state0 = state.copy(deep=True)

            if step_hook is not None:
                step_hook(problem, ts, state)

            if save_results and (is_save[ii] == ts.step):
                filename = problem.get_output_name(suffix=suffix % ts.step)
                problem.save_state(filename, state,
                                   post_process_hook=post_process_hook,
                                   file_per_var=None,
                                   ts=ts)
                ii += 1

            problem.advance(ts)

        return state

    def solve_step(self, ts, state0, nls_status=None):
        """
        Solve a single time step.
        """
        state = make_implicit_step(ts, state0, self.problem,
                                   nls_status=nls_status)

        return state

class ExplicitTimeSteppingSolver(SimpleTimeSteppingSolver):
    """
    Explicit time stepping solver with a fixed time step.
    """
    name = 'ts.explicit'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Process configuration options.
        """
        get = make_get_conf(conf, kwargs)
        common = SimpleTimeSteppingSolver.process_conf(conf, kwargs)

        return Struct(mass=get('mass', None,
                               'missing "mass" in options!'),
                      lumped=get('lumped', False)) + common

    def __init__(self, conf, **kwargs):
        SimpleTimeSteppingSolver.__init__(self, conf, **kwargs)

        self.mass = MassOperator(self.problem, self.conf)

    def solve_step(self, ts, state0, nls_status=None):
        """
        Solve a single time step.
        """
        state = make_explicit_step(ts, state0, self.problem, self.mass,
                                   nls_status=nls_status)

        return state

class AdaptiveTimeSteppingSolver(SimpleTimeSteppingSolver):
    """
    Implicit time stepping solver with an adaptive time step.

    Either the built-in or user supplied function can be used to adapt the time
    step.
    """
    name = 'ts.adaptive'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Process configuration options.
        """
        get = make_get_conf(conf, kwargs)
        common = SimpleTimeSteppingSolver.process_conf(conf, kwargs)

        adt = Struct(red_factor=get('dt_red_factor', 0.2),
                     red_max=get('dt_red_max', 1e-3),
                     inc_factor=get('dt_inc_factor', 1.25),
                     inc_on_iter=get('dt_inc_on_iter', 4),
                     inc_wait=get('dt_inc_wait', 5),
                     red=1.0, wait=0, dt0=0.0)

        return Struct(adapt_fun=get('adapt_fun', adapt_time_step),
                      adt=adt) + common

    def __init__(self, conf, **kwargs):
        TimeSteppingSolver.__init__(self, conf, **kwargs)

        self.ts = VariableTimeStepper.from_conf(self.conf)

        self.adt = adt = self.conf.adt
        adt.dt0 = self.ts.get_default_time_step()
        self.ts.set_n_digit_from_min_dt(get_min_dt(adt))

        self.format = '====== time %e (dt %e, wait %d, step %d of %d) ====='

        if isinstance(self.conf.adapt_fun, basestr):
            self.adapt_time_step = self.problem.functions[self.conf.adapt_fun]

        else:
            self.adapt_time_step = self.conf.adapt_fun

    def __call__(self, state0=None, save_results=True, step_hook=None,
                 post_process_hook=None, nls_status=None):
        """
        Solve the time-dependent problem.
        """
        problem = self.problem
        ts = self.ts

        if state0 is None:
            state0 = get_initial_state(problem)

        ii = 0
        for step, time in ts:
            output(self.format % (time, ts.dt, self.adt.wait,
                                  step + 1, ts.n_step))

            state = self.solve_step(ts, state0, nls_status=nls_status)
            state0 = state.copy(deep=True)

            if step_hook is not None:
                step_hook(problem, ts, state)

            if save_results:
                filename = problem.get_output_name(suffix=ts.suffix % ts.step)
                problem.save_state(filename, state,
                                   post_process_hook=post_process_hook,
                                   file_per_var=None,
                                   ts=ts)
                ii += 1

            problem.advance(ts)

        return state

    def solve_step(self, ts, state0, nls_status=None):
        """
        Solve a single time step.
        """
        status = IndexedStruct(n_iter=0, condition=0)
        while 1:
            state = make_implicit_step(ts, state0, self.problem,
                                       nls_status=status)

            is_break = self.adapt_time_step(ts, status, self.adt, self.problem)
            if is_break:
                break

        if nls_status is not None:
            nls_status.update(status)

        return state
