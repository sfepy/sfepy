import os
import numpy as nm

from sfepy.base.base import output, dict_to_struct, Struct
from sfepy.base.conf import ProblemConf, get_standard_keywords
import sfepy.base.ioutils as io
from sfepy.fem import ProblemDefinition
from sfepy.fem.meshio import MeshIO
from sfepy.fem.mass_operator import MassOperator
from application import Application

def solve_pde(conf, options=None, nls_status=None, **app_options):
    """
    Solve a system of partial differential equations (PDEs).

    This function is a convenience wrapper that creates and runs an instance of
    :class:`PDESolverApp`.

    Parameters
    ----------
    conf : str or ProblemConf instance
        Either the name of the problem description file defining the PDEs,
        or directly the ProblemConf instance.
    options : options
        The command-line options.
    nls_status : dict-like
        The object for storing the nonlinear solver return status.
    app_options : kwargs
        The keyword arguments that can override application-specific options.
    """
    if not isinstance(conf, ProblemConf):
        required, other = get_standard_keywords()
        conf = ProblemConf.from_file(conf, required, other)

    opts = conf.options = dict_to_struct(app_options, flag=(1,)) + conf.options

    output_prefix = opts.get_default_attr('output_prefix', None)
    if output_prefix is None:
        output_prefix = output.prefix

    if options is None:
        options = Struct(output_filename_trunk=None,
                         save_ebc=False,
                         save_ebc_nodes=False,
                         save_regions=False,
                         save_field_meshes=False,
                         save_regions_as_groups=False,
                         solve_not=False)

    app = PDESolverApp(conf, options, output_prefix)
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)

    return app(nls_status=nls_status)

def save_only(conf, save_names, problem=None):
    """
    Save information available prior to setting equations and
    solving them.
    """
    if problem is None:
        problem = ProblemDefinition.from_conf(conf, init_equations=False)

    if save_names.regions is not None:
        problem.save_regions(save_names.regions)

    if save_names.regions_as_groups is not None:
        problem.save_regions_as_groups(save_names.regions_as_groups)

    if save_names.field_meshes is not None:
        problem.save_field_meshes(save_names.field_meshes)

    if save_names.ebc is not None:
        problem.save_ebc(save_names.ebc, force=False)

    if save_names.ebc_nodes is not None:
        problem.save_ebc(save_names.ebc_nodes, force=True)

def solve_stationary(problem, save_results=True, ts=None,
                     post_process_hook=None,
                     nls_status=None):
    """
    Solve a stationary problem.
    """
    if ts is None:
        try:
            ts = problem.get_time_solver().ts
        except ValueError:
            pass

    problem.time_update(ts)
    state = problem.solve(nls_status=nls_status)

    if save_results:
        problem.save_state(problem.get_output_name(), state,
                           post_process_hook=post_process_hook,
                           file_per_var=None)

    return state

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

def prepare_save_data(ts, conf):
    try:
        save_steps = conf.options.save_steps
    except:
        save_steps = -1

    if save_steps == -1:
        save_steps = ts.n_step

    is_save = nm.linspace(0, ts.n_step - 1, save_steps).astype(nm.int32)
    is_save = nm.unique(is_save)

    return ts.suffix, is_save

def make_implicit_step(ts, state0, problem, nls_status=None):
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
            state = problem.solve(state0=state0)

    else:
        if (ts.step == 1) and ts.is_quasistatic and problem.is_linear():
            mtx = prepare_matrix(problem, state0)
            problem.init_solvers(nls_status=nls_status, mtx=mtx)

        state = problem.solve(state0=state0)

    return state

def make_explicit_step(ts, state0, problem, mass, nls_status=None):
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

def solve_evolutionary(problem, time_solver=None,
                       save_results=True, return_history=False,
                       step_hook=None, post_process_hook=None,
                       nls_status=None):
    """
    Solve an evolutionary problem.

    TODO:  return_history
    """
    if time_solver is None:
        time_solver = problem.get_time_solver(step_fun=make_implicit_step,
                                              step_args=(problem, nls_status))

    suffix, is_save = prepare_save_data(time_solver.ts, problem.conf)

    state0 = problem.create_state()
    problem.setup_ic()
    state0.apply_ic()

    ii = 0
    for ts, state in time_solver(state0):

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

def assign_standard_hooks(obj, get, conf):
    """
    Set standard hook function attributes from `conf` to `obj` using the
    `get` function.
    """
    hook_names = ['step_hook', 'post_process_hook',
                  'post_process_hook_final', 'pre_process_hook']
    for hook_name in hook_names:
        setattr(obj, hook_name, conf.get_function(get(hook_name, None)))

class PDESolverApp(Application):

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.
        """
        get = options.get_default_attr

        output_dir = get('output_dir', '.')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        return Struct(save_results=get('save_results', True),
                      # Save each variable into a separate file, using
                      # the region of its definition only.
                      linearization=dict_to_struct(get('linearization',
                                                       {'kind' : 'strip'})),
                      file_per_var=get('file_per_var', False),
                      output_format=get('output_format', 'vtk'),
                      output_dir=output_dir,
                      # Called after each time step, can do anything, no
                      # return value.
                      step_hook=get('step_hook', None),
                      # Called after each time step.
                      post_process_hook=get('post_process_hook', None),
                      # Called after all time steps, or in the
                      # stationary case.
                      post_process_hook_final=get('post_process_hook_final',
                                                  None),
                      # Called in init process.
                      pre_process_hook=get('pre_process_hook', None),
                      use_equations=get('use_equations', 'equations'))

    def __init__(self, conf, options, output_prefix,
                 init_equations=True, **kwargs):
        """`kwargs` are passed to  ProblemDefinition.from_conf()

        Command-line options have precedence over conf.options."""
        Application.__init__( self, conf, options, output_prefix )
        self.setup_options()

        is_eqs = init_equations
        if hasattr(options, 'solve_not') and options.solve_not:
            is_eqs = False
        self.problem = ProblemDefinition.from_conf(conf,
                                                   init_equations=is_eqs,
                                                   **kwargs)

        self.setup_output_info( self.problem, self.options )

    def setup_options( self ):
        self.app_options = PDESolverApp.process_options(self.conf.options)

        assign_standard_hooks(self, self.app_options.get_default_attr,
                              self.conf)

        # Override default equations, if use_equations is set.
        if hasattr(self.conf, 'equations'):
            self.conf.equations = getattr(self.conf,
                                          self.app_options.use_equations)

    def setup_output_info(self, problem, options):
        """Modifies both problem and options!"""
        if options.output_filename_trunk is None:
            filename_mesh = self.conf.filename_mesh
            if isinstance(filename_mesh, MeshIO):
                ofn_trunk = filename_mesh.get_filename_trunk()

            else:
                ofn_trunk = io.get_trunk(filename_mesh)

            options.output_filename_trunk = ofn_trunk

        else:
            ofn_trunk = options.output_filename_trunk

        if hasattr(options, 'output_format') \
               and (options.output_format is not None):
            output_format = options.output_format

        else:
            output_format = self.app_options.output_format

        problem.setup_output(output_filename_trunk=ofn_trunk,
                             output_dir=self.app_options.output_dir,
                             output_format=output_format,
                             file_per_var=self.app_options.file_per_var,
                             linearization=self.app_options.linearization)

    def call(self, nls_status=None):
        problem = self.problem
        options = self.options
        opts = self.app_options

        if self.pre_process_hook is not None: # User pre_processing.
            self.pre_process_hook(problem)

        ofn_trunk = problem.ofn_trunk
        self.save_names = Struct(ebc=ofn_trunk + '_ebc.vtk'
                                 if options.save_ebc else None,

                                 ebc_nodes=ofn_trunk + '_ebc_nodes.vtk'
                                 if options.save_ebc_nodes else None,

                                 regions=ofn_trunk + '_region'
                                 if options.save_regions else None,

                                 regions_as_groups=ofn_trunk + '_regions'
                                 if options.save_regions_as_groups else None,

                                 field_meshes=ofn_trunk + '_field'
                                 if options.save_field_meshes else None)

        if any(self.save_names.to_dict().values()):
            save_only(self.conf, self.save_names, problem=problem)

        if options.solve_not:
            return None, None, None

        if hasattr(self.conf.options, 'ts'):
            time_solver = problem.get_time_solver()
            if time_solver.name == 'ts.simple': # Implicit time stepping.
                time_solver.set_step_fun(make_implicit_step,
                                         (problem, nls_status))

            else: # Explicit time stepping.
                mass = MassOperator(problem, time_solver.conf)

                time_solver.set_step_fun(make_explicit_step,
                                         (problem, mass, nls_status))

            state = solve_evolutionary(problem, time_solver,
                                       save_results=opts.save_results,
                                       step_hook=self.step_hook,
                                       post_process_hook=self.post_process_hook)

        else: # Stationary problem.
            state = solve_stationary(problem,
                                     save_results=opts.save_results,
                                     post_process_hook=self.post_process_hook,
                                     nls_status=nls_status)

        if self.post_process_hook_final is not None: # User postprocessing.
            self.post_process_hook_final(problem, state)

        return problem, state

    def save_dict(self, filename, data):
        """
        Utility function to save a dictionary `data` to a HDF5 file
        `filename`.
        """
        io.write_dict_hdf5(filename, data)

    def load_dict(self, filename):
        """
        Utility function to load a dictionary `data` from a HDF5 file
        `filename`.
        """
        data = io.read_dict_hdf5(filename)

        return data
