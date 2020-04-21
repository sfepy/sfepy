from __future__ import absolute_import
import os

from sfepy.base.base import output, dict_to_struct, Struct
from sfepy.base.conf import ProblemConf, get_standard_keywords
import sfepy.base.ioutils as io
from sfepy.discrete import Problem
from sfepy.discrete.fem.meshio import MeshIO
from .application import Application

def solve_pde(conf, options=None, status=None, **app_options):
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
    status : dict-like
        The object for storing the solver return status.
    app_options : kwargs
        The keyword arguments that can override application-specific options.
    """
    if not isinstance(conf, ProblemConf):
        required, other = get_standard_keywords()
        conf = ProblemConf.from_file(conf, required, other)

    opts = conf.options = (dict_to_struct(app_options, flag=(1,),
                                          constructor=type(conf.options))
                           + conf.options)

    output_prefix = opts.get('output_prefix', None)
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

    if conf.options.get('evps') is None:
        app = PDESolverApp(conf, options, output_prefix)

    else:
        from .evp_solver_app import EVPSolverApp

        app = EVPSolverApp(conf, options, output_prefix)

    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)

    return app(status=status)

def save_only(conf, save_names, problem=None):
    """
    Save information available prior to setting equations and
    solving them.
    """
    if problem is None:
        problem = Problem.from_conf(conf, init_equations=False)

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
        get = options.get

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
                      format_variant=get('format_variant', None),
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
        """`kwargs` are passed to  Problem.from_conf()

        Command-line options have precedence over conf.options."""
        Application.__init__( self, conf, options, output_prefix )
        self.setup_options()

        is_eqs = init_equations
        if hasattr(options, 'solve_not') and options.solve_not:
            is_eqs = False
        self.problem = Problem.from_conf(conf, init_equations=is_eqs, **kwargs)

        self.setup_output_info( self.problem, self.options )

    def setup_options( self ):
        self.app_options = PDESolverApp.process_options(self.conf.options)

        assign_standard_hooks(self, self.app_options.get, self.conf)

        # Override default equations, if use_equations is set.
        if hasattr(self.conf, 'equations'):
            self.conf.equations = getattr(self.conf,
                                          self.app_options.use_equations)

    def setup_output_info(self, problem, options):
        """Modifies both problem and options!"""
        if options.output_filename_trunk is None:
            if self.conf.get('filename_mesh') is not None:
                filename_mesh = self.conf.filename_mesh
                if isinstance(filename_mesh, MeshIO):
                    ofn_trunk = filename_mesh.get_filename_trunk()

                else:
                    ofn_trunk = io.get_trunk(filename_mesh)

            elif self.conf.get('filename_domain') is not None:
                ofn_trunk = io.get_trunk(self.conf.filename_domain)

            else:
                raise ValueError('missing filename_mesh or filename_domain!')

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
                             format_variant=self.app_options.format_variant,
                             file_per_var=self.app_options.file_per_var,
                             linearization=self.app_options.linearization)

    def call(self, status=None):
        problem = self.problem
        options = self.options

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

        state = problem.solve(
            status=status, save_results=self.app_options.save_results,
            step_hook=self.step_hook,
            post_process_hook=self.post_process_hook,
            post_process_hook_final=self.post_process_hook_final)

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
