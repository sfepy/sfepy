import os.path as op

from sfepy.base.base import *
import sfepy.base.ioutils as io
from sfepy.fem import ProblemDefinition
from sfepy.solvers.generic import solve_direct
from application import Application

def assign_standard_hooks(obj, get, mod):
    """
    Set standard hook function attributes from `mod` to `obj` using the
    `get` function.
    """
    hook = get('step_hook', None)
    if hook is not None:
        hook = getattr(mod, hook)
    obj.step_hook = hook

    hook = get('post_process_hook', None)
    if hook is not None:
        hook = getattr(mod, hook)
    obj.post_process_hook = hook

    hook = get('post_process_hook_final', None)
    if hook is not None:
        hook = getattr(mod, hook)
    obj.post_process_hook_final = hook

class SimpleApp( Application ):

    def process_options( options ):
        """Application options setup. Sets default values for missing
        non-compulsory options."""
        get = options.get_default_attr

        save_results = get( 'save_results', True )
        # Save each variable into a separate file, using the region of its
        # definition only.
        file_per_var = get( 'file_per_var', False )
        output_format = get( 'output_format', 'vtk' )

        output_dir = get( 'output_dir', '.' )
        if not os.path.exists( output_dir ):
            os.makedirs( output_dir )

        # Called after each time step, can do anything, no return value.
        step_hook = get( 'step_hook', None )
        # Called after each time step.
        post_process_hook = get( 'post_process_hook', None )
        # Called after all time steps, or in the stationary case.
        post_process_hook_final = get( 'post_process_hook_final', None )

        use_equations = get('use_equations', 'equations')

        return Struct( **locals() )
    process_options = staticmethod( process_options )

    def __init__(self, conf, options, output_prefix,
		 init_equations=False, **kwargs):
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
        self.app_options = SimpleApp.process_options( self.conf.options )

        assign_standard_hooks(self, self.app_options.get_default_attr,
                              self.conf.funmod)

        # Override default equations, if use_equations is set.
        if hasattr(self.conf, 'equations'):
            self.conf.equations = getattr(self.conf,
                                          self.app_options.use_equations)

    def setup_output_info(self, problem, options):
        """Modifies both problem and options!"""
        if options.output_filename_trunk is None:
            ofn_trunk = op.join(self.app_options.output_dir,
                                 io.get_trunk(self.conf.filename_mesh))
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
                             output_format=output_format)


    def call( self ):
	out = solve_direct( self.conf, self.options,
                            problem=self.problem,
                            step_hook=self.step_hook,
                            post_process_hook=self.post_process_hook,
                            post_process_hook_final=self.post_process_hook_final)

	return out
