import os.path as op

from sfepy.base.base import *
import sfepy.base.ioutils as io
from sfepy.fem import ProblemDefinition
from sfepy.solvers.generic import solve_direct
from application import Application

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

    def __init__( self, conf, options, output_prefix, **kwargs ):
        """`kwargs` are passed to  ProblemDefinition.from_conf()

        Command-line options have precedence over conf.options."""
        Application.__init__( self, conf, options, output_prefix )
        self.setup_options()

        is_vars = True
        if hasattr(options, 'solve_not') and options.solve_not:
            is_vars = False
        self.problem = ProblemDefinition.from_conf( conf,
                                                    init_variables=is_vars,
                                                    **kwargs )

        self.setup_output_info( self.problem, self.options )

    def setup_options( self ):
        self.app_options = SimpleApp.process_options( self.conf.options )
        funmod = self.conf.funmod
        
        hook = self.app_options.step_hook
        if hook is not None:
            hook = getattr( funmod, hook )
        self.step_hook = hook

        hook = self.app_options.post_process_hook
        if hook is not None:
            hook = getattr( funmod, hook )
        self.post_process_hook = hook

        hook = self.app_options.post_process_hook_final
        if hook is not None:
            hook = getattr( funmod, hook )
        self.post_process_hook_final = hook

        # Override default equations, if use_equations is set.
        if hasattr(self.conf, 'equations'):
            self.conf.equations = getattr(self.conf,
                                          self.app_options.use_equations)

    def setup_output_info( self, problem, options ):
        """Modifies both problem and options!"""
        if options.output_filename_trunk is None:
            ofn_trunk = op.join( self.app_options.output_dir,
                                 io.get_trunk( self.conf.filename_mesh ) )
	    options.output_filename_trunk = ofn_trunk
	else:
            ofn_trunk = options.output_filename_trunk

        problem.ofn_trunk = ofn_trunk
        problem.output_dir = self.app_options.output_dir

        if hasattr(options, 'output_format') \
               and (options.output_format is not None):
            problem.output_format = options.output_format
        else:
            problem.output_format = self.app_options.output_format

    def call( self ):
	out = solve_direct( self.conf, self.options,
                            problem=self.problem,
                            step_hook=self.step_hook,
                            post_process_hook=self.post_process_hook,
                            post_process_hook_final=self.post_process_hook_final)

	return out
