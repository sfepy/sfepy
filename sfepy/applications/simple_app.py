import os.path as op

from sfepy.base.base import *
import sfepy.base.ioutils as io
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.solvers.generic import solve_direct
from application import Application

class SimpleApp( Application ):

    def __init__( self, conf, options, output_prefix ):
        Application.__init__( self, conf, options, output_prefix )

        opts = conf.options
        output_dir = get_default_attr( opts, 'output_dir', '.' )
        if not os.path.exists( output_dir ):
            os.makedirs( output_dir )

        if options.output_filename_trunk is None:
            ofn_trunk = op.join( output_dir, io.get_trunk( conf.filename_mesh ) )
	    options.output_filename_trunk = ofn_trunk
	else:
            ofn_trunk = options.output_filename_trunk

        self.problem = ProblemDefinition.from_conf( conf )
        self.problem.ofn_trunk = ofn_trunk
        self.problem.output_dir = output_dir

    def call( self ):
	dpb, vec_dp, data = solve_direct( self.conf, self.options,
					  problem = self.problem )

	opts = self.conf.options
	if hasattr( opts, 'post_process_hook_final' ): # User postprocessing.
	    hook = getattr( conf, opts.post_process_hook_final )
	    hook( dpb, vec_dp, data )

	return dpb, vec_dp, data
