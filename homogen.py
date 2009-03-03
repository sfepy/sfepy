#!/usr/bin/env python
# 11.07.2006, c 
import os.path as op
import shutil
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.fem.evaluate import eval_term_op
import sfepy.homogenization.pfdpm as pfdpm
from sfepy.homogenization.coefs_base import MiniAppBase
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.applications import SimpleApp
from sfepy.base.conf import get_standard_keywords
from sfepy.solvers.generic import solve_evolutionary_op, solve_stationary_op

class Volume( MiniAppBase ):

    def __call__( self, problem = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        volume = eval_term_op( None, self.expression, problem )

        return volume

class PorousMediaHomogenizationApp( SimpleApp ):

    def process_options( options ):
        """Application options setup. Sets default values for missing
        non-compulsory options."""
        get = options.get_default_attr
        
        print_digits = get( 'print_digits', 3 )

        float_format = get( 'float_format', '%8.3e' )
        coef_save_name = get( 'coef_save_name', 'coefs' )

        coefs = get( 'coefs', None, 'missing "coefs" in options!' )
        requirements = get( 'requirements', None,
                            'missing "requirements" in options!' )
        volume = get( 'volume', None, 'missing "volume" in options!' )

        return Struct( **locals() )
    process_options = staticmethod( process_options )

    def __init__( self, conf, options, output_prefix, **kwargs ):
        SimpleApp.__init__( self, conf, options, output_prefix,
                            init_equations = False )

        self.setup_options()
        self.cached_coefs = None

        output_dir = self.problem.output_dir
        shutil.copyfile( conf._filename,
                         op.join( output_dir, op.basename( conf._filename ) ) )

    def setup_options( self ):
        SimpleApp.setup_options( self )
        po = PorousMediaHomogenizationApp.process_options
        self.app_options += po( self.conf.options )
    
    def call( self, ret_all = False ):
        opts = self.app_options
        
        volume = Volume( 'volume', self.problem, opts.volume )()
        output( 'volume: %.2f' % volume )
        
        he = HomogenizationEngine( self.problem, self.options, volume = volume )

        aux = he( ret_all = ret_all)
        if ret_all:
            coefs, dependencies = aux
        else:
            coefs = aux

        coefs = pfdpm.Coefficients( **coefs.to_dict() )
        coefs.volume = volume
        
        prec = nm.get_printoptions()[ 'precision']
        if hasattr( opts, 'print_digits' ):
            nm.set_printoptions( precision = opts.print_digits )
        print coefs
        nm.set_printoptions( precision = prec )
##        pause()

        coef_save_name = op.join( opts.output_dir, opts.coef_save_name )
        coefs.to_file_hdf5( coef_save_name + '.h5' )
        coefs.to_file_txt( coef_save_name + '.txt',
                           self.conf.options.tex_names,
                           opts.float_format )

        if ret_all:
            return coefs, dependencies
        else:
            return coefs

##
# c: 13.06.2008, r: 22.06.2008
def verify_steady_solution( conf, options ):
    if not hasattr( conf, 'equations_steady' ):
        output( 'set "equations_steady" in the input!' )
        return False

    ok = True
    pb = ProblemDefinition.from_conf( conf )
    opts = conf.options
    if hasattr( opts, 'post_process_hook' ) and opts.post_process_hook is not None:
        # User postprocessing.
        pph = getattr( conf.funmod, opts.post_process_hook )
    else:
        pph = None
    state_t, data_t = solve_evolutionary_op( pb, options, post_process_hook = pph )

    pb.set_equations( conf.equations_steady )
    tsc = pb.ts_conf
    ts = pb.get_default_ts( tsc.t0, tsc.t1, tsc.dt, tsc.n_step, tsc.n_step - 1 )
    # To reassemble matrix with new equations.
    pb.set_linear( False )
    state_s, data_s = solve_stationary_op( pb, options, ts = ts,
                                       post_process_hook = pph )

    err = nla.norm( state_s - state_t )
    eps = get_default_attr( conf.options, 'steady_state_error', 1e-8 )
    if err > eps:
        ok = False

    output( 'error: %.2e' % err )

    aux = eval_term_op( state_s, 'dw_volume_wdot.i1.Omega( m.BiotMHI, q, p )',
                      pb, dw_mode = 'vector' )
    print nla.norm( aux )
    aux = eval_term_op( state_t, 'dw_volume_wdot.i1.Omega( m.BiotMHI, q, p )',
                      pb, dw_mode = 'vector' )
    print nla.norm( aux )

    import pylab
    pylab.plot( state_s )
    pylab.plot( state_t )
    pylab.show()

    return ok

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'verify_steady' :
    'verify steady state',
}

##
# c: 11.07.2006, r: 18.06.2008
def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "-v", "--verify-steady",
                       action = "store_true", dest = "verify_steady",
                       default = False, help = help['verify_steady'] )
    (options, args) = parser.parse_args()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    required, other = get_standard_keywords()


    if not options.verify_steady:
        required.remove( 'equations' )

    conf = ProblemConf.from_file( filename_in, required, other )
##     print conf
##     pause()

    if options.verify_steady:
        ok = verify_steady_solution( conf, options )
        if not ok:
            output( 'failed!' )
        else:
            output( 'ok!' )
    else:
        app = PorousMediaHomogenizationApp( conf, options, 'homogen:' )
        opts = conf.options
        if hasattr( opts, 'parametric_hook' ): # Parametric study.
            parametric_hook = getattr( conf, opts.parametric_hook )
            app.parametrize( parametric_hook )
        app()

if __name__ == '__main__':
    main()
