#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
import shutil
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.la import eig
from sfepy.fem import eval_term_op, ProblemDefinition
from sfepy.fem.evaluate import assemble_by_blocks
from sfepy.homogenization.phono import process_options, get_method,\
     transform_plot_data, plot_logs, plot_gaps, detect_band_gaps, compute_cat
from sfepy.homogenization.utils import create_pis
from sfepy.homogenization.coefs import CorrectorsRS, CoefE
from sfepy.applications import SimpleApp
from sfepy.solvers import Solver
from sfepy.base.plotutils import pylab

def make_save_hook( base_name, post_process_hook = None, file_per_var = None ):
    def save_phono_correctors( state, problem, ir, ic ):
        get_state = problem.variables.get_state_part_view

        problem.save_state( (base_name % (ir, ic)) + '.vtk', state,
                            post_process_hook = post_process_hook,
                            file_per_var = file_per_var )
    return save_phono_correctors

class AcousticBandGapsApp( SimpleApp ):

    def __init__( self, conf, options, output_prefix, **kwargs ):
        SimpleApp.__init__( self, conf, options, output_prefix,
                            init_equations = False )

        opts = conf.options
        post_process_hook = get_default_attr( opts, 'post_process_hook', None )
        if post_process_hook is not None:
            post_process_hook = getattr( conf.funmod, post_process_hook )

        self.post_process_hook = post_process_hook

        output_dir = self.problem.output_dir
        shutil.copyfile( conf._filename,
                         op.join( output_dir, op.basename( conf._filename ) ) )
        
    def call( self ):
        options = self.options
        evp = self.solve_eigen_problem()

        if options.detect_band_gaps:
            bg = detect_band_gaps( self.problem, evp.eigs, evp.eig_vectors,
                                   self.conf.options, self.conf.funmod )

            if options.plot:
                plot_range, tlogs = transform_plot_data( bg.logs,
                                                         bg.opts.plot_tranform,
                                                         self.conf.funmod )

                plot_rsc = bg.opts.plot_rsc
                plot_opts =  bg.opts.plot_options
                
                pylab.rcParams.update( plot_rsc['params'] )

                fig = plot_gaps( 1, plot_rsc, bg.gaps, bg.kinds,
                                 bg.freq_range_margins, plot_range,
                                 clear = True )
                fig = plot_logs( 1, plot_rsc, tlogs, bg.valid[bg.eig_range],
                                 bg.freq_range_initial,
                                 plot_range, bg.opts.squared,
                                 show = plot_opts['show'],
                                 show_legend = plot_opts['legend'],
                                 new_axes = True )

                fig_name = bg.opts.fig_name
                if fig_name is not None:
                    fig.savefig( fig_name )

        elif options.analyze_dispersion:
            dim = self.problem.domain.mesh.dim

            opts = self.conf.options
            iw_dir = nm.array( opts.incident_wave_dir, dtype = nm.float64 )
            assert_( dim == iw_dir.shape[0] )

            iw_dir = iw_dir / nla.norm( iw_dir )

            mtx_d = self.eval_coef_e( )
            print mtx_d

            christoffel = compute_cat( mtx_d, iw_dir )
            print iw_dir
            print christoffel
            
            bg = detect_band_gaps( self.problem, evp.eigs, evp.eig_vectors,
                                   self.conf.options, self.conf.funmod,
                                   christoffel = christoffel )
            if options.plot:
                # This will change...
                plot_range, tlogs = transform_plot_data( bg.logs,
                                                         None,
                                                         self.conf.funmod )

                plot_rsc = bg.opts.plot_rsc
                plot_opts =  bg.opts.plot_options
                
                pylab.rcParams.update( plot_rsc['params'] )

                fig = plot_gaps( 1, plot_rsc, bg.gaps, bg.kinds,
                                 bg.freq_range_margins, plot_range,
                                 clear = True )
                fig = plot_logs( 1, plot_rsc, tlogs, bg.valid[bg.eig_range],
                                 bg.freq_range_initial,
                                 plot_range, bg.opts.squared,
                                 show = plot_opts['show'],
                                 show_legend = plot_opts['legend'],
                                 new_axes = True )

                fig_name = bg.opts.fig_name
                if fig_name is not None:
                    fig.savefig( fig_name )

        return evp, bg
    
    def solve_eigen_problem( self, ofn_trunk = None, post_process_hook = None ):

        problem = self.problem
        ofn_trunk = get_default( ofn_trunk, problem.ofn_trunk,
                                 'output file name trunk missing!' )
        post_process_hook = get_default( post_process_hook,
                                         self.post_process_hook )
 
        conf = self.conf
        
        eig_problem = get_default_attr( conf.options, 'eig_problem', 'simple' )
        if eig_problem == 'simple':
            problem.set_equations( conf.equations )

            dim = problem.domain.mesh.dim
            problem.time_update()

            dummy = problem.create_state_vector()
            mtx_a = eval_term_op( dummy, conf.equations['lhs'], problem,
                                  dw_mode = 'matrix',
                                  tangent_matrix = problem.mtx_a )

            mtx_m = eval_term_op( dummy, conf.equations['rhs'], problem,
                                  dw_mode = 'matrix',
                                  tangent_matrix = problem.mtx_a.copy() )

        elif eig_problem == 'schur':
            # A = K + B^T D^{-1} B.
            mtx = assemble_by_blocks( conf.equations, self.problem )
            problem.set_equations( conf.equations )
            problem.time_update()

            ls = Solver.any_from_conf( problem.ls_conf,
                                       presolve = True, mtx = mtx['D'] )

            mtx_b, mtx_m = mtx['B'], mtx['M']
            mtx_dib = nm.empty( mtx_b.shape, dtype = mtx_b.dtype )
            for ic in xrange( mtx_b.shape[1] ):
                mtx_dib[:,ic] = ls( mtx_b[:,ic].toarray().squeeze() )
            mtx_a = mtx['K'] + mtx_b.T * mtx_dib

        else:
            raise NotImplementedError

##     from sfepy.base.plotutils import spy, pylab
##     spy( mtx_b, eps = 1e-12 )
##     pylab.show()
##     mtx_a.save( 'a.txt', format='%d %d %.12f\n' )
##     mtx_b.save( 'b.txt', format='%d %d %.12f\n' )
##     pause()

        output( 'computing resonance frequencies...' )
        tt = [0]

        if isinstance( mtx_a, sc.sparse.spmatrix ):
            mtx_a = mtx_a.toarray()
        if isinstance( mtx_m, sc.sparse.spmatrix ):
            mtx_m = mtx_m.toarray()

        eigs, mtx_s_phi = eig( mtx_a, mtx_m, return_time = tt,
                               method = get_method( conf.options ) )
        output( '...done in %.2f s' % tt[0] )
        output( eigs )
        output( 'number of frequencies: %d' % eigs.shape[0] )

        try:
            assert_( nm.isfinite( eigs ).all() )
        except ValueError:
            debug()

        # B-orthogonality check.
##         print nm.dot( mtx_s_phi[:,5], nm.dot( mtx_m, mtx_s_phi[:,5] ) )
##         print nm.dot( mtx_s_phi[:,5], nm.dot( mtx_m, mtx_s_phi[:,0] ) )
##         debug()

        n_eigs = eigs.shape[0]
        opts = process_options( conf.options, n_eigs )

        mtx_phi = nm.empty( (problem.variables.di.ptr[-1], mtx_s_phi.shape[1]),
                           dtype = nm.float64 )

        make_full = problem.variables.make_full_vec
        if eig_problem == 'simple':
            for ii in xrange( n_eigs ):
                mtx_phi[:,ii] = make_full( mtx_s_phi[:,ii] )
            eig_vectors = mtx_phi
            
        elif eig_problem == 'schur':
            # Update also eliminated variables.
            schur = conf.options.schur
            primary_var = schur['primary_var']
            eliminated_var = schur['eliminated_var']

            mtx_s_phi_schur = - sc.dot( mtx_dib, mtx_s_phi )
            aux = nm.empty( (problem.variables.adi.ptr[-1],),
                            dtype = nm.float64 )
            set = problem.variables.set_state_part
            for ii in xrange( n_eigs ):
                set( aux, mtx_s_phi[:,ii], primary_var, stripped = True )
                set( aux, mtx_s_phi_schur[:,ii], eliminated_var,
                     stripped = True )

                mtx_phi[:,ii] = make_full( aux )

            indx = problem.variables.get_indx( primary_var )
            eig_vectors = mtx_phi[indx,:]

        out = {}
        for ii in xrange( n_eigs ):
            if (ii > opts.save[0]) and (ii < (n_eigs - opts.save[1])): continue
            aux = problem.state_to_output( mtx_phi[:,ii] )
            for name, val in aux.iteritems():
                out[name+'%03d' % ii] = val

        if post_process_hook is not None:
            out = post_process_hook( out, problem, mtx_phi )

        problem.domain.mesh.write( ofn_trunk + '.vtk', io = 'auto', out = out )

        fd = open( ofn_trunk + '_eigs.txt', 'w' )
        eigs.tofile( fd, ' ' )
        fd.close()

        return Struct( eigs = eigs, eig_vectors = eig_vectors )


    def eval_coef_e( self ):
        dconf_raw = self.conf.options.dispersion_conf
        dconf = ProblemConf.from_dict( *dconf_raw )

        dconf.materials = self.conf.materials
        dconf.fe = self.conf.fe
        dconf.regions.update( self.conf.regions )

        dproblem = ProblemDefinition.from_conf( dconf, init_variables = False )

        req = dconf.requirements['pis']
        pis = create_pis( dproblem, req['variables'][0] )

        req = dconf.requirements['corrs_phono_rs']
        solve_corrs = CorrectorsRS( dproblem, req )

        fc = self.conf.options.file_conf
        save_hook = make_save_hook( dproblem.ofn_trunk + fc['corrs_rs'],
                                    self.post_process_hook )
        corrs_rs = solve_corrs( pis, save_hook = save_hook )

        volume = eval_term_op( None, 'd_volume.i1.Y( uy )', self.problem )

        cargs = dconf.coefs['E']
        get_coef = CoefE( dproblem, cargs )
        mtx = get_coef( pis, corrs_rs, volume )

        return mtx

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'detect_band_gaps' :
    'detect frequency band gaps',
    'analyze_dispersion' :
    'analyze dispersion properties (low frequency domain)',
    'plot' :
    'plot frequency band gaps, assumes -b',
}

def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "-b", "--band-gaps",
                       action = "store_true", dest = "detect_band_gaps",
                       default = False, help = help['detect_band_gaps'] )
    parser.add_option( "-d", "--dispersion",
                       action = "store_true", dest = "analyze_dispersion",
                       default = False, help = help['analyze_dispersion'] )
    parser.add_option( "-p", "--plot",
                       action = "store_true", dest = "plot",
                       default = False, help = help['plot'] )

    options, args = parser.parse_args()
    if options.plot:
        if pylab is None:
            output( 'pylab cannot be imported, ignoring option -p!' )
            options.plot = False
        elif options.analyze_dispersion == False:
            options.detect_band_gaps = True

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    required, other = get_standard_keywords()
    required.remove( 'solver_[0-9]+|solvers' )
    conf = ProblemConf.from_file( filename_in, required, other )

    app = AcousticBandGapsApp( conf, options, 'eigen:' )
    opts = conf.options
    if hasattr( opts, 'parametric_hook' ): # Parametric study.
        parametric_hook = getattr( conf, opts.parametric_hook )
        app.parametrize( parametric_hook )
    app()

if __name__ == '__main__':
##     mtx_k = io.read_sparse_matrix_hdf5( '1todo/K.h5', output_format = 'csr' )
##     print mtx_k.__repr__()
##     mtx_m = io.read_sparse_matrix_hdf5( '1todo/M.h5', output_format = 'csr' )
##     print mtx_m.__repr__()
##     mtx_k.save( 'k.txt', format='%d %d %.12f\n' )
##     mtx_m.save( 'm.txt', format='%d %d %.12f\n' )
##     eigs, mtx_s_phi = eig( mtx_k.toarray(), mtx_m.toarray(),
##                          print_time = True )
##     print eigs
##     eigs, aux = eig( mtx_m.toarray(),
##                      print_time = True )
##     print eigs
##     pause()
    main()
