#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.la import eig
from sfepy.fem.evaluate import eval_term_op
import sfepy.base.ioutils as io
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.homogenization.phono import process_options, get_method,\
     transform_plot_data, plot_logs, plot_gaps, detect_band_gaps

##
# c: 25.09.2007, r: 08.04.2008
def solve_eigen_problem( conf, options ):

    if options.output_filename_trunk:
        ofn_trunk = options.output_filename_trunk
    else:
        ofn_trunk = io.get_trunk( conf.filename_mesh )

    pb = ProblemDefinition.from_conf( conf )
    dim = pb.domain.mesh.dim

    pb.time_update()

    dummy = pb.create_state_vector()
    mtx_a = eval_term_op( dummy, conf.equations['lhs'], pb,
                       dw_mode = 'matrix', tangent_matrix = pb.mtx_a )
    mtx_b = eval_term_op( dummy, conf.equations['rhs'], pb,
                       dw_mode = 'matrix', tangent_matrix = pb.mtx_a.copy() )

##     mtx_a.save( 'a.txt', format='%d %d %.12f\n' )
##     mtx_b.save( 'b.txt', format='%d %d %.12f\n' )
    output( 'computing resonance frequencies...' )
    tt = [0]
    eigs, mtx_s_phi = eig( mtx_a.toarray(), mtx_b.toarray(), return_time = tt,
                         method = get_method( conf.options ) )
    output( '...done in %.2f s' % tt[0] )
    output( eigs )
##     import sfepy.base.plotutils as plu
##     plu.spy( mtx_b, eps = 1e-12 )
##     plu.pylab.show()
##     pause()
    n_eigs = eigs.shape[0]
    opts = process_options( conf.options, n_eigs )

    mtx_phi = nm.empty( (pb.variables.di.ptr[-1], mtx_s_phi.shape[1]),
                       dtype = nm.float64 )
    for ii in xrange( n_eigs ):
        mtx_phi[:,ii] = pb.variables.make_full_vec( mtx_s_phi[:,ii] )

    out = {}
    for ii in xrange( n_eigs ):
        if (ii > opts.save[0]) and (ii < (n_eigs - opts.save[1])): continue
        aux = pb.state_to_output( mtx_phi[:,ii] )
        key = aux.keys()[0]
        out[key+'%03d' % ii] = aux[key]

    pb.domain.mesh.write( ofn_trunk + '.vtk', io = 'auto', out = out )

    fd = open( ofn_trunk + '_eigs.txt', 'w' )
    eigs.tofile( fd, ' ' )
    fd.close()

    return Struct( pb = pb, eigs = eigs, mtx_phi = mtx_phi )


usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'detect_band_gaps' :
    'detect frequency band gaps [default: %default]',
    'plot' :
    'plot frequency band gaps [default: %default], assumes -b',
}

##
# c: 25.09.2007, r: 08.04.2008
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
    parser.add_option( "-p", "--plot",
                       action = "store_true", dest = "plot",
                       default = False, help = help['plot'] )

    options, args = parser.parse_args()
    if options.plot:
        options.detect_band_gaps = True

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    set_output_prefix( 'eigen:' )

    required, other = get_standard_keywords()
    required.remove( 'solver_[0-9]+|solvers' )
    conf = ProblemConf.from_file( filename_in, required, other )
##     print conf
##     pause()

    evp = solve_eigen_problem( conf, options )

    if options.detect_band_gaps:
        bg = detect_band_gaps( evp.pb, evp.eigs, evp.mtx_phi, conf, options )

        if options.plot:
            plot_range, tlogs = transform_plot_data( bg.logs,
                                                  bg.opts.plot_tranform,
                                                  conf.funmod )
            plot_gaps( 1, bg.gaps, bg.kinds, bg.freq_range_margins,
                      plot_range, show = False )
            plot_logs( 1, tlogs, bg.freq_range,
                      plot_range, bg.opts.squared, show = True )

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
