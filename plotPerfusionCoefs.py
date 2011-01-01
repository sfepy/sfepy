#!/usr/bin/env python
# 13.07.2007, c
import os
import os.path as op
from optparse import OptionParser
try:
    import pylab
except:
    pylab = None

import numpy as nm

from sfepy.base.base import import_file, get_default_attr
from sfepy.homogenization.coefficients import Coefficients

##
# c: 13.07.2007, r: 03.04.2008
def show_array( arr, name ):
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot( 111 )
    ax.set_title( name )
    ax.xaxis.tick_top()
    ax.title.set_y( 1.05 ) # raise it up a bit for tick top
    if arr.ndim == 1:
        arr = arr[:,nm.newaxis]
    pylab.imshow( arr, aspect = 'auto',
                  interpolation = 'nearest', origin = 'lower' )
    ax.set_ylim( ax.get_ylim()[::-1] )
    pylab.colorbar()
    fig.canvas.draw()
    pylab.ioff()
    return fig

##
# c: 13.07.2007, r: 03.04.2008
def plot_history( arr, name ):
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot( 111 )
    ax.set_title( name )
    ax.hold( True )
    legend = []
    if arr.ndim == 2:
        arr = arr[...,nm.newaxis]
    for ir in range( arr.shape[1] ):
        for ic in range( arr.shape[2] ):
            ax.plot( arr[0:,ir,ic] )
            legend.append( '%d%d' % (ir, ic) )
    ax.legend( legend )
    fig.canvas.draw()
    pylab.ioff()
    return fig

##
# 16.07.2007, c
def plot_volume_fractions( vfs, name ):
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot( 111 )
    ax.set_title( name )
    labels = []
    vals = []
    for key, val in vfs.iteritems():
        labels.append( r'$|%s|$' % key[6:] )
        vals.append( val )
    ax.pie( vals, labels = labels, autopct = '%.2f%%' )
    ax.axis( 'image' )
    fig.canvas.draw()
    pylab.ioff()
    return fig

##
# c: 20.06.2008, r: 23.06.2008
def load_coefs( filename ):
    coefs = Coefficients.from_file_hdf5( filename )

    try:
        options = import_file( coefs.filename ).options
    except:
        options = None

    if options is not None:
        plot_info = options['plot_info']
        tex_names = options['tex_names']

    else:
        plot_info = get_default_attr(coefs, 'plot_info', {})
        tex_names = get_default_attr(coefs, 'tex_names', {})

    return coefs, plot_info, tex_names

##
# c: 28.08.2007, r: 23.06.2008
def plot_and_save( filename_in, dir_name, fname = None, keys = None,
                 interactive = True ):
    """If keys == None take all."""
    if fname is None:
        def get_filename( name, format ):
            return op.join( dir_name, name + '.%s' % format )
    else:
        def get_filename( name, format ):
            d_name = op.join( dir_name, name )
            try:
                os.mkdir( d_name )
            except:
                pass
            return op.join( d_name, '%s.%s' % (fname, format) )

    coefs, plot_info, tex_names = load_coefs( filename_in )
    all_keys = coefs.to_dict().keys()
    print 'all:', all_keys
    if keys is None:
        keys = all_keys
    else:
        keys = keys.split()
    print 'selected:', keys

    format = 'png'

    vfs = {}
    for key, val in coefs.to_dict().iteritems():
        if key.startswith( 'volume' ):
            if key == 'volume_y': continue
            if key not in keys: continue
            vfs[key] = val
        print key, nm.asarray( val ).shape

    if vfs:
        fig = plot_volume_fractions( vfs, r'volume fractions' )
        fig.savefig( get_filename( 'coef_vf', format ) )

    pylab.rcParams['lines.linewidth'] = 3.0
    pylab.rcParams['legend.fontsize'] = 20.0
    pylab.rcParams['axes.labelsize'] = 20.0
    pylab.rcParams['xtick.labelsize'] = 20.0
    pylab.rcParams['ytick.labelsize'] = 20.0

    for key, val in coefs.to_dict().iteritems():
        if not plot_info.has_key( key ): continue
        if key not in keys: continue

        val = nm.asarray( val )

        kind = plot_info[key][0]
        name = tex_names[key]
        fig_name = 'coef' + key

        print kind, name, fig_name
        
        for indx in plot_info[key][1]:
            if indx == '...':
                suffix = ''
            else:
                suffix = indx
            if kind == 'mtx':
                fig = show_array( eval( 'val[%s]' % indx ), name )
            elif kind == 'th':
                fig = plot_history( eval( 'val[%s]' % indx ), name )
            fig.savefig( get_filename( fig_name + suffix, format ) )
    if interactive:
        pylab.show()
    else:
        return coefs, plot_info, tex_names

usage = """%prog [options] filename_in"""

help = {
    'dir_name' :
    'output directory [default: %default]',
    'keys' :
    'plot only selected keys [default: %default]',
}

##
# c: 13.07.2007, r: 19.06.2008
def main():
    install_dir = op.dirname( __file__ )
    version = open( op.join( install_dir, 'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-d", "", metavar = 'dir_name',
                       action = "store", dest = "dir_name",
                       default = '.', help = help['dir_name'] )
    parser.add_option( "-k", "--keys", metavar = 'keys',
                       action = "store", dest = "keys",
                       default = None, help = help['keys'] )
    (options, args) = parser.parse_args()

    if (len( args ) >= 1):
        filenames_in = args[0:];
    else:
        parser.print_help(),
        return

    if len( filenames_in ) == 1:
        plot_and_save( filenames_in[0], options.dir_name, keys = options.keys )
    else:
        for iseq, filename_in in enumerate( filenames_in ):
            try: 
                fname = op.splitext( op.basename( filename_in ) )[0][6:]
            except:
                fname = '%d' % iseq
            plot_and_save( filename_in, options.dir_name,
                         fname = fname, keys = options.keys )

if __name__ == '__main__':
    main()
#########
# Example: ./plot_perfusion_coefs.py --keys="BiotH H BiotMHR" coefs\:input_bones_bones_micro.h5
