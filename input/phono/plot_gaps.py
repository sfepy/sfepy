#!/usr/bin/env python
import sys
sys.path.append( '.' )

import numpy as nm
import os.path as op
from optparse import OptionParser
import pylab

from sfepy.homogenization.phono import plot_gaps, plot_eigs

usage = """%prog [options] filename [filename, ...]"""
help = {
    'fig_name' :
    'output figure name [default: %default]',
    'same_dir' :
    'put figure into the directory of input files',
}

plot_rsc = { # Resources for all plots.
    'resonance' : {'linewidth' : 0.5, 'color' : 'k', 'linestyle' : '-' },
    'masked' : {'linewidth' : 0.2, 'color' : 'k', 'linestyle' : ':' },
    'x_axis' : {'linewidth' : 1, 'color' : 'k', 'linestyle' : '-' },
    'eig_min' : {'linewidth' : 1, 'color' : 'k', 'linestyle' : '--' },
    'eig_max' : {'linewidth' : 1, 'color' : 'k', 'linestyle' : '-' },
    'strong_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 0.5) },
    'weak_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 1) },
    'propagation' : {'linewidth' : 0, 'facecolor' : (0.5, 1, 0.5) },
##     'strong_gap' : {'linewidth' : 0, 'facecolor' : (0.6, 0.6, 0.6) },
##     'weak_gap' : {'linewidth' : 0, 'facecolor' : (0.8, 0.8, 0.8) },
##     'propagation' : {'linewidth' : 0, 'facecolor' : (1, 1, 1) },
    'params' : {'axes.labelsize': 'x-large',
                'text.fontsize': 'large',
                'legend.fontsize': 'large',
                'xtick.labelsize': 'large',
                'ytick.labelsize': 'large',
                'text.usetex': False},
}

def main():
    """
    Log file format:
      par_name: par
      squared: <bool>
      header
      header
      f0 f1 flag_min f_min v_min flag_max f_max v_max kind
      desc
    """
    parser = OptionParser( usage = usage, version = "%prog" )
    parser.add_option( "-o", "", metavar = 'fig_name',
                       action = "store", dest = "fig_name",
                       default = 'band_gaps.png', help = help['fig_name'] )
    parser.add_option( "-d", "--same-dir", metavar = 'same_dir',
                       action = "store_true", dest = "same_dir",
                       default = False, help = help['same_dir'] )
    options, args = parser.parse_args()

    n_strips = len( args )

    pylab.rcParams.update( plot_rsc['params'] )

    fig = pylab.figure( 1 )
    fig.clf()

    pars = []
    for ii, arg in enumerate( args ):
        print arg
        
        fd = open( arg, 'r' )
        par_name, par = fd.readline().split( ':' )
        pars.append( '%.2f' % float( par ) )
        squared = eval( fd.readline().split( ':' )[1] )
        n_zeroed = eval( fd.readline().split( ':' )[1] )
        n_eigs = eval( fd.readline().split( ':' )[1] )
        print 'zeroed eigenmomenta: %d of %d' % (n_zeroed, n_eigs)

        # Skip gaps header.
        fd.readline()
        fd.readline()
        n_interval = int( fd.readline() )
        
        freq_range = []
        gaps = []
        kinds = []
        for ir in xrange( n_interval ):
            vv = fd.readline().split()
            desc = fd.readline()

            freq_range.append( float( vv[0] ) )
            last_f = float( vv[1] )
            gmin = (int( vv[2] ), float( vv[3] ), float( vv[4] ))
            gmax = (int( vv[5] ), float( vv[6] ), float( vv[7] ))
            gaps.append( (gmin, gmax) )
            kinds.append( (vv[8], desc.strip()) )
        freq_range.append( last_f )
        freq_range = nm.array( freq_range )

        # Skip resonances header.
        fd.readline()
        n_resonance = int( fd.readline() )
        valid = []
        resonances = []
        n_valid = 0
        for ir in xrange( n_resonance ):
            vv = fd.readline().split()
            valid.append( int( vv[0] ) )
            if valid[-1]:
                n_valid += 1
            resonances.append( float( vv[1] ) )
        fd.close()

##         print freq_range
##         print gaps
##         print kinds

        plot_gaps( 1, plot_rsc, gaps, kinds, freq_range, [ii-0.4, ii+0.4] )
        plot_eigs( 1, plot_rsc, valid, resonances, [ii-0.4, ii+0.4] )

    ax = fig.gca()
    ax.set_xlim( [freq_range[0], freq_range[-1]] )
    ax.set_yticks( range( n_strips ) )
    ax.set_yticklabels( pars )
    ax.set_ylabel( par_name )
    if squared:
        ax.set_xlabel( r'$\lambda$, $\omega^2$' )
    else:
        ax.set_xlabel( r'$\sqrt{\lambda}$, $\omega$' )
    pylab.axis( 'tight' )

    print par_name

    fig_name = options.fig_name
    if options.same_dir:
        in_dir = op.split( args[0] )[0]
        fig_name = op.join( in_dir, op.split( fig_name )[1] )

    fig.savefig( fig_name )
    pylab.show()

if __name__ == '__main__':
    main()
