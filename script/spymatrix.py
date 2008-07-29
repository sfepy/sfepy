#!/usr/bin/env python
# 23.11.2005, c 
from pylab import *
import scipy as nm
from optparse import OptionParser

usage = """%prog [options] filename"""

# 23.11.2005, c
# 24.11.2005
# 27.11.2005
# 10.08.2006
def main():
    parser = OptionParser( usage = usage )
    parser.add_option( "-e", "--epsilon", type = float,
                       dest = "eps", default = None,
                       help = "set drop-off tolerance [default: %default]" )
    (options, args) = parser.parse_args()
    if len( args ) < 1:
        print usage
        return
    filename = args[0]
    print filename + ':'

    fd = open( filename, "r" )
    n_row, n_col = map( int, fd.readline().split() )
    n_item = int( fd.readline() )
    print n_row, n_col, n_item

    ij = nm.zeros( (n_item,2), nm.int32 )
    val = nm.zeros( (n_item,), nm.float64 )
    for ii, row in enumerate( fd.readlines() ):
        aux = row.split()
        ij[ii] = int( aux[0] ), int( aux[1] )
        val[ii] = float( aux[2] )

    if options.eps is not None:
        print 'using', options.eps
        ij = nm.compress( nm.absolute( val ) > options.eps, ij, 0 )
        n_item = ij.shape[0]
    else:
        print 'showing all'

    print n_item
    if n_item:
        plot( ij[:,1] + 0.5, ij[:,0] + 0.5, linestyle = 'None',
              marker = ',', markersize = 0.5, markeredgewidth = 0.1 )
    axis( [-0.5, n_row+0.5, -0.5, n_col+0.5] )
    axis( 'image' )
    xlabel( '%d x %d: %d nnz, %.2f\%% fill'
            % (n_row, n_col, n_item, 100. * n_item / float( n_row * n_col )) )
    gca().set_ylim( gca().get_ylim()[::-1] )
    show()
    
if (__name__ == '__main__'):
    main()
