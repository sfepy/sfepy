#!/usr/bin/env python
# 23.11.2005, c 
from pylab import *
import scipy as nm
from optparse import OptionParser

usage = """%prog [options] fileName"""

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
    fileName = args[0]
    print fileName + ':'

    fd = open( fileName, "r" )
    nRow, nCol = map( int, fd.readline().split() )
    nItem = int( fd.readline() )
    print nRow, nCol, nItem

    ij = nm.zeros( (nItem,2), nm.int32 )
    val = nm.zeros( (nItem,), nm.float64 )
    for ii, row in enumerate( fd.readlines() ):
        aux = row.split()
        ij[ii] = int( aux[0] ), int( aux[1] )
        val[ii] = float( aux[2] )

    if options.eps is not None:
        print 'using', options.eps
        ij = nm.compress( nm.absolute( val ) > options.eps, ij, 0 )
        nItem = ij.shape[0]
    else:
        print 'showing all'

    print nItem
    if nItem:
        plot( ij[:,1] + 0.5, ij[:,0] + 0.5, linestyle = 'None',
              marker = ',', markersize = 0.5, markeredgewidth = 0.1 )
    axis( [-0.5, nRow+0.5, -0.5, nCol+0.5] )
    axis( 'image' )
    xlabel( '%d x %d: %d nnz, %.2f\%% fill'
            % (nRow, nCol, nItem, 100. * nItem / float( nRow * nCol )) )
    gca().set_ylim( gca().get_ylim()[::-1] )
    show()
    
if (__name__ == '__main__'):
    main()
