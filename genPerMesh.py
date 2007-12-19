#!/usr/bin/env python
from sfe.base.base import *
from sfe.fem.mesh import Mesh, findMap, mergeMesh, makeMesh
from sfe.base.la import cycle
from optparse import OptionParser

##
# 24.05.2007, c
def test():
    x1 = nm.fix( nm.random.rand( 5, 3 ) * 5 ) / 10
    x2 = nm.fix( nm.random.rand( 10, 3 ) * 5 ) / 10
    cmap = findMap( x1, x2 )
    conn1 = nm.array( [[range( 5 )], [range( 5 )]] ).squeeze()
    conn2 = nm.array( [[range( 10 )], [range( 10 )]] ).squeeze()
    conn2.shape = (4, 5)

    print x1
    print conn1
    print x2
    print conn2
    print cmap
    xx, conns = mergeMesh( x1, [conn1], x2, [conn2], cmap )
    print xx
    print conns

##
# 25.05.2007, c
# 28.05.2007
def getMinEdgeSize( coor, conns ):

    mes = 1e16
    for conn in conns:
        nEP = conn.shape[1]
        for ir in range( nEP ):
            x1 = coor[conn[:,ir]]
            for ic in range( ir + 1, nEP ):
                x2 = coor[conn[:,ic]]
                aux = nm.sum( (x2 - x1)**2.0, axis = 1 ).min()
                mes = min( mes, aux )

    return mes

##
# 25.05.2007, c
def getMinVertexDistance( coor, guess ):
    """Can miss the minimum, but is enough for our purposes."""
    # Sort by x.
    ix = nm.argsort( coor[:,0] )
    scoor = coor[ix]

    mvd = 1e16
    
    # Get mvd in chunks potentially smaller than guess.
    nCoor = coor.shape[0]
    print nCoor
    
    i0 = i1 = 0
    x0 = scoor[i0,0]
    while 1:
        while ((scoor[i1,0] - x0) < guess) and (i1 < (nCoor - 1)):
            i1 += 1

#        print i0, i1, x0, scoor[i1,0]
        aim, aa1, aa2, aux = getMinVertexDistanceNaive( scoor[i0:i1+1] )
        if aux < mvd:
            im, a1, a2 = aim, aa1 + i0, aa2 + i0
        mvd = min( mvd, aux )
        i0 = i1 = int( 0.5 * (i1 + i0 ) ) + 1
#        i0 += 1
        x0 = scoor[i0,0]
#        print '-', i0

        if i1 == nCoor - 1: break

    print im, ix[a1], ix[a2], a1, a2, scoor[a1], scoor[a2]

    return mvd
        
##
# 25.05.2007, c
# 28.05.2007
def getMinVertexDistanceNaive( coor ):

    ii = nm.arange( coor.shape[0] )
    i1, i2 = nm.meshgrid( ii, ii )
    i1 = i1.flatten()
    i2 = i2.flatten()

    ii = nm.where( i1 < i2 )
    aux = coor[i1[ii]] - coor[i2[ii]]
    aux = nm.sum( aux**2.0, axis = 1 )

    im = aux.argmin()

    return im, i1[ii][im], i2[ii][im], aux[im]

usage = """%prog [options] fileNameIn fileNameOut

The program scales a periodic input mesh (a unit square or cube) in fileNameIn
by a scale factor and generates a new mesh by repeating the scaled original
mesh in a regular grid (scale x scale [x scale]), producing again a periodic
unit square or cube mesh."""

help = {
    'scale' : 'scale factor [default: %default]',
    'eps'   : 'coordinate precision [default: %default]',
    'test'  : 'test the code',
    'nomvd' : 'omit mesh periodicity test using minimum vertex distance'\
    + ' (it is demanding in cpu time and memory) [default: %default]',
}

##
# 23.05.2007, c
# 24.05.2007
# 25.05.2007
# 28.05.2007
def main():

    parser = OptionParser( usage = usage, version = "%prog 42" )
    parser.add_option( "-s", "--scale", type = int, metavar = 'scale',
                       action = "store", dest = "scale",
                       default = 2, help = help['scale'] )
    parser.add_option( "-e", "--eps", type = float, metavar = 'eps',
                       action = "store", dest = "eps",
                       default = 1e-8, help = help['eps'] )
    parser.add_option( "-t", "--test",
                       action = "store_true", dest = "test",
                       default = False, help = help['test'] )
    parser.add_option( "-n", "--no-mvd",
                       action = "store_true", dest = "nomvd",
                       default = False, help = help['nomvd'] )
    (options, args) = parser.parse_args()

    if options.test:
        test()
        return

    if (len( args ) == 2):
        fileNameIn = args[0]
        fileNameOut = args[1]
    else:
        parser.print_help()
        return

    print 'scale:', options.scale
    print 'eps:', options.eps

    meshIn = Mesh.fromFile( fileNameIn )
    bbox = meshIn.getBoundingBox()
    print 'bbox:\n', bbox
    mscale = bbox[1] - bbox[0]
    centre0 = 0.5 * (bbox[1] + bbox[0])
    print 'centre:\n', centre0

    scale = nm.array( options.scale, dtype = nm.float64 )

    # Normalize original coordinates.
    coor0 = (meshIn.nod0[:,:-1] - centre0) / (mscale[0])
    dim = meshIn.dim

    if not options.nomvd:
        mes0 = getMinEdgeSize( coor0, meshIn.conns )
        mvd0 = getMinVertexDistance( coor0, mes0 )
    
    for indx in cycle( [options.scale] * dim ):
        aindx = nm.array( indx, dtype = nm.float64 )
        centre = 0.5 * (2.0 * aindx - scale + 1.0)
        print indx, centre

        if aindx.sum() == 0:
            coor = coor0 + centre
            conns = meshIn.conns
        else:
            coor1 = coor0 + centre
            conns1 = meshIn.conns

            cmap = findMap( coor, coor1, eps = options.eps )
            if not cmap.size:
                print 'non-periodic mesh!'
#                raise ValueError
            else:
                print cmap.size / 2
            coor, conns = mergeMesh( coor, conns, coor1, conns1, cmap,
                                     eps = options.eps )

    if not options.nomvd:
        mes = getMinEdgeSize( coor, conns )
        mvd = getMinVertexDistance( coor, mes0 )

        print '          original min. "edge" length: %.5e' % mes0
        print '             final min. "edge" length: %.5e' % mes
        print 'original approx. min. vertex distance: %.5e' % mvd0
        print '   final approx. min. vertex distance: %.5e' % mvd
        if mvd < 0.99999 * mvd0:
            print '-> probably non-periodic input mesh!'
            print '   ... adjacent sides were not connected!'
        else:
            print '-> input mesh looks periodic'
    else:
        print 'non-periodic input mesh detection skipped!'

    print 'renormalizing...'
    coor = coor / scale
    print 'saving...'
    meshOut = makeMesh( coor, conns, meshIn )
    meshOut.write( fileNameOut )
    print 'done.'
    
if __name__ == '__main__':
    main()

##     import sfe.fem.extmods.meshutils as mu
##     import time
##     n = 3000000
##     x = nm.fix( nm.random.rand( n, 3 ) * 10 ) / 10
##     x = nm.array( 10 * x, dtype = nm.int32 )
##     xt = nm.transpose( x )

##     tt = time.clock()
##     ii = nm.lexsort( keys = (xt[2], xt[1], xt[0]) )
##     print time.clock() - tt
##     tt = time.clock()
##     x2 = x[ii]
##     print time.clock() - tt
##     print x2

##     tt = time.clock()
##     mu.sortRows( x, nm.array( [0,1,2], nm.int32 ) )
##     print time.clock() - tt
##     print x

##     assert nm.all( x == x2 )
