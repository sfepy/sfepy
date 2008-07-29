#!/usr/bin/env python
from sfepy.base.base import *
from sfepy.fem.mesh import Mesh, find_map, merge_mesh, make_mesh
from sfepy.base.la import cycle
from optparse import OptionParser

##
# 24.05.2007, c
def test():
    x1 = nm.fix( nm.random.rand( 5, 3 ) * 5 ) / 10
    x2 = nm.fix( nm.random.rand( 10, 3 ) * 5 ) / 10
    cmap = find_map( x1, x2 )
    conn1 = nm.array( [[range( 5 )], [range( 5 )]] ).squeeze()
    conn2 = nm.array( [[range( 10 )], [range( 10 )]] ).squeeze()
    conn2.shape = (4, 5)

    print x1
    print conn1
    print x2
    print conn2
    print cmap
    xx, conns = merge_mesh( x1, [conn1], x2, [conn2], cmap )
    print xx
    print conns

##
# c: 05.05.2008, r: 05.05.2008
def fix_double_nodes( coor, conns, eps ):
    n_nod, dim = coor.shape
    cmap = find_map( coor, nm.zeros( (0,dim) ), eps = eps, allow_double = True )
    if cmap.size:
        print 'double nodes in input mesh!'
        print 'trying to fix...'

        while cmap.size:
            print cmap.size

            # Just like in Variable.equation_mapping()...
            ii = nm.argsort( cmap[:,1] )
            scmap = cmap[ii]

            eq = nm.arange( n_nod )
            eq[scmap[:,1]] = -1
            eqi = eq[eq >= 0]
            eq[eqi] = nm.arange( eqi.shape[0] )
            remap = eq.copy()
            remap[scmap[:,1]] = eq[scmap[:,0]]
            print coor.shape
            coor = coor[eqi]
            print coor.shape
            ccs = []
            for conn in conns:
                ccs.append( remap[conn] )
            conns = ccs
            cmap = find_map( coor, nm.zeros( (0,dim) ), eps = eps,
                            allow_double = True )
        print '...done'
    return coor, conns

##
# c: 25.05.2007, r: 05.05.2008
# 28.05.2007
def get_min_edge_size( coor, conns ):

    mes = 1e16
    for conn in conns:
        n_ep = conn.shape[1]
        for ir in range( n_ep ):
            x1 = coor[conn[:,ir]]
            for ic in range( ir + 1, n_ep ):
                x2 = coor[conn[:,ic]]
                aux = nm.sqrt( nm.sum( (x2 - x1)**2.0, axis = 1 ).min() )
                mes = min( mes, aux )

    return mes

##
# 25.05.2007, c
def get_min_vertex_distance( coor, guess ):
    """Can miss the minimum, but is enough for our purposes."""
    # Sort by x.
    ix = nm.argsort( coor[:,0] )
    scoor = coor[ix]

    mvd = 1e16
    
    # Get mvd in chunks potentially smaller than guess.
    n_coor = coor.shape[0]
    print n_coor
    
    i0 = i1 = 0
    x0 = scoor[i0,0]
    while 1:
        while ((scoor[i1,0] - x0) < guess) and (i1 < (n_coor - 1)):
            i1 += 1

#        print i0, i1, x0, scoor[i1,0]
        aim, aa1, aa2, aux = get_min_vertex_distance_naive( scoor[i0:i1+1] )
        if aux < mvd:
            im, a1, a2 = aim, aa1 + i0, aa2 + i0
        mvd = min( mvd, aux )
        i0 = i1 = int( 0.5 * (i1 + i0 ) ) + 1
#        i0 += 1
        x0 = scoor[i0,0]
#        print '-', i0

        if i1 == n_coor - 1: break

    print im, ix[a1], ix[a2], a1, a2, scoor[a1], scoor[a2]

    return mvd
        
##
# c: 25.05.2007, r: 05.05.2008
def get_min_vertex_distance_naive( coor ):

    ii = nm.arange( coor.shape[0] )
    i1, i2 = nm.meshgrid( ii, ii )
    i1 = i1.flatten()
    i2 = i2.flatten()

    ii = nm.where( i1 < i2 )
    aux = coor[i1[ii]] - coor[i2[ii]]
    aux = nm.sum( aux**2.0, axis = 1 )

    im = aux.argmin()

    return im, i1[ii][im], i2[ii][im], nm.sqrt( aux[im] )

usage = """%prog [options] filename_in filename_out

The program scales a periodic input mesh (a rectangle or box) in filename_in
by a scale factor and generates a new mesh by repeating the scaled original
mesh in a regular grid (scale x scale [x scale]), producing again a periodic
rectagle or box mesh."""

help = {
    'scale' : 'scale factor [default: %default]',
    'eps'   : 'coordinate precision [default: %default]',
    'test'  : 'test the code',
    'nomvd' : 'omit mesh periodicity test using minimum vertex distance'\
    + ' (it is demanding in cpu time and memory) [default: %default]',
}

##
# c: 23.05.2007, r: 06.05.2008
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
        filename_in = args[0]
        filename_out = args[1]
    else:
        parser.print_help()
        return

    print 'scale:', options.scale
    print 'eps:', options.eps

    mesh_in = Mesh.from_file( filename_in )
    bbox = mesh_in.get_bounding_box()
    print 'bbox:\n', bbox
    mscale = bbox[1] - bbox[0]
    centre0 = 0.5 * (bbox[1] + bbox[0])
    print 'centre:\n', centre0

    scale = nm.array( options.scale, dtype = nm.float64 )

    # Normalize original coordinates.
    coor0 = (mesh_in.nod0[:,:-1] - centre0) / (mscale)
    dim = mesh_in.dim

    coor0, mesh_in.conns = fix_double_nodes( coor0, mesh_in.conns, options.eps )
    if not options.nomvd:
        mes0 = get_min_edge_size( coor0, mesh_in.conns )
        mvd0 = get_min_vertex_distance( coor0, mes0 )
        if mes0 > (mvd0 + options.eps):
            print '          original min. "edge" length: %.5e' % mes0
            print 'original approx. min. vertex distance: %.5e' % mvd0
            print '-> still double nodes in input mesh!'
            print 'try increasing eps...'
            raise ValueError

    for indx in cycle( [options.scale] * dim ):
        aindx = nm.array( indx, dtype = nm.float64 )
        centre = 0.5 * (2.0 * aindx - scale + 1.0)
        print indx, centre

        if aindx.sum() == 0:
            coor = coor0 + centre
            conns = mesh_in.conns
        else:
            coor1 = coor0 + centre
            conns1 = mesh_in.conns

            cmap = find_map( coor, coor1, eps = options.eps )
            if not cmap.size:
                print 'non-periodic mesh!'
#                raise ValueError
            else:
                print cmap.size / 2
            coor, conns = merge_mesh( coor, conns, coor1, conns1, cmap,
                                     eps = options.eps )

    if not options.nomvd:
        mes = get_min_edge_size( coor, conns )
        mvd = get_min_vertex_distance( coor, mes0 )

        print '          original min. "edge" length: %.5e' % mes0
        print '             final min. "edge" length: %.5e' % mes
        print 'original approx. min. vertex distance: %.5e' % mvd0
        print '   final approx. min. vertex distance: %.5e' % mvd
        if mvd < 0.99999 * mvd0:
            if mvd0 < (mes0 - options.eps):
                print '-> probably non-periodic input mesh!'
                print '   ... adjacent sides were not connected!'
                print '   try increasing eps...'
            else:
                print '-> input mesh might be periodic'
                print '   try increasing eps...'
        else:
            print '-> input mesh looks periodic'
    else:
        print 'non-periodic input mesh detection skipped!'

    print 'renormalizing...'
    coor = (coor * mscale) / scale
    print 'saving...'
    mesh_out = make_mesh( coor, conns, mesh_in )
    mesh_out.write( filename_out, io = 'auto' )
    print 'done.'
    
if __name__ == '__main__':
    main()

##     import sfepy.fem.extmods.meshutils as mu
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
##     mu.sort_rows( x, nm.array( [0,1,2], nm.int32 ) )
##     print time.clock() - tt
##     print x

##     assert nm.all( x == x2 )
