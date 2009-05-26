#!/usr/bin/env python
import sys
sys.path.append( '.' )
from optparse import OptionParser
from sfepy.base.base import *
from sfepy.fem import Mesh
from sfepy.base.la import cycle
from sfepy.base.progressbar import MyBar

usage = """%prog [options]

Block mesh generator.
"""
help = {
    'filename' :
    'output file name [default: %default]',
    'dims' :
    'dimension of the block [default: %default]',
    'shape' :
    'shape (counts of nodes in x, y, z) of the block [default: %default]',
    'centre' :
    'centre of the block [default: %default]',
}

##
# c: 19.06.2008, r: 19.06.2008
def main():
    parser = OptionParser( usage = usage, version = "%prog" )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename",
                       default = 'out.vtk', help = help['filename'] )
    parser.add_option( "-d", "--dims", metavar = 'dims',
                       action = "store", dest = "dims",
                       default = '[1.0, 1.0, 1.0]', help = help['dims'] )
    parser.add_option( "-s", "--shape", metavar = 'shape',
                       action = "store", dest = "shape",
                       default = '[11, 11, 11]', help = help['shape'] )
    parser.add_option( "-c", "--centre", metavar = 'centre',
                       action = "store", dest = "centre",
                       default = '[0.0, 0.0, 0.0]', help = help['centre'] )
    (options, args) = parser.parse_args()

    dims = eval( "nm.array( %s, dtype = nm.float64 )" % options.dims )
    shape = eval( "nm.array( %s, dtype = nm.int32 )" % options.shape )
    centre = eval( "nm.array( %s, dtype = nm.float64 )" % options.centre )
    
    print dims
    print shape
    print centre

    dim = shape.shape[0]

    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    grid = nm.zeros( shape, dtype = nm.int32 )
    n_nod = nm.prod( shape )
    coors = nm.zeros( (n_nod, dim), dtype = nm.float64 )

    # This is 3D only...
    bar = MyBar( "       nodes:" )
    bar.init( n_nod )
    for ii, ic in enumerate( cycle( shape ) ):
        ix, iy, iz = ic
        grid[ix,iy,iz] = ii
        coors[ii] = x0 + ic * dd
        if not (ii % 100):
            bar.update( ii )
    print
    n_el = nm.prod( shape - 1 )
    conn = nm.zeros( (n_el, 8), dtype = nm.int32 )
    bar = MyBar( "       elements:" )
    bar.init( n_el )
    for ii, (ix, iy, iz) in enumerate( cycle( shape - 1 ) ):
        conn[ii,:] = [grid[ix  ,iy  ,iz  ], grid[ix+1,iy  ,iz  ],
                      grid[ix+1,iy+1,iz  ], grid[ix  ,iy+1,iz  ],
                      grid[ix  ,iy  ,iz+1], grid[ix+1,iy  ,iz+1],
                      grid[ix+1,iy+1,iz+1], grid[ix  ,iy+1,iz+1]]
        if not (ii % 100):
            bar.update( ii )
    print
    mat_id = nm.zeros( (n_el,), dtype = nm.int32 )
    desc = '3_8'

    mesh = Mesh.from_data( options.output_filename,
                          coors, None, [conn], [mat_id], [desc] )
    mesh.write( options.output_filename, io = 'auto' )

if __name__ == '__main__':
    main()
