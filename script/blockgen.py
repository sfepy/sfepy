#!/usr/bin/env python
import sys
sys.path.append( '.' )
from optparse import OptionParser
from sfepy.base.base import *
from sfepy.fem import gen_block_mesh

usage = """%prog [options]

Block mesh generator.
"""
help = {
    'filename' :
    'output file name [default: %default]',
    'dims' :
    'dimensions  of the block [default: %default]',
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

    mesh = gen_block_mesh(dims, shape, centre, name=options.output_filename)
    mesh.write( options.output_filename, io = 'auto' )

if __name__ == '__main__':
    main()
