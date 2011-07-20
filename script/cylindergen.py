#!/usr/bin/env python
import sys
sys.path.append( '.' )
from optparse import OptionParser
from sfepy.fem import gen_cylinder_mesh

usage = """%prog [options]

Cylinder mesh generator.
"""
help = {
    'filename' :
    'output file name [default: %default]',
    'axis' :
    'axis of the cylinder, one of x, y, z [default: %default]',
    'dims' :
    'dimensions of the cylinder: inner surface semi-axes a1, b1, outer'\
    ' surface semi-axes a2, b2, length [default: %default]',
    'shape' :
    'shape (counts of nodes in radial, circumferential and longitudinal'\
    ' directions) of the cylinder mesh [default: %default]',
    'centre' :
    'centre of the cylinder [default: %default]',
    'force_hollow' :
    'force hollow mesh even if inner radii a1 = b1 = 0',
    'is_open' :
    'generate an open cylinder segment',
    'open_angle' :
    'opening angle in radians [default: %default]',
    'non_uniform' :
    'space the mesh nodes in radial direction so that the element'\
    ' volumes are (approximately) the same, making thus the elements towards'\
    ' the outer surface thinner',
}

def main():
    parser = OptionParser( usage = usage, version = "%prog" )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename",
                       default = 'out.vtk', help = help['filename'] )
    parser.add_option( "-a", "--axis", metavar = 'axis',
                       action = "store", dest = "axis",
                       default = 'x', help = help['axis'] )
    parser.add_option( "-d", "--dims", metavar = 'dims',
                       action = "store", dest = "dims",
                       default = '[1.0, 1.0, 2.0, 2.0, 3.0]',
                       help = help['dims'] )
    parser.add_option( "-s", "--shape", metavar = 'shape',
                       action = "store", dest = "shape",
                       default = '[11, 11, 11]', help = help['shape'] )
    parser.add_option( "-c", "--centre", metavar = 'centre',
                       action = "store", dest = "centre",
                       default = '[0.0, 0.0, 0.0]', help = help['centre'] )
    parser.add_option( "", "--force-hollow",
                       action = "store_true", dest = "force_hollow",
                       default = False, help = help['force_hollow'] )
    parser.add_option( "", "--is-open",
                       action = "store_true", dest = "is_open",
                       default = False, help = help['is_open'] )
    parser.add_option( "", "--open-angle", metavar = 'angle', type='float',
                       action = "store", dest = "open_angle",
                       default = '0.0', help = help['open_angle'] )
    parser.add_option( "", "--non-uniform",
                       action = "store_true", dest = "non_uniform",
                       default = False, help = help['non_uniform'] )
    (options, args) = parser.parse_args()

    import numpy as nm
    dims = eval( "nm.array( %s, dtype = nm.float64 )" % options.dims )
    shape = eval( "nm.array( %s, dtype = nm.int32 )" % options.shape )
    centre = eval( "nm.array( %s, dtype = nm.float64 )" % options.centre )

    print dims
    print shape
    print centre

    mesh = gen_cylinder_mesh(dims, shape, centre,
                             axis=options.axis,
                             force_hollow=options.force_hollow,
                             is_open=options.is_open,
                             open_angle=options.open_angle,
                             non_uniform=options.non_uniform,
                             name=options.output_filename)
    mesh.write( options.output_filename, io = 'auto' )

if __name__ == '__main__':
    main()
