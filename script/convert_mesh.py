#!/usr/bin/env python
import sys
sys.path.append('.')

from optparse import OptionParser
import numpy as nm 
from sfepy.fem import Mesh

usage = """%prog [options] filename_in filename_out

Convert a mesh file from one SfePy-supported format to another.
"""

help = {
    'scale' : 'scale factor [default: %default]',
}

def main():
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--scale", type=int, metavar='scale',
                      action="store", dest="scale",
                      default=None, help=help['scale'])
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        sys.exit(1)

    filename_in, filename_out = args
    
    mesh = Mesh.from_file(filename_in)

    if options.scale is not None:
        tr = nm.eye(mesh.dim, dtype=nm.float64) * options.scale
        mesh.transform_coords(tr)

    mesh.write(filename_out, io='auto')

if __name__ == '__main__':
    main()
