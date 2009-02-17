#!/usr/bin/env python
import sys
sys.path.append('.')

from optparse import OptionParser
import numpy as nm 
from sfepy.fem import Mesh

usage = """%prog [options] filename_in filename_out

Convert a mesh file from one SfePy-supported format to another.

Examples:

$script/convert_mesh.py database/simple.mesh new.vtk
$script/convert_mesh.py database/simple.mesh new.vtk -s2.5
$script/convert_mesh.py database/simple.mesh new.vtk -s0.5,2,1
"""

help = {
    'scale' : 'scale factor [default: %default]',
}

def main():
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--scale", metavar='scale',
                      action="store", dest="scale",
                      default=None, help=help['scale'])
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        sys.exit(1)

    scale = options.scale
    if scale is not None:
        try:
            scale = float(scale)
        except ValueError:
            scale = [float(ii) for ii in scale.split(',')]
        scale = nm.array(scale, dtype=nm.float64, ndmin=1)
        
    filename_in, filename_out = args
    
    mesh = Mesh.from_file(filename_in)

    if scale is not None:
        if len(scale) == 1:
            tr = nm.eye(mesh.dim, dtype=nm.float64) * scale
        elif len(scale) == mesh.dim:
            tr = nm.diag(scale)
        mesh.transform_coords(tr)

    mesh.write(filename_out, io='auto')

if __name__ == '__main__':
    main()
