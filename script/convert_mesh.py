#!/usr/bin/env python
import sys
sys.path.append('.')

from optparse import OptionParser
from sfepy.base.base import nm, output
from sfepy.fem import Mesh

usage = """%prog [options] filename_in filename_out

Convert a mesh file from one SfePy-supported format to another.

Examples:

$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk
$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s2.5
$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1
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
            try:
                scale = float(scale)
            except ValueError:
                scale = [float(ii) for ii in scale.split(',')]
            scale = nm.array(scale, dtype=nm.float64, ndmin=1)
        except:
            output('bad scale! (%s)' % scale)
            parser.print_help()
            sys.exit(1)
        
    filename_in, filename_out = args
    
    mesh = Mesh.from_file(filename_in)

    if scale is not None:
        if len(scale) == 1:
            tr = nm.eye(mesh.dim, dtype=nm.float64) * scale
        elif len(scale) == mesh.dim:
            tr = nm.diag(scale)
        else:
            raise ValueError('bad scale! (%s)' % scale)
        mesh.transform_coors(tr)

    output('writing %s...' % filename_out)
    mesh.write(filename_out, io='auto')
    output('...done')

if __name__ == '__main__':
    main()
