#!/usr/bin/env python
"""
Convert a mesh file from one SfePy-supported format to another.

Examples:

$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk
$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s2.5
$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1
"""
import sys
sys.path.append('.')

from optparse import OptionParser
from sfepy.base.base import nm, output
from sfepy.fem import Mesh
from sfepy.fem.meshio import output_writable_meshes, MeshIO

usage = '%prog [options] filename_in filename_out\n' + __doc__.rstrip()

help = {
    'scale' : 'scale factor [default: %default]',
    'format' : 'output mesh format (overrides filename_out extension)',
    'list' : 'list supported writable output mesh formats',
}

def main():
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--scale', metavar='scale',
                      action='store', dest='scale',
                      default=None, help=help['scale'])
    parser.add_option('-f', '--format', metavar='format',
                      action='store', type='string', dest='format',
                      default=None, help=help['format'])
    parser.add_option('-l', '--list', action='store_true',
                      dest='list', help=help['list'])
    (options, args) = parser.parse_args()

    if options.list:
        output_writable_meshes()
        sys.exit(0)

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

    io = MeshIO.for_format(filename_out, format=options.format,
                           writable=True)

    output('writing %s...' % filename_out)
    mesh.write(filename_out, io=io)
    output('...done')

if __name__ == '__main__':
    main()
