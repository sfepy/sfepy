#!/usr/bin/env python
import sys
sys.path.append('.')

from optparse import OptionParser
from sfepy.base.base import nm, output
from sfepy.fem import Mesh
from sfepy.fem.meshio import io_table, supported_capabilities

usage = """%prog [options] filename_in filename_out

Convert a mesh file from one SfePy-supported format to another.

Examples:

$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk
$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s2.5
$ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1
"""

help = {
    'scale' : 'scale factor [default: %default]',
    'format' : 'output mesh format (overrides filename_out extension)',
    'list' : 'list supported writable output mesh formats',
}

def output_writable_meshes():
    output('Supported writable mesh formats are:')
    for key, val in supported_capabilities.iteritems():
        if 'w' in val:
            output(key)

def main():
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--scale", metavar='scale',
                      action="store", dest="scale",
                      default=None, help=help['scale'])
    parser.add_option("-f", "--format", metavar='format',
                      action="store", type='string', dest="format",
                      default=None, help=help['format'])
    parser.add_option("-l", "--list", action="store_true", 
                      dest="list", help=help['list'])
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

    io = 'auto'
    if options.format:
        try:
            io = io_table[options.format](filename_out)
        except KeyError:
            output('unknown output mesh format! (%s)' % options.format)
            output_writable_meshes()
            sys.exit(1)
            
        if 'w' not in supported_capabilities[options.format]:
            output('write support not implemented for output mesh format! (%s)'
                    % options.format)
            output_writable_meshes()
            sys.exit(1)

    output('writing %s...' % filename_out)
    mesh.write(filename_out, io=io)
    output('...done')

if __name__ == '__main__':
    main()
