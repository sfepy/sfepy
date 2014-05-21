#!/usr/bin/env python
"""
Convert a mesh file from one SfePy-supported format to another.

Examples::

  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s2.5
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1 -c 0
"""
import sys
sys.path.append('.')

from optparse import OptionParser
from sfepy.base.base import nm, output
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import (output_writable_meshes, MeshIO,
                              supported_cell_types)

usage = '%prog [options] filename_in filename_out\n' + __doc__.rstrip()

help = {
    'scale' : 'scale factor (float or comma-separated list for each axis)'
    ' [default: %default]',
    'center' : 'center of the output mesh (0 for origin or'
    ' comma-separated list for each axis) applied after scaling'
    ' [default: %default]',
    'refine' : 'uniform refinement level [default: %default]',
    'format' : 'output mesh format (overrides filename_out extension)',
    'list' : 'list supported writable output mesh formats',
}

def _parse_val_or_vec(option, name, parser):
    if option is not None:
        try:
            try:
                option = float(option)
            except ValueError:
                option = [float(ii) for ii in option.split(',')]
            option = nm.array(option, dtype=nm.float64, ndmin=1)
        except:
            output('bad %s! (%s)' % (name, option))
            parser.print_help()
            sys.exit(1)

    return option

def main():
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--scale', metavar='scale',
                      action='store', dest='scale',
                      default=None, help=help['scale'])
    parser.add_option('-c', '--center', metavar='center',
                      action='store', dest='center',
                      default=None, help=help['center'])
    parser.add_option('-r', '--refine', metavar='level',
                      action='store', type=int, dest='refine',
                      default=0, help=help['refine'])
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

    scale = _parse_val_or_vec(options.scale, 'scale', parser)
    center = _parse_val_or_vec(options.center, 'center', parser)

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

    if center is not None:
        cc = 0.5 * mesh.get_bounding_box().sum(0)
        shift = center - cc
        tr = nm.c_[nm.eye(mesh.dim, dtype=nm.float64), shift[:, None]]
        mesh.transform_coors(tr)

    if options.refine > 0:
        domain = FEDomain(mesh.name, mesh)
        output('initial mesh: %d nodes %d elements'
               % (domain.shape.n_nod, domain.shape.n_el))

        for ii in range(options.refine):
            output('refine %d...' % ii)
            domain = domain.refine()
            output('... %d nodes %d elements'
                   % (domain.shape.n_nod, domain.shape.n_el))

        mesh = domain.mesh

    io = MeshIO.for_format(filename_out, format=options.format,
                           writable=True)

    cell_types = ', '.join(supported_cell_types[io.format])
    output('writing [%s] %s...' % (cell_types, filename_out))
    mesh.write(filename_out, io=io)
    output('...done')

if __name__ == '__main__':
    main()
