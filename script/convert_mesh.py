#!/usr/bin/env python
"""
Convert a mesh file from one SfePy-supported format to another.

Examples::

  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s2.5
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.vtk -s0.5,2,1 -c 0
  $ ./script/convert_mesh.py meshes/3d/cylinder.mesh new.mesh --remesh='q2/0 a1e-8'
"""
from __future__ import absolute_import
import sys
import os.path as op
from six.moves import range
sys.path.append('.')

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from sfepy.base.base import nm, output
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import (output_mesh_formats, MeshIO,
                                       supported_cell_types)
from sfepy.discrete.fem.mesh import fix_double_nodes
from sfepy.mesh.mesh_tools import elems_q2t

helps = {
    'scale' : 'scale factor (float or comma-separated list for each axis)'
    ' [default: %(default)s]',
    'center' : 'center of the output mesh (0 for origin or'
    ' comma-separated list for each axis) applied after scaling'
    ' [default: %(default)s]',
    'refine' : 'uniform refinement level [default: %(default)s]',
    'format' : 'output mesh format (overrides filename_out extension)',
    'list' : 'list supported readable/writable output mesh formats',
    'merge' : 'remove duplicate vertices',
    'tri-tetra' : 'convert elements: quad->tri, hexa->tetra',
    '2d' : 'force a 2D mesh by removing the z coordinates - assumes a 3D mesh'
    ' in the xy plane',
    'save-per-mat': 'extract cells by material id and save them into'
    ' separate mesh files with a name based on filename_out and the id'
    ' numbers (preserves original mesh vertices)',
    'remesh' : 'when given, extract the mesh surface and remesh its volume'
    ' using tetgen. The option value define the required mesh quality, for'
    ' example: --remesh="q2/0 a0.1" - consult tetgen documentation.',
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
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--scale', metavar='scale',
                        action='store', dest='scale',
                        default=None, help=helps['scale'])
    parser.add_argument('-c', '--center', metavar='center',
                        action='store', dest='center',
                        default=None, help=helps['center'])
    parser.add_argument('-r', '--refine', metavar='level',
                        action='store', type=int, dest='refine',
                        default=0, help=helps['refine'])
    parser.add_argument('-f', '--format', metavar='format',
                        action='store', type=str, dest='format',
                        default=None, help=helps['format'])
    parser.add_argument('-l', '--list', action='store_true',
                        dest='list', help=helps['list'])
    parser.add_argument('-m', '--merge', action='store_true',
                        dest='merge', help=helps['merge'])
    parser.add_argument('-t', '--tri-tetra', action='store_true',
                        dest='tri_tetra', help=helps['tri-tetra'])
    parser.add_argument('-2', '--2d', action='store_true',
                        dest='force_2d', help=helps['2d'])
    parser.add_argument('--save-per-mat', action='store_true',
                        dest='save_per_mat', help=helps['save-per-mat'])
    parser.add_argument('--remesh', metavar=r'q<float>/<float>[a<float>]',
                        action='store', dest='remesh',
                        default=None, help=helps['remesh'])
    parser.add_argument('filename_in')
    parser.add_argument('filename_out')
    options = parser.parse_args()

    if options.list:
        output('Supported readable mesh formats:')
        output('--------------------------------')
        output_mesh_formats('r')
        output('')
        output('Supported writable mesh formats:')
        output('--------------------------------')
        output_mesh_formats('w')
        sys.exit(0)

    scale = _parse_val_or_vec(options.scale, 'scale', parser)
    center = _parse_val_or_vec(options.center, 'center', parser)

    filename_in = options.filename_in
    filename_out = options.filename_out

    mesh = Mesh.from_file(filename_in)

    if options.force_2d:
        data = list(mesh._get_io_data())
        data[0] = data[0][:, :2]
        mesh = Mesh.from_data(mesh.name, *data)

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

    if options.tri_tetra > 0:
        conns = None
        for k, new_desc in [('3_8', '3_4'), ('2_4', '2_3')]:
            if k in mesh.descs:
                conns = mesh.get_conn(k)
                break

        if conns is not None:
            nelo = conns.shape[0]
            output('initial mesh: %d elements' % nelo)

            new_conns = elems_q2t(conns)
            nn = new_conns.shape[0] // nelo
            new_cgroups = nm.repeat(mesh.cmesh.cell_groups, nn)

            output('new mesh: %d elements' % new_conns.shape[0])
            mesh = Mesh.from_data(mesh.name, mesh.coors,
                                  mesh.cmesh.vertex_groups,
                                  [new_conns], [new_cgroups], [new_desc])

    if options.merge:
        desc = mesh.descs[0]
        coor, ngroups, conns = fix_double_nodes(mesh.coors,
                                                mesh.cmesh.vertex_groups,
                                                mesh.get_conn(desc), 1e-9)
        mesh = Mesh.from_data(mesh.name + '_merged',
                              coor, ngroups,
                              [conns], [mesh.cmesh.cell_groups], [desc])

    if options.save_per_mat:
        desc = mesh.descs[0]
        conns, cgroups = mesh.get_conn(desc), mesh.cmesh.cell_groups
        coors, ngroups = mesh.coors, mesh.cmesh.vertex_groups
        mat_ids = nm.unique(cgroups)

        for mat_id in mat_ids:
            idxs = nm.where(cgroups == mat_id)[0]
            imesh = Mesh.from_data(mesh.name + '_matid_%d' % mat_id,
                                   coors, ngroups,
                                   [conns[idxs]], [cgroups[idxs]], [desc])

            fbase, fext = op.splitext(filename_out)
            ifilename_out = '%s_matid_%d%s' % (fbase, mat_id, fext)
            io = MeshIO.for_format(ifilename_out, format=options.format,
                                   writable=True)
            output('writing %s...' % ifilename_out)
            imesh.write(ifilename_out, io=io)
            output('...done')

    if options.remesh:
        import tempfile
        import shlex
        import subprocess
        from sfepy.base.ioutils import remove_files

        domain = FEDomain(mesh.name, mesh)
        region = domain.create_region('surf', 'vertices of surface', 'facet')
        surf_mesh = Mesh.from_region(region, mesh,
                                     localize=True, is_surface=True)

        dirname = tempfile.mkdtemp()

        surf_filename = op.join(dirname, 'surf.mesh')
        surf_mesh.write(surf_filename, io='auto')

        qopts = ''.join(options.remesh.split()) # Remove spaces.
        command = 'tetgen -BFENkACpO%s %s' % (qopts, surf_filename)
        args = shlex.split(command)
        subprocess.call(args)

        mesh = Mesh.from_file(op.join(dirname, 'surf.1.vtk'))

        remove_files(dirname)

    io = MeshIO.for_format(filename_out, format=options.format,
                           writable=True)

    cell_types = ', '.join(supported_cell_types[io.format])
    output('writing [%s] %s...' % (cell_types, filename_out))
    mesh.write(filename_out, io=io)
    output('...done')

if __name__ == '__main__':
    main()
