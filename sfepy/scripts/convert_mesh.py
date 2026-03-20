#!/usr/bin/env python
"""
Convert a mesh file from one SfePy-supported format to another.
"""
import sys
import os.path as op
from ast import literal_eval
sys.path.append('.')

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from sfepy.base.base import nm, output
from sfepy.base.ioutils import remove_files
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import output_mesh_formats
from sfepy.discrete.fem.mesh import fix_double_nodes
from sfepy.discrete.fem.utils import prepare_translate
from sfepy.linalg import make_axis_rotation_matrix
import sfepy.mesh.mesh_tools as mt
from sfepy.mesh.mesh_generators import gen_tiled_mesh

helps = {
    'scale' : 'scale factor (float or comma-separated list for each axis)'
    ' [default: %(default)s]',
    'center' : 'center of the output mesh (0 for origin or'
    ' comma-separated list for each axis) applied after scaling'
    ' [default: %(default)s]',
    'rot_axis' : """rotation axis (x, y, z or comma-separated list of
      coordinates) [default: %(default)s]""",
    'rot_angle' : """rotation angle in degrees around rotation axis
      [default: %(default)s]""",
    'refine' : 'uniform refinement level [default: %(default)s]',
    'format' : 'output mesh format (overrides filename_out extension)',
    'list' : 'list supported readable/writable output mesh formats',
    'merge' : 'remove duplicate vertices',
    'tri-tetra' : 'convert elements: quad->tri, hexa->tetra',
    '2d' : 'force a 2D mesh by removing the z coordinates - assumes a 3D mesh'
    ' in the xy plane',
    '3d' : 'force a 3D mesh by adding zero (y,) z coordinates',
    'cell-dim' : 'write only cells of a given dimension,'
    ' use a comma-separated list for several values',
    'cell-vertices-only': 'remove vertices not used in any cell',
    'save-per-mat': 'extract cells by material id and save them into'
    ' separate mesh files with a name based on filename_out and the id'
    ' numbers (preserves original mesh vertices)',
    'remap-vertex-groups' : 'remap vertex groups to a contiguous range from 0',
    'remap-cell-groups' : 'remap cell groups (materials) to a contiguous'
    ' range from 0',
    'tile' : """
      scale a periodic input mesh (a rectangle or box) by a scale factor
      (--scale with a single integer value) and generate a new mesh by
      repeating the scaled original mesh in a regular grid (scale x scale [x
      scale]) if --tile=scale, or in a grid nx x ny[ x nz] for
      --tile=<nx,ny[,nz]>, producing again a periodic rectangle or box mesh.
      Use the --eps option to set the coordinate precision for periodicity
      checks.
    """,
    'extract_edges' : """
      extract outline edges of a given mesh and save it into a meshio-supported
      file. The outline edge is an edge for which norm(nvec1 - nvec2) < eps,
      where nvec1 and nvec2 are the normal vectors of the incident facets and
      eps is given by the --eps option.
    """,
    'eps' : """
      tolerance parameter [default: %(default)s]
    """,
    'extract_surface' : 'extract surface facets of a mesh',
    'print_surface' : """
       extract surface facets of a mesh and print it to a given file (use '-'
       for stdout) in form of a list where each row is [element, face,
       component]. A component corresponds to a contiguous surface region - for
       example, a cubical mesh with a spherical hole has two surface
       components. Two surface faces sharing a single node belong to one
       component.
    """,
    'remesh' : """when given, remesh the given mesh using tetgen.
      The options can be the following, separated by spaces, in this order: 1.
      "r" causes remeshing of the mesh volume - if not present the mesh surface
      is extracted and used for the volume mesh generation. 2.
      "q[<float>/<float>]" (required) - the two numbers after "q" are a maximum
      radius-edge ratio bound and a minimum dihedral angle bound. 3. "a<float>"
      (optional) - the number imposes a maximum volume constraint on all
      tetrahedra. 4. O[<0-9>/<0-7>] - the two numbers correspond to a mesh
      optimization level and a choice of optimizing operations. 5. "V"
      (optional) - if present, mesh statistics are printed. Consult the tetgen
      documentation for details. Examples: --remesh='rq2/0 a1e-8 O9/7 V'""",
    'extract_by_matid': """extract mesh cells by a given "mat_id".
      The "math_id" can be an integer, range '(a,b)', or a list '[a,b,c]'.""",
    'extract_by_nodegroup': """extract mesh cells where all nodes have
      the "node_groups" values equal to a given "node_group". The
      "node_group" can be an integer, range '(a,b)', or a list '[a,b,c]'.""",
    'extrude': """extrude the given 2D mesh. The options can be the following:
      "c[<list of [<float>,<float>,<float>]>]" - defines the points of the
      extrusion path, each point means one element layer, or
      "c[<float>,<int>]" - sets the extrusion path in the z-direction with
      defining the element thickness and the number of element layers, the
      extrusion path, each point means one element layer, or
      "t" (optional) - defines the twist angle in each layer,
      "s" (optional) - sets scaling in each layer.
      Example: --extrude='c0.1,10 t0.2 s0.9'""",
    'revolve': """revolve the given 2D mesh around the axis. The point on
      the axis is defined by options "p[<float>,<float>,<float>]" and the axis
      direction by "v[<float>,<float>,<float>]". The revolution angle can be
      specified in degrees by "m[<float>]", it is equal to 360 deg, if not
      given. The number of divisions is given by "n[<int>]".
      Example: --revolve='p[-1,-2,0] v[1,0,0] n12 m180'""",
    'mirror': """mirror the given mesh using a plane defined by a point and
      a normal vector. Example: --mirror='p[-0.5,-0.2,0] v[0,1,0]'""",
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


def _parse_fun_args(s, args_tab):
    args, dargs = {}, None

    for arg in s.split():
        for a, n in args_tab:
            if arg.startswith(a):
                args[n] = nm.array(literal_eval(arg[1:]))
                break
        else:
            dargs = nm.array(literal_eval(arg))

    return args, dargs


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--scale', metavar='scale',
                        action='store', dest='scale',
                        default=None, help=helps['scale'])
    parser.add_argument('-c', '--center', metavar='center',
                        action='store', dest='center',
                        default=None, help=helps['center'])
    parser.add_argument('--rot-axis', metavar='axis',
                        action='store', dest='rot_axis',
                        default=None, help=helps['rot_axis'])
    parser.add_argument('--rot-angle', metavar='angle',
                        action='store', type=float, dest='rot_angle',
                        default=0.0, help=helps['rot_angle'])
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
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-2', '--2d', action='store_true',
                       dest='force_2d', help=helps['2d'])
    group.add_argument('-3', '--3d', action='store_true',
                       dest='force_3d', help=helps['3d'])
    parser.add_argument('-d', '--cell-dim', metavar='cell_dim',
                        action='store', dest='cell_dim',
                        default=None, help=helps['cell-dim'])
    parser.add_argument('--cell-vertices-only', action='store_true',
                        dest='cell_vertices_only',
                        help=helps['cell-vertices-only'])
    parser.add_argument('--save-per-mat', action='store_true',
                        dest='save_per_mat', help=helps['save-per-mat'])
    parser.add_argument('--remap-vertex-groups', action='store_true',
                        dest='remap_vertex_groups',
                        help=helps['remap-vertex-groups'])
    parser.add_argument('--remap-cell-groups', action='store_true',
                        dest='remap_cell_groups',
                        help=helps['remap-cell-groups'])
    parser.add_argument('--tile', metavar='nx,ny[,nz]', dest='tile',
                        default=None, help=helps['tile'])
    parser.add_argument('--extract-edges', action='store_true',
                        dest='extract_edges', help=helps['extract_edges'])
    parser.add_argument('--eps', action='store', type=float, dest='eps',
                        default=1e-12, help=helps['eps'])
    parser.add_argument('--extract-surface', action='store_true',
                        dest='extract_surface', help=helps['extract_surface'])
    parser.add_argument('--print-surface', action='store', metavar='filename',
                        dest='print_surface', default=None,
                        help=helps['print_surface'])
    parser.add_argument('--remesh', metavar='options',
                        action='store', dest='remesh',
                        default=None, help=helps['remesh'])
    parser.add_argument('-i', '--extract-by-matid', metavar='matid',
                        action='store', dest='extract_by_matid',
                        default=None, help=helps['extract_by_matid'])
    parser.add_argument('-n', '--extract-by-nodegroup', metavar='node_group',
                        action='store', dest='extract_by_nodegroup',
                        default=None, help=helps['extract_by_nodegroup'])
    parser.add_argument('--extrude', metavar='options',
                        action='store', dest='extrude',
                        default=None, help=helps['extrude'])
    parser.add_argument('--revolve', metavar='options',
                        action='store', dest='revolve',
                        default=None, help=helps['revolve'])
    parser.add_argument('--mirror', metavar='options',
                        action='store', dest='mirror',
                        default=None, help=helps['mirror'])
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

    if options.rot_axis in ('x', 'y', 'z'):
        rot_axis = dict(
            x=[1.0, 0.0, 0.0],
            y=[0.0, 1.0, 0.0],
            z=[0.0, 0.0, 1.0],
        )[options.rot_axis]

    else:
        rot_axis = _parse_val_or_vec(options.rot_axis, 'rot_axis', parser)

    cell_dim = _parse_val_or_vec(options.cell_dim, 'cell_dim', parser)

    filename_in = options.filename_in
    filename_out = options.filename_out

    if options.remesh:
        import tempfile
        import shlex
        import subprocess

        dirname = tempfile.mkdtemp()

        is_surface = options.remesh.startswith('q')
        if is_surface:
            mesh = Mesh.from_file(filename_in)
            domain = FEDomain(mesh.name, mesh)
            region = domain.create_region('surf', 'vertices of surface',
                                          'facet')
            surf_mesh = Mesh.from_region(region, mesh,
                                         localize=True, is_surface=True)

            filename = op.join(dirname, 'surf.mesh')
            surf_mesh.write(filename, io='auto')

        else:
            import shutil

            shutil.copy(filename_in, dirname)
            filename = op.join(dirname, op.basename(filename_in))

        qopts = ''.join(options.remesh.split()) # Remove spaces.
        command = 'tetgen -BFENkACp%s %s' % (qopts, filename)
        args = shlex.split(command)
        subprocess.call(args)

        root, ext = op.splitext(filename)
        mesh = Mesh.from_file(root + '.1.vtk')

        remove_files(dirname)

    else:
        mesh = Mesh.from_file(filename_in)

    if cell_dim is not None:
        cell_dim = [int(ii) for ii in cell_dim]
        data = list(mesh._get_io_data(cell_dim_only=cell_dim))
        mesh = Mesh.from_data(mesh.name, *data)

    if options.force_2d:
        data = list(mesh._get_io_data(cell_dim_only=2))
        data[0] = data[0][:, :2]
        mesh = Mesh.from_data(mesh.name, *data)

    elif options.force_3d:
        data = list(mesh._get_io_data())
        data[0] = nm.pad(data[0], [(0, 0), (0, 3 - data[0].shape[1])])
        mesh = Mesh.from_data(mesh.name, *data)

    if options.remap_vertex_groups:
        vgs = nm.unique(mesh.cmesh.vertex_groups)
        remap = prepare_translate(vgs, nm.arange(len(vgs)))
        mesh.cmesh.vertex_groups[:] = remap[mesh.cmesh.vertex_groups]

    if options.remap_cell_groups:
        cgs = nm.unique(mesh.cmesh.cell_groups)
        remap = prepare_translate(cgs, nm.arange(len(cgs)))
        mesh.cmesh.cell_groups[:] = remap[mesh.cmesh.cell_groups]

    if scale is not None:
        if len(scale) == 1:
            tr = nm.eye(mesh.dim, dtype=nm.float64) * scale
        elif len(scale) == mesh.dim:
            tr = nm.diag(scale)
        else:
            raise ValueError('bad scale! (%s)' % scale)
        mesh.transform_coors(tr)

    if center is not None:
        mesh.set_centre(center)

    if rot_axis is not None:
        tr = make_axis_rotation_matrix(rot_axis,
                                       nm.pi * options.rot_angle / 180.0)
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
        mesh = mt.triangulate(mesh, verbose=True)

    if options.merge:
        desc = mesh.descs[0]
        coor, ngroups, conns = fix_double_nodes(mesh.coors,
                                                mesh.cmesh.vertex_groups,
                                                mesh.get_conn(desc))
        mesh = Mesh.from_data(mesh.name + '_merged',
                              coor, ngroups,
                              [conns], [mesh.cmesh.cell_groups], [desc])

    if options.cell_vertices_only:
        mesh = mt.get_cell_vertices_only(mesh)

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
            output('writing %s...' % ifilename_out)
            imesh.write(ifilename_out, file_format=options.format)
            output('...done')

    if options.tile:
        _scale = (scale[0] if isinstance(scale, nm.ndarray)
                  else scale if scale is not None else 1)
        repeat = (list(literal_eval(options.tile)) if options.tile != 'scale'
                  else None)
        output('tiling paramaters:')
        output('scale:', scale)
        output('repeat:', repeat)
        output('eps:', options.eps)

        mesh = gen_tiled_mesh(mesh, repeat, 1./_scale, options.eps)

    if options.extract_surface or options.print_surface:
        domain = FEDomain(mesh.name, mesh)

    if options.extract_surface:
        region = domain.create_region('surf', 'vertices of surface', 'facet')
        surf_mesh = Mesh.from_region(region, mesh,
                                     localize=True, is_surface=True)
        mesh = surf_mesh

    if options.extract_by_matid is not None:
        cgroup = options.extract_by_matid
        if isinstance(cgroup, str):
            cgroup = literal_eval(cgroup)

        mesh = mt.get_mesh_by(mesh, cgroup)

    if options.extract_by_nodegroup is not None:
        ngroup = options.extract_by_nodegroup
        if isinstance(ngroup, str):
            ngroup = literal_eval(ngroup)

        mesh = mt.get_mesh_by_ngroup(mesh, ngroup)

    if options.extrude is not None:
        args, darg = _parse_fun_args(options.extrude,
                                     [('c', 'cline'), ('t', 'twist'),
                                      ('s', 'scale'), ('v', 'nvec')])
        if darg is not None and 'cline' not in args:
            args['cline'] = darg

        mesh = mt.extrude(mesh, **args)

    if options.revolve is not None:
        args, _ = _parse_fun_args(options.revolve,
                                  [('p', 'p'), ('v', 'v'),
                                   ('n', 'nphi'), ('m', 'phi_max')])
        mesh = mt.revolve(mesh, **args)

    if options.mirror is not None:
        args, _ = _parse_fun_args(options.mirror, [('p', 'p'), ('v', 'v')])
        mesh = mt.mirror(mesh, **args)

    if options.print_surface:
        domain.fix_element_orientation()

        lst, surf_faces = mt.get_surface_faces(domain)

        gr_s = mt.surface_graph(surf_faces, domain.mesh.n_nod)

        n_comp, comps = mt.surface_components(gr_s, surf_faces)
        output('number of surface components:', n_comp)

        ccs, comps = comps, nm.zeros((0,1), nm.int32)
        for cc in ccs:
            comps = nm.concatenate((comps, cc[:,nm.newaxis]), 0)

        out = nm.concatenate((lst, comps), 1)
        if (options.print_surface == '-'):
            file_out = sys.stdout
        else:
            file_out = open(options.print_surface, 'w')
        for row in out:
            file_out.write('%d %d %d\n' % (row[0], row[1], row[2]))
        if (options.print_surface != '-'):
            file_out.close()

    if options.extract_edges:
        mesh_out = mt.extract_edges(mesh, eps=options.eps)
        mesh_out = mt.merge_lines(mesh_out)

        import meshio
        emesh = meshio.Mesh(mesh_out[0], [('line', mesh_out[2][0])],
                            cell_data={'mat_id' : mesh_out[3]})
        emesh.write(filename_out)

    else:
        output('writing %s...' % filename_out)
        mesh.write(filename_out, file_format=options.format, binary=False)
        output('...done')

if __name__ == '__main__':
    main()
