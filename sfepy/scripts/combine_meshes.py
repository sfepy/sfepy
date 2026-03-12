#!/usr/bin/env python
"""
Create a mesh by combining one or more meshes centered at specified coordinates.

Usage Examples
--------------

- Create a cylinder mesh::

    sfepy-mesh cylinder -a z -d 0.3,0.3,0.34,0.34,0.45 -s 2,8,8 -o cyl1.vtk
    sfepy-mesh info cyl1.vtk

- Place the mesh at two sites::

    sfepy-mesh combine 'cyl1.vtk=[0.0,0.0,-0.275]' 'cyl1.vtk=[0.0,0.0,+0.275]' -o cut-cyl1.vtk
    sfepy-mesh info cut-cyl1.vtk
"""
import sys
sys.path.append('.')
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from sfepy.base.conf import dict_from_string as parse_as_dict
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import MeshIO

helps = {
    'filename' :
    'output file name [default: %(default)s]',
    'format' : 'output mesh format (overrides output file name extension)',
    'meshes' :
    'mesh filenames and coordinates where their centres should be placed',
}

def combine_meshes(options):
    infos = [list(parse_as_dict(mesh_info, free_word=True).items())[0]
             for mesh_info in options.meshes]

    combined_mesh = None
    meshes = {}
    for info in infos:
        mesh = meshes.get(info[0])
        if mesh is None:
            meshes[info[0]] = mesh = Mesh.from_file(info[0])

        mesh.set_centre(info[1])

        if combined_mesh is None:
            combined_mesh = mesh.copy()

        else:
            combined_mesh = combined_mesh + mesh

    io = MeshIO.any_from_filename(options.output_filename,
                                  file_format=options.format, mode='w')
    combined_mesh.write(options.output_filename, io=io)

def add_args(parser):
    parser.add_argument('-o', metavar='filename',
                        action='store', dest='output_filename',
                        default='out.vtk', help=helps['filename'])
    parser.add_argument('-f', '--format', metavar='format',
                        action='store', type=str, dest='format',
                        default=None, help=helps['format'])
    parser.add_argument(metavar='filename=[x,y,z]', nargs='+',
                        action='store', dest='meshes',
                        help=helps['meshes'])

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    add_args(parser)

    options = parser.parse_args()
    combine_meshes(options)

if __name__ == '__main__':
    main()
