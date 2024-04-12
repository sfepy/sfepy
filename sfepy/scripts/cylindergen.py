#!/usr/bin/env python
"""
Cylinder mesh generator.
"""
import sys
sys.path.append('.')
from argparse import ArgumentParser
import os.path as op

import numpy as nm

from sfepy.base.base import output
from sfepy.mesh.mesh_generators import gen_cylinder_mesh
from sfepy.discrete.fem.meshio import check_format_suffix, MeshIO

helps = {
    'filename' :
    'output file name [default: %(default)s]',
    'format' : 'output mesh format (overrides output file name extension)',
    'axis' :
    'axis of the cylinder, one of x, y, z [default: %(default)s]',
    'dims' :
    'dimensions of the cylinder: inner surface semi-axes a1, b1, outer'
    ' surface semi-axes a2, b2, length. If the length is 0, a planar annular'
    ' mesh topology is generated embedded in 3D space or, with --2d option,'
    ' in 2D space, forcing all z coordinates to be zero [default: %(default)s]',
    'shape' :
    'shape (counts of nodes in radial, circumferential and longitudinal'\
    ' directions) of the cylinder mesh [default: %(default)s]',
    'centre' :
    'centre of the cylinder [default: %(default)s]',
    'force_hollow' :
    'force hollow mesh even if inner radii a1 = b1 = 0',
    'is_open' :
    'generate an open cylinder segment',
    'open_angle' :
    'opening angle in radians [default: %(default)s]',
    'non_uniform' :
    'space the mesh nodes in radial direction so that the element'\
    ' volumes are (approximately) the same, making thus the elements towards'\
    ' the outer surface thinner',
    '2d' :
    'set axis to z, length and all z coordinates to zero and generate'
    ' an annular mesh in 2D space',
}

def gen_cylinder(options):
    dims = nm.array(eval(options.dims), dtype=nm.float64)
    shape = nm.array(eval(options.shape), dtype=nm.int32)
    centre = nm.array(eval(options.centre), dtype=nm.float64)

    output.prefix = 'cylindergen:'
    output('dimensions:', dims)
    output('shape:', shape)
    output('centre:', centre)
    output('output file:', options.output_filename)

    check_format_suffix(options.format,
                        op.splitext(options.output_filename)[1][1:])

    mesh = gen_cylinder_mesh(dims, shape, centre,
                             axis=options.axis,
                             force_hollow=options.force_hollow,
                             is_open=options.is_open,
                             open_angle=options.open_angle,
                             non_uniform=options.non_uniform,
                             make_2d=options.is_2d,
                             name=options.output_filename)

    io = MeshIO.any_from_filename(options.output_filename,
                                  file_format=options.format, mode='w')

    mesh.write(options.output_filename, io=io)

def add_args(parser):
    parser.add_argument('-o', metavar = 'filename',
                        action = "store", dest = "output_filename",
                        default = 'out.vtk', help = helps['filename'])
    parser.add_argument('-f', '--format', metavar='format',
                        action='store', type=str, dest='format',
                        default=None, help=helps['format'])
    parser.add_argument("-a", "--axis", metavar = 'axis',
                        action = "store", dest = "axis",
                        default = 'x', help = helps['axis'])
    parser.add_argument("-d", "--dims", metavar = 'dims',
                        action = "store", dest = "dims",
                        default = '[1.0, 1.0, 2.0, 2.0, 3.0]',
                        help = helps['dims'])
    parser.add_argument("-s", "--shape", metavar = 'shape',
                        action = "store", dest = "shape",
                        default = '[11, 11, 11]', help = helps['shape'])
    parser.add_argument("-c", "--centre", metavar = 'centre',
                        action = "store", dest = "centre",
                        default = '[0.0, 0.0, 0.0]', help = helps['centre'])
    parser.add_argument("--force-hollow",
                        action = "store_true", dest = "force_hollow",
                        default = False, help = helps['force_hollow'])
    parser.add_argument("--is-open",
                        action = "store_true", dest = "is_open",
                        default = False, help = helps['is_open'])
    parser.add_argument("--open-angle", metavar = 'angle', type=float,
                        action = "store", dest = "open_angle",
                        default = '0.0', help = helps['open_angle'])
    parser.add_argument("--non-uniform",
                        action = "store_true", dest = "non_uniform",
                        default = False, help = helps['non_uniform'])
    parser.add_argument("-2", "--2d",
                        action = "store_true", dest = "is_2d",
                        default = False, help = helps['2d'])

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version = "%(prog)s")
    add_args(parser)

    options = parser.parse_args()
    gen_cylinder(options)

if __name__ == '__main__':
    main()
