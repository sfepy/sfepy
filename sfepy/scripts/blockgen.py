#!/usr/bin/env python
"""
Block mesh generator.
"""
import sys
sys.path.append('.')
from argparse import ArgumentParser
import os.path as op

import numpy as nm

from sfepy.base.base import output
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import check_format_suffix, MeshIO

helps = {
    'filename' :
    'output file name [default: %(default)s]',
    'format' : 'output mesh format (overrides output file name extension)',
    'dims' :
    'dimensions  of the block [default: %(default)s]',
    'shape' :
    'shape (counts of nodes in x, y, z) of the block [default: %(default)s]',
    'centre' :
    'centre of the block [default: %(default)s]',
    '2d' :
    'generate a 2D rectangular mesh, the third components of the above'
    ' options are ignored',
}

def gen_block(options):
    dim = 2 if options.is_2d else 3

    dims = nm.array(eval(options.dims), dtype=nm.float64)[:dim]
    shape = nm.array(eval(options.shape), dtype=nm.int32)[:dim]
    centre = nm.array(eval(options.centre), dtype=nm.float64)[:dim]

    output.prefix = 'blockgen:'
    output('dimensions:', dims)
    output('shape:', shape)
    output('centre:', centre)
    output('output file:', options.output_filename)

    check_format_suffix(options.format,
                        op.splitext(options.output_filename)[1][1:])

    mesh = gen_block_mesh(dims, shape, centre, name=options.output_filename)

    io = MeshIO.any_from_filename(options.output_filename,
                                  file_format=options.format, mode='w')

    mesh.write(options.output_filename, io=io)

def add_args(parser):
    parser.add_argument('-o', metavar='filename',
                        action='store', dest='output_filename',
                        default='out.vtk', help=helps['filename'])
    parser.add_argument('-f', '--format', metavar='format',
                        action='store', type=str, dest='format',
                        default=None, help=helps['format'])
    parser.add_argument('-d', '--dims', metavar='dims',
                        action='store', dest='dims',
                        default='[1.0, 1.0, 1.0]', help=helps['dims'])
    parser.add_argument('-s', '--shape', metavar='shape',
                        action='store', dest='shape',
                        default='[11, 11, 11]', help=helps['shape'])
    parser.add_argument('-c', '--centre', metavar='centre',
                        action='store', dest='centre',
                        default='[0.0, 0.0, 0.0]', help=helps['centre'])
    parser.add_argument('-2', '--2d',
                        action='store_true', dest='is_2d',
                        default=False, help=helps['2d'])

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    add_args(parser)

    options = parser.parse_args()
    gen_block(options)

if __name__ == '__main__':
    main()
