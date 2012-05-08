#!/usr/bin/env python
"""
Block mesh generator.
"""
import sys
sys.path.append('.')
from optparse import OptionParser

import numpy as nm; nm

from sfepy.base.base import output
from sfepy.fem import gen_block_mesh
from sfepy.fem.meshio import MeshIO

usage = '%prog [options]\n' + __doc__.rstrip()

help = {
    'filename' :
    'output file name [default: %default]',
    'format' : 'output mesh format (overrides output file name extension)',
    'dims' :
    'dimensions  of the block [default: %default]',
    'shape' :
    'shape (counts of nodes in x, y, z) of the block [default: %default]',
    'centre' :
    'centre of the block [default: %default]',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-o', '', metavar='filename',
                      action='store', dest='output_filename',
                      default='out.vtk', help=help['filename'])
    parser.add_option('-f', '--format', metavar='format',
                      action='store', type='string', dest='format',
                      default=None, help=help['format'])
    parser.add_option('-d', '--dims', metavar='dims',
                      action='store', dest='dims',
                      default='[1.0, 1.0, 1.0]', help=help['dims'])
    parser.add_option('-s', '--shape', metavar='shape',
                      action='store', dest='shape',
                      default='[11, 11, 11]', help=help['shape'])
    parser.add_option('-c', '--centre', metavar='centre',
                      action='store', dest='centre',
                      default='[0.0, 0.0, 0.0]', help=help['centre'])
    (options, args) = parser.parse_args()

    dims = eval('nm.array(%s, dtype=nm.float64)' % options.dims)
    shape = eval('nm.array(%s, dtype=nm.int32)' % options.shape)
    centre = eval('nm.array(%s, dtype=nm.float64)' % options.centre)

    output.prefix = 'blockgen:'
    output('dimensions:', dims)
    output('shape:', shape)
    output('centre:', centre)

    mesh = gen_block_mesh(dims, shape, centre, name=options.output_filename)

    io = MeshIO.for_format(options.output_filename, format=options.format,
                           writable=True)

    mesh.write(options.output_filename, io=io)

if __name__ == '__main__':
    main()
