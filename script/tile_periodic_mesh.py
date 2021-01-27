#!/usr/bin/env python
"""
The program scales a periodic input mesh (a rectangle or box) in
filename_in by a scale factor and generates a new mesh by repeating the
scaled original mesh in a regular grid (scale x scale [x scale]) if
repeat option is None, or in a grid nx x ny x nz for repeat 'nx,ny,nz',
producing again a periodic rectangle or box mesh.
"""
from __future__ import absolute_import
import sys
sys.path.append('.')

from sfepy.base.base import Output
from sfepy.discrete.fem.mesh import Mesh
from sfepy.mesh.mesh_generators import gen_tiled_mesh
from argparse import ArgumentParser, Action

helps = {
    'scale' : 'scale factor [default: %(default)s]',
    'repeat' : 'repetition counts in each axial direction'
               ' [default: %(default)s]',
    'eps'   : 'coordinate precision [default: %(default)s]',
}

class ParseRepeat(Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, [int(r) for r in values.split(',')])

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="%(prog)s 42")
    parser.add_argument("-s", "--scale", type=int, metavar='scale',
                        action="store", dest="scale",
                        default=2, help=helps['scale'])
    parser.add_argument("-r", "--repeat", type=str, metavar='nx,ny[,nz]',
                        action=ParseRepeat, dest="repeat",
                        default=None, help=helps['repeat'])
    parser.add_argument("-e", "--eps", type=float, metavar='eps',
                        action="store", dest="eps",
                        default=1e-8, help=helps['eps'])
    parser.add_argument('filename_in')
    parser.add_argument('filename_out')
    options = parser.parse_args()

    filename_in = options.filename_in
    filename_out = options.filename_out

    output = Output('tpm:')
    output('scale:', options.scale)
    output('repeat:', options.repeat)
    output('eps:', options.eps)

    mesh_in = Mesh.from_file(filename_in)
    mesh_out = gen_tiled_mesh(mesh_in, options.repeat, 1./options.scale,
                              options.eps)
    mesh_out.write(filename_out, io='auto')
    output('done.')

if __name__ == '__main__':
    main()
