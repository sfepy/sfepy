#!/usr/bin/env python
import numpy as nm

from sfepy.base.base import Output
from sfepy.fem.mesh import Mesh
from sfepy.fem.mesh_generators import compose_periodic_mesh
from optparse import OptionParser

usage = """%prog [options] filename_in filename_out

The program scales a periodic input mesh (a rectangle or box) in
filename_in by a scale factor and generates a new mesh by repeating the
scaled original mesh in a regular grid (scale x scale [x scale]) if
repeat option is None, or in a grid nx x ny x nz for repeat 'nx,ny,nz',
producing again a periodic rectangle or box mesh.
"""

help = {
    'scale' : 'scale factor [default: %default]',
    'repeat' : 'repetition counts in each axial direction'
               ' [default: %default]',
    'eps'   : 'coordinate precision [default: %default]',
    'nomvd' : 'omit mesh periodicity test using minimum vertex distance'\
    + ' (it is demanding in cpu time and memory) [default: %default]',
}

def parse_repeat(option, opt, value, parser):
    if value is not None:
        setattr(parser.values, option.dest, [int(r) for r in value.split(',')])

def main():
    parser = OptionParser(usage=usage, version="%prog 42")
    parser.add_option("-s", "--scale", type=int, metavar='scale',
                      action="store", dest="scale",
                      default=2, help=help['scale'])
    parser.add_option("-r", "--repeat", type='str', metavar='nx,ny[,nz]',
                      action="callback", dest="repeat",
                      callback=parse_repeat, default=None, help=help['repeat'])
    parser.add_option("-e", "--eps", type=float, metavar='eps',
                      action="store", dest="eps",
                      default=1e-8, help=help['eps'])
    parser.add_option("-n", "--no-mvd",
                      action="store_true", dest="nomvd",
                      default=False, help=help['nomvd'])
    (options, args) = parser.parse_args()

    if (len( args ) == 2):
        filename_in = args[0]
        filename_out = args[1]
    else:
        parser.print_help()
        return

    output = Output('genPerMesh:')
    output('scale:', options.scale)
    output('repeat:', options.repeat)
    output('eps:', options.eps)

    mesh_in = Mesh.from_file(filename_in)
    mesh_out = compose_periodic_mesh(mesh_in, options.scale, options.repeat,
                                     options.eps, check_mvd=not options.nomvd)
    mesh_out.write(filename_out, io='auto')
    output('done.')
    
if __name__ == '__main__':
    main()
