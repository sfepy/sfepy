#!/usr/bin/env python
"""
Generate a single IGA patch block in 2D or 3D of given degrees and continuity
using igakit.

The grid has equally-spaced knot vectors.
"""
from __future__ import absolute_import
import sys
sys.path.append('.')
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as nm

from sfepy.base.base import output
from sfepy.discrete.iga.domain_generators import gen_patch_block_domain
import sfepy.discrete.iga.plot_nurbs as pn
import sfepy.discrete.iga.io as io

helps = {
    'filename' :
    'output file name [default: block%%dd.iga]',
    'dims' :
    'dimensions of the block [default: %(default)s]',
    'centre' :
    'centre of the block [default: %(default)s]',
    'shape' :
    'numbers of unique knot values along each axis [default: %(default)s]',
    'degrees' :
    'NURBS degrees along each axis [default: %(default)s]',
    'continuity' :
    'NURBS continuity along each axis [default: degrees-1]',
    'cp_mode' :
    "The control points mode. The default 'greville' results in a uniform"
    " Bezier mesh, while the 'uniform' mode results in a uniform grid of"
    " control points a finer Bezier mesh inside the block and a coarser"
    " Bezier mesh near the block boundary.",
    '2d' :
    'generate a 2D block, the third components of the above'
    ' options are ignored',
    'plot' :
    'plot parametric, control and Bezier meshes as well as 1D basis slices',
    'label' :
    'label control and Bezier mesh points in figures',
}

def gen_iga_patch(options):
    dim = 2 if options.is_2d else 3

    filename = options.filename
    if filename is None:
        filename = 'block%dd.iga' % dim

    dims = nm.array(eval(options.dims), dtype=nm.float64)[:dim]
    centre = nm.array(eval(options.centre), dtype=nm.float64)[:dim]
    shape = nm.array(eval(options.shape), dtype=nm.int32)[:dim]
    degrees = nm.array(eval(options.degrees), dtype=nm.int32)[:dim]

    if options.continuity is None:
        continuity = degrees - 1

    else:
        continuity = nm.array(eval(options.continuity), dtype=nm.int32)[:dim]

    n_dofs = degrees + 1 + (shape - 2) * (degrees - continuity)

    output('dimensions:', dims)
    output('centre:    ', centre)
    output('shape:     ', shape)
    output('degrees:   ', degrees)
    output('continuity:', continuity)
    output('cp_mode:   ', options.cp_mode)
    output('->        :', filename)
    output('number of DOFs per axis:', n_dofs)

    nurbs, bmesh, regions = gen_patch_block_domain(dims, shape, centre,
                                                   degrees,
                                                   continuity=continuity,
                                                   cp_mode=options.cp_mode,
                                                   name='block', verbose=True)

    io.write_iga_data(filename, None, nurbs.knots, nurbs.degrees, nurbs.cps,
                      nurbs.weights, nurbs.cs, nurbs.conn,
                      bmesh.cps, bmesh.weights, bmesh.conn,
                      regions)

    if options.plot:
        pn.plt.rcParams['lines.linewidth'] = 2

        block = nurbs.nurbs
        ax = pn.plot_parametric_mesh(None, block.knots)
        ax.set_title('parametric mesh')
        ax.axis('equal')

        points = block.points[..., :dim]
        ax = pn.plot_control_mesh(None, points, label=options.label)
        ax = pn.plot_iso_lines(ax, block)
        ax.set_title('control mesh and iso lines (blue)'
                     ' in Greville abscissae coordinates')
        ax.axis('equal')

        points = bmesh.cps
        ax = pn.plot_bezier_mesh(None, points, bmesh.conn, block.degree,
                                 label=options.label)
        ax = pn.plot_iso_lines(ax, block)
        ax.set_title('Bezier mesh and iso lines (blue)'
                     ' in Greville abscissae coordinates')
        ax.axis('equal')

        pn.plt.rcParams['lines.linewidth'] = 3

        line = block.extract(0, 0.5)
        if dim == 3:
            line = line.extract(0, 0.5)
        ax = pn.plot_nurbs_basis_1d(None, line, n_points=1000,
                                    legend=options.label)
        ax.set_xlabel('last parametric coordinate')
        ax.set_title('1D NURBS basis')

        ax = pn.plot_nurbs_basis_1d(None, line, n_points=1000, x_axis=dim-1,
                                    legend=options.label)
        ax.set_xlabel('last physical coordinate')
        ax.set_title('1D NURBS basis')

        pn.plt.show()

def add_args(parser):
    parser.add_argument('-o', metavar='filename',
                        action='store', dest='filename',
                        default=None, help=helps['filename'])
    parser.add_argument('-d', '--dims', metavar='dims',
                        action='store', dest='dims',
                        default='[1.0, 1.0, 1.0]', help=helps['dims'])
    parser.add_argument('-c', '--centre', metavar='centre',
                        action='store', dest='centre',
                        default='[0.0, 0.0, 0.0]', help=helps['centre'])
    parser.add_argument('-s', '--shape', metavar='shape',
                        action='store', dest='shape',
                        default='[5, 5, 5]', help=helps['shape'])
    parser.add_argument('--degrees', metavar='degrees',
                        action='store', dest='degrees',
                        default='[2, 2, 2]', help=helps['degrees'])
    parser.add_argument('--continuity', metavar='continuity',
                        action='store', dest='continuity',
                        default=None, help=helps['continuity'])
    parser.add_argument('--cp-mode', metavar="'greville' or 'uniform'",
                        action='store', dest='cp_mode',
                        choices=['greville', 'uniform'],
                        default='greville', help=helps['cp_mode'])
    parser.add_argument('-2', '--2d',
                        action='store_true', dest='is_2d',
                        default=False, help=helps['2d'])
    parser.add_argument('-p', '--plot',
                        action='store_true', dest='plot',
                        default=False, help=helps['plot'])
    parser.add_argument('-l', '--label',
                        action='store_true', dest='label',
                        default=False, help=helps['label'])

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    add_args(parser)

    options = parser.parse_args()
    gen_iga_patch(options)

if __name__ == '__main__':
    main()
