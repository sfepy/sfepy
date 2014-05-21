#!/usr/bin/env python
"""
Generate a single IGA patch block in 2D or 3D of given degrees and continuity
using igakit.

The grid has equally-spaced knot vectors. The generated control points form a
regular grid as well - this prevents coarser resolution inside the block.
"""
from optparse import OptionParser
import numpy as nm

import igakit.cad as cad

from sfepy.base.base import output
import sfepy.discrete.iga as iga
import sfepy.discrete.iga.plot_nurbs as pn
import sfepy.discrete.iga.io as io

usage = '%prog [options]\n' + __doc__.rstrip()

helps = {
    'filename' :
    'output file name [default: block%dd.iga]',
    'dims' :
    'dimensions of the block [default: %default]',
    'centre' :
    'centre of the block [default: %default]',
    'shape' :
    'numbers of unique knot values along each axis [default: %default]',
    'degrees' :
    'NURBS degrees along each axis [default: %default]',
    'continuity' :
    'NURBS continuity along each axis [default: degrees-1]',
    '2d' :
    'generate a 2D block, the third components of the above'
    ' options are ignored',
    'plot' :
    'plot parametric, control and Bezier meshes as well as 1D basis slices',
    'label' :
    'label control and Bezier mesh points in figures',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-o', '', metavar='filename',
                      action='store', dest='filename',
                      default=None, help=helps['filename'])
    parser.add_option('-d', '--dims', metavar='dims',
                      action='store', dest='dims',
                      default='[1.0, 1.0, 1.0]', help=helps['dims'])
    parser.add_option('-c', '--centre', metavar='centre',
                      action='store', dest='centre',
                      default='[0.0, 0.0, 0.0]', help=helps['centre'])
    parser.add_option('-s', '--shape', metavar='shape',
                      action='store', dest='shape',
                      default='[5, 5, 5]', help=helps['shape'])
    parser.add_option('', '--degrees', metavar='degrees',
                      action='store', dest='degrees',
                      default='[2, 2, 2]', help=helps['degrees'])
    parser.add_option('', '--continuity', metavar='continuity',
                      action='store', dest='continuity',
                      default=None, help=helps['continuity'])
    parser.add_option('-2', '--2d',
                      action='store_true', dest='is_2d',
                      default=False, help=helps['2d'])
    parser.add_option('-p', '--plot',
                      action='store_true', dest='plot',
                      default=False, help=helps['plot'])
    parser.add_option('-l', '--label',
                      action='store_true', dest='label',
                      default=False, help=helps['label'])
    (options, args) = parser.parse_args()

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

    output('dimensions:', dims)
    output('centre:    ', centre)
    output('shape:     ', shape)
    output('degrees:   ', degrees)
    output('continuity:', continuity)
    output('->        :', filename)

    dd = centre - 0.5 * dims
    block = cad.grid(shape - 1, degree=degrees, continuity=continuity)

    for ia in xrange(dim):
        block.scale(dims[ia], ia)

    for ia in xrange(dim):
        block.translate(dd[ia], ia)

    # Force uniform control points. This prevents coarser resolution inside the
    # block.
    shape = nm.asarray(block.points.shape[:-1])
    n_nod = nm.prod(shape)
    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    ngrid = nm.mgrid[[slice(ii) for ii in shape]]
    ngrid.shape = (dim, n_nod)

    coors = x0 + ngrid.T * dd
    coors.shape = shape.tolist() + [dim]

    block.array[..., :dim] = coors

    # Compute Bezier extraction data.
    cs = iga.compute_bezier_extraction(block.knots, block.degree)
    n_els = [len(ii) for ii in cs]
    conn, bconn = iga.create_connectivity(n_els, block.knots, block.degree)

    ccs = iga.combine_bezier_extraction(cs)

    cps = block.points[..., :dim].copy()
    cps = cps.reshape((-1, dim))
    bcps, bweights = iga.compute_bezier_control(cps, block.weights.ravel(), ccs,
                                                conn, bconn)

    regions = iga.get_patch_box_regions(n_els, block.degree)

    io.write_iga_data(filename, block.knots, block.degree, cps,
                      block.weights.ravel(), cs, conn, bcps, bweights, bconn,
                      regions)

    if options.plot:
        pn.plt.rcParams['lines.linewidth'] = 2

        ax = pn.plot_parametric_mesh(None, block.knots)
        ax.set_title('parametric mesh')
        ax.axis('equal')

        points = block.points[..., :dim]
        ax = pn.plot_control_mesh(None, points, label=options.label)
        ax = pn.plot_iso_lines(ax, block)
        ax.set_title('control mesh and iso lines (blue)'
                     ' in Greville abscissae coordinates')
        ax.axis('equal')

        points = bcps
        ax = pn.plot_bezier_mesh(None, points, bconn, block.degree,
                                 label=options.label)
        ax = pn.plot_iso_lines(ax, block)
        ax.set_title('Bezier mesh and iso lines (blue)'
                     ' in Greville abscissae coordinates')
        ax.axis('equal')

        pn.plt.rcParams['lines.linewidth'] = 3

        line = block.extract(0, 0)
        if dim == 3:
            line = line.extract(0, 0)
        ax = pn.plot_nurbs_basis_1d(None, line, n_points=1000,
                                    legend=options.label)
        ax.set_xlabel('last parametric coordinate')
        ax.set_title('1D NURBS basis')

        ax = pn.plot_nurbs_basis_1d(None, line, n_points=1000, x_axis=dim-1,
                                    legend=options.label)
        ax.set_xlabel('last physical coordinate')
        ax.set_title('1D NURBS basis')

        pn.plt.show()

if __name__ == '__main__':
    main()
