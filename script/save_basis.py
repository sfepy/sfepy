#!/usr/bin/env python
"""
Save polynomial basis on reference elements for visualization.
"""
from optparse import OptionParser
import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import get_print_info
from sfepy.fem import Mesh
from sfepy.fem.geometry_element import GeometryElement
from sfepy.fem.poly_spaces import PolySpace
from sfepy.fem.linearizer import create_output

usage = '%prog [options]\n' + __doc__.rstrip()

help = {
    'basis' :
    'name of the FE basis [default: %default]',
    'derivative' :
    'save d-th derivative of FE basis, can be 0 or 1 [default: %default]',
    'max_order' :
    'maximum order of polynomials [default: %default]',
    'geometry' :
    'reference element geometry, one of "2_3", "2_4", "3_4", "3_8"'
    ' [default: %default]',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-b', '--basis', metavar='name',
                      action='store', dest='basis',
                      default='lagrange', help=help['basis'])
    parser.add_option('-d', '--derivative', metavar='d', type=int,
                      action='store', dest='derivative',
                      default=0, help=help['derivative'])
    parser.add_option('-n', '--max-order', metavar='order', type=int,
                      action='store', dest='max_order',
                      default=2, help=help['max_order'])
    parser.add_option('-g', '--geometry', metavar='name',
                      action='store', dest='geometry',
                      default='2_4', help=help['geometry'])
    options, args = parser.parse_args()

    dim, n_ep = int(options.geometry[0]), int(options.geometry[2])
    output('reference element geometry:')
    output('  dimension: %d, vertices: %d' % (dim, n_ep))

    output('polynomial space:', options.basis)

    output('max. order:', options.max_order)

    gel = GeometryElement(options.geometry)
    gps = PolySpace.any_from_args(None, gel, 1,
                                  base=options.basis)
    ps = PolySpace.any_from_args(None, gel, options.max_order,
                                 base=options.basis)

    n_digit, _format = get_print_info(ps.n_nod, fill='0')
    name_template = 'bf_%s.vtk' % _format
    for ip in range(ps.n_nod):
        output('shape function %d...' % ip)

        def eval_dofs(iels, rx, bf):
            if options.derivative == 0:
                rvals = bf[None, :, ip:ip+1]

            else:
                bfg = ps.eval_base(rx, diff=True)
                rvals = bfg[None, ..., ip]

            return rvals

        def eval_coors(iels, rx):
            bf = gps.eval_base(rx).squeeze()
            coors = nm.dot(bf, gel.coors)[None, ...]
            return coors

        (level, coors, conn,
         vdofs, mat_ids) = create_output(eval_dofs, eval_coors, 1,
                                         ps, min_level=2, max_level=5,
                                         eps=1e-3)
        out = {
            'bf' : Struct(name='output_data',
                          mode='vertex', data=vdofs,
                          var_name='bf', dofs=None)
        }

        mesh = Mesh.from_data('bf_mesh', coors, None, [conn], [mat_ids],
                              [options.geometry])

        name = name_template % ip
        mesh.write(name, out=out)

        output('...done (%s)' % name)

if __name__ == '__main__':
    main()
