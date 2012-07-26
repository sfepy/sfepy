#!/usr/bin/env python
"""
Save polynomial basis on reference elements or on a mesh for visualization.
"""
from optparse import OptionParser
import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import get_print_info
from sfepy.fem import Mesh, Domain, Field, FieldVariable
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
    'mesh' :
    'name of the mesh file - alternative to --geometry [default: %default]',
    'lin_options' :
    'linearizer options [default: %default]',
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
    parser.add_option('-m', '--mesh', metavar='mesh',
                      action='store', dest='mesh',
                      default=None, help=help['mesh'])
    parser.add_option('-l', '--lin-options', metavar='options',
                      action='store', dest='lin_options',
                      default='min_level=2,max_level=5,eps=1e-3',
                      help=help['lin_options'])
    options, args = parser.parse_args()

    output('polynomial space:', options.basis)
    output('max. order:', options.max_order)

    lin = Struct(kind='adaptive', min_level=2, max_level=5, eps=1e-3)
    for opt in options.lin_options.split(','):
        key, val = opt.split('=')
        setattr(lin, key, eval(val))

    if options.mesh is None:
        dim, n_ep = int(options.geometry[0]), int(options.geometry[2])
        output('reference element geometry:')
        output('  dimension: %d, vertices: %d' % (dim, n_ep))

        gel = GeometryElement(options.geometry)
        gps = PolySpace.any_from_args(None, gel, 1,
                                      base=options.basis)
        ps = PolySpace.any_from_args(None, gel, options.max_order,
                                     base=options.basis)

        n_digit, _format = get_print_info(ps.n_nod, fill='0')
        name_template = 'bf_%s.vtk' % _format
        for ip in range(ps.n_nod):
            output('shape function %d...' % ip)

            def eval_dofs(iels, rx):
                if options.derivative == 0:
                    bf = ps.eval_base(rx).squeeze()
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
                                             ps, min_level=lin.min_level,
                                             max_level=lin.max_level,
                                             eps=lin.eps)
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

    else:
        mesh = Mesh.from_file(options.mesh)
        output('mesh geometry:')
        output('  dimension: %d, vertices: %d, elements: %d'
               % (mesh.dim, mesh.n_nod, mesh.n_el))

        domain = Domain('domain', mesh)
        omega = domain.create_region('Omega', 'all')
        field = Field.from_args('f', nm.float64, shape=1, region=omega,
                                approx_order=options.max_order,
                                poly_space_base=options.basis)
        var = FieldVariable('u', 'unknown', field, 1)

        output('dofs: %d' % var.n_dof)

        vec = nm.empty(var.n_dof, dtype=var.dtype)
        n_digit, _format = get_print_info(var.n_dof, fill='0')
        name_template = 'dof_%s.vtk' % _format
        for ip in range(var.n_dof):
            output('dof %d...' % ip)

            vec.fill(0.0)
            vec[ip] = 1.0

            var.data_from_any(vec)

            out = var.create_output(vec, linearization=lin)

            name = name_template % ip
            out['u'].mesh.write(name, out=out)

            output('...done (%s)' % name)

if __name__ == '__main__':
    main()
