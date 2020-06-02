#!/usr/bin/env python
"""
Save polynomial basis on reference elements or on a mesh for visualization into
a given output directory.
"""
from __future__ import absolute_import
import sys
from six.moves import range
sys.path.append('.')
import os
from argparse import ArgumentParser
import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import get_print_info, ensure_path
from sfepy.discrete import FieldVariable, Variables, PolySpace
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete.fem.geometry_element import GeometryElement
from sfepy.discrete.fem.linearizer import create_output
from sfepy.discrete.fem.fields_base import create_expression_output

helps = {
    'basis' :
    'name of the FE basis [default: %(default)s]',
    'derivative' :
    'save d-th derivative of FE basis, can be 0 or 1 [default: %(default)s]',
    'max_order' :
    'maximum order of polynomials [default: %(default)s]',
    'geometry' :
    'reference element geometry, one of "2_3", "2_4", "3_4", "3_8"'
    ' [default: %(default)s]',
    'mesh' :
    ('name of the mesh file - alternative to --geometry '
     '[default: %(default)s]'),
    'permutations' :
    'list of geometry element permutations for each element, e.g. 0,1 is a'
    ' single permutation for two elements, 0,1,0,2,1,0 are three permutations'
    ' for two elements. Special value "all" can be used to save all possible'
    ' permutations for given reference element. Works only with --mesh option'
    ' [default: %(default)s]',
    'dofs' :
    'if given, save only the DOFs specified as a comma-separated list'
    ' [default: %(default)s]',
    'lin_options' :
    'linearizer options [default: %(default)s]',
    'plot_dofs' :
    'plot local and global DOF numberings, with --mesh option',
}

def get_dofs(dofs, n_total):
    if dofs is None:
        dofs = list(range(n_total))

    else:
        dofs = [int(ii) for ii in dofs.split(',')]

    return dofs

def save_basis_on_mesh(mesh, options, output_dir, lin,
                       permutations=None, suffix=''):
    if permutations is not None:
        mesh = mesh.copy()
        gel = GeometryElement(mesh.descs[0])
        perms = gel.get_conn_permutations()[permutations]
        conn = mesh.cmesh.get_cell_conn()
        n_el, n_ep = conn.num, gel.n_vertex
        offsets = nm.arange(n_el) * n_ep

        conn.indices[:] = conn.indices.take((perms + offsets[:, None]).ravel())

    domain = FEDomain('domain', mesh)

    omega = domain.create_region('Omega', 'all')
    field = Field.from_args('f', nm.float64, shape=1, region=omega,
                            approx_order=options.max_order,
                            poly_space_base=options.basis)
    var = FieldVariable('u', 'unknown', field)

    if options.plot_dofs:
        import sfepy.postprocess.plot_dofs as pd
        import sfepy.postprocess.plot_cmesh as pc
        ax = pc.plot_wireframe(None, mesh.cmesh)
        ax = pd.plot_global_dofs(ax, field.get_coor(), field.econn)
        ax = pd.plot_local_dofs(ax, field.get_coor(), field.econn)
        if options.dofs is not None:
            ax = pd.plot_nodes(ax, field.get_coor(), field.econn,
                               field.poly_space.nodes,
                               get_dofs(options.dofs, var.n_dof))
        pd.plt.show()

    output('dofs: %d' % var.n_dof)

    vec = nm.empty(var.n_dof, dtype=var.dtype)
    n_digit, _format = get_print_info(var.n_dof, fill='0')
    name_template = os.path.join(output_dir,
                                 'dof_%s%s.vtk' % (_format, suffix))
    for ip in get_dofs(options.dofs, var.n_dof):
        output('dof %d...' % ip)

        vec.fill(0.0)
        vec[ip] = 1.0

        var.set_data(vec)

        if options.derivative == 0:
            out = var.create_output(vec, linearization=lin)

        else:
            out = create_expression_output('ev_grad.ie.Elements(u)',
                                           'u', 'f', {'f' : field}, None,
                                           Variables([var]),
                                           mode='qp', verbose=False,
                                           min_level=lin.min_level,
                                           max_level=lin.max_level,
                                           eps=lin.eps)

        name = name_template % ip
        ensure_path(name)
        out['u'].mesh.write(name, out=out)

        output('...done (%s)' % name)

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-b', '--basis', metavar='name',
                        action='store', dest='basis',
                        default='lagrange', help=helps['basis'])
    parser.add_argument('-d', '--derivative', metavar='d', type=int,
                        action='store', dest='derivative',
                        default=0, help=helps['derivative'])
    parser.add_argument('-n', '--max-order', metavar='order', type=int,
                        action='store', dest='max_order',
                        default=2, help=helps['max_order'])
    parser.add_argument('-g', '--geometry', metavar='name',
                        action='store', dest='geometry',
                        default='2_4', help=helps['geometry'])
    parser.add_argument('-m', '--mesh', metavar='mesh',
                        action='store', dest='mesh',
                        default=None, help=helps['mesh'])
    parser.add_argument('--permutations', metavar='permutations',
                        action='store', dest='permutations',
                        default=None, help=helps['permutations'])
    parser.add_argument('--dofs', metavar='dofs',
                        action='store', dest='dofs',
                        default=None, help=helps['dofs'])
    parser.add_argument('-l', '--lin-options', metavar='options',
                        action='store', dest='lin_options',
                        default='min_level=2,max_level=5,eps=1e-3',
                        help=helps['lin_options'])
    parser.add_argument('--plot-dofs',
                        action='store_true', dest='plot_dofs',
                        default=False, help=helps['plot_dofs'])
    parser.add_argument('output_dir')
    options = parser.parse_args()

    output_dir = options.output_dir

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
        name_template = os.path.join(output_dir, 'bf_%s.vtk' % _format)
        for ip in get_dofs(options.dofs, ps.n_nod):
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
            ensure_path(name)
            mesh.write(name, out=out)

            output('...done (%s)' % name)

    else:
        mesh = Mesh.from_file(options.mesh)
        output('mesh geometry:')
        output('  dimension: %d, vertices: %d, elements: %d'
               % (mesh.dim, mesh.n_nod, mesh.n_el))

        if options.permutations:
            if options.permutations == 'all':
                from sfepy.linalg import cycle
                gel = GeometryElement(mesh.descs[0])
                n_perms = gel.get_conn_permutations().shape[0]
                all_permutations = [ii for ii in cycle(mesh.n_el * [n_perms])]

            else:
                all_permutations = [int(ii)
                                    for ii in options.permutations.split(',')]
                all_permutations = nm.array(all_permutations)
                np = len(all_permutations)
                all_permutations.shape = (np // mesh.n_el, mesh.n_el)

            output('using connectivity permutations:\n', all_permutations)

        else:
            all_permutations = [None]

        for ip, permutations in enumerate(all_permutations):
            if permutations is None:
                suffix = ''

            else:
                suffix = '_' + '_'.join('%d' % ii for ii in permutations)

            save_basis_on_mesh(mesh, options, output_dir, lin, permutations,
                               suffix)

if __name__ == '__main__':
    main()
