#!/usr/bin/env python
"""
Plot conditions numbers w.r.t. polynomial approximation order of reference
element matrices for various FE polynomial spaces (bases).
"""
from __future__ import absolute_import
import sys
sys.path.append('.')
from argparse import ArgumentParser
import os.path as op
from functools import partial
import numpy as nm
import matplotlib.pyplot as plt

from sfepy import data_dir
from sfepy.base.base import output, assert_
from sfepy.base.timing import Timer
from sfepy.discrete import FieldVariable, Material, Integral
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.terms import Term
from sfepy.solvers import eig
from sfepy.mechanics.matcoefs import stiffness_from_lame

helps = {
    'basis' :
    'name of the FE basis [default: %(default)s]',
    'max_order' :
    'maximum order of polynomials [default: %(default)s]',
    'matrix_type' :
    'matrix type [default: %(default)s]',
    'geometry' :
    'reference element geometry, one of "2_3", "2_4", "3_4", "3_8"'
    ' [default: %(default)s]',
    'output_dir' :
    'output directory',
    'no_show' :
    'do not show the figures',
}

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-b', '--basis', metavar='name',
                        action='store', dest='basis',
                        default='lagrange', help=helps['basis'])
    parser.add_argument('-n', '--max-order', metavar='order', type=int,
                        action='store', dest='max_order',
                        default=10, help=helps['max_order'])
    parser.add_argument('-m', '--matrix', action='store', dest='matrix_type',
                        choices=['laplace', 'elasticity', 'smass', 'vmass'],
                        default='laplace', help=helps['matrix_type'])
    parser.add_argument('-g', '--geometry', metavar='name',
                        action='store', dest='geometry',
                        default='2_4', help=helps['geometry'])
    parser.add_argument('-o', '--output-dir', metavar='path',
                        action='store', dest='output_dir',
                        default=None, help=helps['output_dir'])
    parser.add_argument('--no-show',
                        action='store_false', dest='show',
                        default=True, help=helps['no_show'])
    options = parser.parse_args()

    dim, n_ep = int(options.geometry[0]), int(options.geometry[2])
    output('reference element geometry:')
    output('  dimension: %d, vertices: %d' % (dim, n_ep))

    n_c = {'laplace' : 1, 'elasticity' : dim,
           'smass' : 1, 'vmass' : dim}[options.matrix_type]

    output('matrix type:', options.matrix_type)
    output('number of variable components:',  n_c)

    output('polynomial space:', options.basis)

    output('max. order:', options.max_order)

    mesh = Mesh.from_file(data_dir + '/meshes/elements/%s_1.mesh'
                          % options.geometry)
    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')

    orders = nm.arange(1, options.max_order + 1, dtype=nm.int32)
    conds = []

    for order in orders:
        output('order:', order, '...')

        field = Field.from_args('fu', nm.float64, n_c, omega,
                                approx_order=order,
                                space='H1', poly_space_basis=options.basis)

        quad_order = 2 * field.approx_order
        output('quadrature order:', quad_order)

        integral = Integral('i', order=quad_order)
        qp, _ = integral.get_qp(options.geometry)
        output('number of quadrature points:', qp.shape[0])

        u = FieldVariable('u', 'unknown', field)
        v = FieldVariable('v', 'test', field, primary_var_name='u')

        m = Material('m', D=stiffness_from_lame(dim, 1.0, 1.0))

        if options.matrix_type == 'laplace':
            term = Term.new('dw_laplace(v, u)',
                            integral, omega, v=v, u=u)
            n_zero = 1

        elif options.matrix_type == 'elasticity':
            term = Term.new('dw_lin_elastic(m.D, v, u)',
                            integral, omega, m=m, v=v, u=u)
            n_zero = (dim + 1) * dim // 2

        elif options.matrix_type in ('smass', 'vmass'):
            term = Term.new('dw_dot(v, u)',
                            integral, omega, v=v, u=u)
            n_zero = 0

        term.setup()

        output('assembling...')
        timer = Timer(start=True)
        mtx, iels = term.evaluate(mode='weak', diff_var='u')
        output('...done in %.2f s' % timer.stop())
        mtx = mtx[0, 0]

        try:
            assert_(nm.max(nm.abs(mtx - mtx.T)) < 1e-10)

        except:
            from sfepy.base.base import debug; debug()

        output('matrix shape:', mtx.shape)

        eigs = eig(mtx, solver_kind='eig.sgscipy', eigenvectors=False)
        eigs.sort()

        # Zero 'true' zeros.
        eigs[:n_zero] = 0.0

        ii = nm.where(eigs < 0.0)[0]
        if len(ii):
            output('matrix is not positive semi-definite!')

        ii = nm.where(eigs[n_zero:] < 1e-12)[0]
        if len(ii):
            output('matrix has more than %d zero eigenvalues!' % n_zero)

        output('smallest eigs:\n', eigs[:10])

        ii = nm.where(eigs > 0.0)[0]
        emin, emax = eigs[ii[[0, -1]]]

        output('min:', emin, 'max:', emax)

        cond = emax / emin
        conds.append(cond)

        output('condition number:', cond)

        output('...done')

    if options.output_dir is not None:
        indir = partial(op.join, options.output_dir)

    else:
        indir = None

    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 3

    fig, ax = plt.subplots()
    ax.semilogy(orders, conds)
    ax.set_xticks(orders)
    ax.set_xticklabels(orders)
    ax.set_xlabel('polynomial order')
    ax.set_ylabel('condition number')
    ax.set_title(f'{options.basis.capitalize()} basis')
    ax.grid()
    plt.tight_layout()
    if indir is not None:
        fig.savefig(indir(f'{options.basis}-{options.matrix_type}-'
                          f'{options.geometry}-{options.max_order}-xlin.png'),
                    bbox_inches='tight')

    fig, ax = plt.subplots()
    ax.loglog(orders, conds)
    ax.set_xticks(orders)
    ax.set_xticklabels(orders)
    ax.set_xlabel('polynomial order')
    ax.set_ylabel('condition number')
    ax.set_title(f'{options.basis.capitalize()} basis')
    ax.grid()
    plt.tight_layout()
    if indir is not None:
        fig.savefig(indir(f'{options.basis}-{options.matrix_type}-'
                          f'{options.geometry}-{options.max_order}-xlog.png'),
                    bbox_inches='tight')

    if options.show:
        plt.show()

if __name__ == '__main__':
    main()
