#!/usr/bin/env python
"""
Plot conditions numbers w.r.t. polynomial approximation order of reference
element matrices for various FE polynomial spaces (bases).
"""
from optparse import OptionParser
import numpy as nm
import matplotlib.pyplot as plt

from sfepy import data_dir
from sfepy.base.base import output, assert_
from sfepy.fem import Mesh, Domain, Field, FieldVariable, Material, Integral
from sfepy.terms import Term
from sfepy.solvers import eig

usage = '%prog [options]\n' + __doc__.rstrip()

help = {
    'basis' :
    'name of the FE basis [default: %default]',
    'max_order' :
    'maximum order of polynomials [default: %default]',
    'matrix_type' :
    'matrix type, one of "elasticity", "laplace" [default: %default]',
    'geometry' :
    'reference element geometry, one of "2_3", "2_4", "3_4", "3_8"'
    ' [default: %default]',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-b', '--basis', metavar='name',
                      action='store', dest='basis',
                      default='lagrange', help=help['basis'])
    parser.add_option('-n', '--max-order', metavar='order', type=int,
                      action='store', dest='max_order',
                      default=5, help=help['max_order'])
    parser.add_option('-m', '--matrix', metavar='type',
                      action='store', dest='matrix_type',
                      default='laplace', help=help['matrix_type'])
    parser.add_option('-g', '--geometry', metavar='name',
                      action='store', dest='geometry',
                      default='2_3', help=help['geometry'])
    options, args = parser.parse_args()

    dim, n_ep = int(options.geometry[0]), int(options.geometry[2])
    output('reference element geometry:')
    output('  dimension: %d, vertices: %d' % (dim, n_ep))

    n_c = {'laplace' : 1, 'elasticity' : dim}[options.matrix_type]

    output('matrix type:', options.matrix_type)
    output('number of variable components:',  n_c)

    output('polynomial space:', options.basis)

    output('max. order:', options.max_order)

    mesh = Mesh.from_file(data_dir + '/meshes/elements/%s_1.mesh'
                          % options.geometry)
    domain = Domain('domain', mesh)
    omega = domain.create_region('Omega', 'all')

    eps = nm.finfo(nm.float64).eps
    eps100 = 100.0 * eps

    orders = nm.arange(1, options.max_order + 1, dtype=nm.int)
    conds = []

    for order in orders:
        output('order:', order, '...')

        field = Field('fu', nm.float64, n_c, omega,
                      space='H1', poly_space_base=options.basis,
                      approx_order=order)

        to = field.get_true_order()
        quad_order = 2 * (max(to - 1, 0))
        output('quadrature order:', quad_order)


        u = FieldVariable('u', 'unknown', field, n_c)
        v = FieldVariable('v', 'test', field, n_c, primary_var_name='u')

        m = Material('m', lam=1.0, mu=1.0)

        integral = Integral('i', order=quad_order)

        if options.matrix_type == 'laplace':
            term = Term.new('dw_laplace(m.mu, v, u)',
                            integral, omega, m=m, v=v, u=u)

        else:
            assert_(options.matrix_type == 'elasticity')
            term = Term.new('dw_lin_elastic_iso(m.lam, m.mu, v, u)',
                            integral, omega, m=m, v=v, u=u)

        term.setup()

        mtx, iels = term.evaluate(mode='weak', diff_var='u')
        mtx = mtx[0][0, 0]

        try:
            assert_(nm.max(nm.abs(mtx - mtx.T)) < 1e-12)

        except:
            from sfepy.base.base import debug; debug()

        output('matrix shape:', mtx.shape)

        eigs = eig(mtx, method='eig.sgscipy', eigenvectors=False)
        eigs.sort()

        # Zero 'true' zeros.
        eigs[nm.abs(eigs) < eps100] = 0.0

        ii = nm.where(eigs < 0.0)[0]
        if len(ii):
            output('matrix is not positive semi-definite!')

        output('eigs:\n', eigs)

        ii = nm.where(eigs > 0.0)[0]
        emin, emax = eigs[ii[[0, -1]]]

        output('min:', emin, 'max:', emax)

        cond = emax / emin
        conds.append(cond)

        output('condition number:', cond)

        output('...done')

    plt.figure(1)
    plt.semilogy(orders, conds)
    plt.xticks(orders)
    plt.xlabel('polynomial order')
    plt.ylabel('condition number')
    plt.grid()

    plt.figure(2)
    plt.loglog(orders, conds)
    plt.xticks(orders)
    plt.xlabel('polynomial order')
    plt.ylabel('condition number')
    plt.grid()

    plt.show()

if __name__ == '__main__':
    main()
