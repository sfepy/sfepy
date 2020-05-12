#!/usr/bin/env python
"""
Generate simplex DG basis data.
"""
from argparse import ArgumentParser
from sympy import symbols
from sympy import jacobi_poly as jacobi_P
from sympy import expand_mul as expand
from sympy import cancel, simplify

from toolz import reduce, map
from operator import mul

import numpy as np

from sfepy.base.ioutils import InDir
from sfepy.discrete.dg.dg_poly_spaces import iter_by_order

helps = {
    'max_order' :
        'maximum order of polynomials [default: %(default)s]',
    'output_dir':
        'output directory',
}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('-m', '--max-order', metavar='order', type=int,
                        action='store', dest='max_order',
                        default=10, help=helps['max_order'])
    parser.add_argument('output_dir', help=helps['output_dir'])
    options = parser.parse_args()

    indir = InDir(options.output_dir)

    a, b = symbols("a, b")
    r, s = symbols("r, s")
    x, y = symbols("x, y")
    order = options.max_order
    dim = 2
    n_el_nod = int(reduce(mul, map(lambda i: order + i + 1, range(dim))) /
                   reduce(mul, range(1, dim + 1)))  # number of DOFs per element

    simplexP = []
    exponentM = np.zeros((n_el_nod, 3), dtype=np.int32)
    coefM = np.zeros((n_el_nod, n_el_nod))
    exponentList = []
    for m, idx in enumerate(iter_by_order(order, dim)):
        # print(m, idx)
        exponentM[m, :dim] = idx
        pa = jacobi_P(idx[0], 0, 0, a)
        pb = jacobi_P(idx[1], 2 * idx[0] + 1, 0, b) * (1 - b) ** idx[0]
        # print("P_{} = {}".format(m, pa*pb))
        polrs = cancel((pa * pb).subs(b, s).subs(
                a, 2 * (1 + r) / (1 - s) - 1))
        # print("P_{} = {}".format(m, polrs))
        polxy = expand(polrs.subs(r, 2 * x - 1).subs(s, 2 * y - 1))
        # polxy = expand(polrs.subs(r, x).subs(s, y))

        simplexP.append(simplify(polxy))
        exponentList.append(x ** idx[0] * y ** idx[1])
        for j, exponent in enumerate(exponentList):
            coefM[m, j] = simplexP[m].as_coefficients_dict()[exponent]
        print("P_{}{} = {}".format(m, idx, simplexP[m]))
        print()

    np.savetxt(indir("legendre2D_simplex_expos.txt"), exponentM, fmt="%d")
    # TODO are coefs always integers?
    np.savetxt(indir("legendre2D_simplex_coefs.txt"), coefM, fmt="%d")


if __name__ == '__main__':
    main()
