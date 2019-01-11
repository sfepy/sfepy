import sympy as sp
from sympy import var, symbols
from sympy import jacobi_poly as jacobi_P
from sympy import expand_mul as expand
from sympy import cancel, simplify
from sympy import Matrix

from toolz import reduce, map
from operator import mul

import numpy as np

from dg_basis import iter_by_order

a, b = symbols("a, b")
r, s = symbols("r, s")
x, y = symbols("x, y")
order = 10
dim = 2
n_el_nod = int(reduce(mul, map(lambda i: order + i + 1, range(dim))) /
               reduce(mul, range(1, dim+1)))  # number of DOFs per element



simplexP = []
exponentM = np.zeros((n_el_nod, 3), dtype=np.int32)
coefM = np.zeros((n_el_nod, n_el_nod))
exponentList = []
for m, idx in enumerate(iter_by_order(order, dim)):
    # print(m, idx)
    exponentM[m, :dim] = idx
    pa = jacobi_P(idx[0], 0, 0, a)
    pb = jacobi_P(idx[1], 2*idx[0] + 1, 0, b)*(1 - b)**idx[0]
    # print("P_{} = {}".format(m, pa*pb))
    polrs = cancel((pa*pb).subs(b , s).subs(
                    a, 2 * (1 + r) / (1 - s) - 1))
    # print("P_{} = {}".format(m, polrs))
    polxy = expand(polrs.subs(r, 2*x - 1).subs(s, 2*y - 1))
    simplexP.append(simplify(polxy))
    exponentList.append(x**idx[0]*y**idx[1])
    for j, exponent in enumerate(exponentList):
        coefM[m, j] = simplexP[m].as_coefficients_dict()[exponent]
    print("P_{}{} = {}".format(m, idx, simplexP[m]))
    print()
np.savetxt("legendre2D_simplex_expos.txt", exponentM, fmt="%d")
np.savetxt("legendre2D_simplex_coefs.txt", coefM, fmt="%d") # tODO are coefs always integers?


