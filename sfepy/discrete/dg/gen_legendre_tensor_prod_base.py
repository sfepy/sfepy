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

r, s = symbols("r, s")
x, y = symbols("x, y")
order = 10
dim = 2
n_el_nod = int(reduce(mul, map(lambda i: order + i + 1, range(dim))) /
               reduce(mul, range(1, dim+1)))  # number of DOFs per element



tensorP = []
exponentM = np.zeros((n_el_nod, 3), dtype=np.int32)
coefM = np.zeros((n_el_nod, n_el_nod))
exponentList = []
for m, idx in enumerate(iter_by_order(order, dim)):
    # print(m, idx)
    exponentM[m, :dim] = idx
    # print("P_{} = {}".format(m, pa*pb))
    polrs = jacobi_P(idx[0], 0, 0, r) * jacobi_P(idx[1], 0, 0, s)
    # print("P_{} = {}".format(m, polrs))
    # polxy = expand(polrs.subs(r, 2*x - 1).subs(s, 2*y - 1))
    polxy = expand(polrs.subs(r, x).subs(s, y))

    tensorP.append(simplify(polxy))
    exponentList.append(x**idx[0]*y**idx[1])
    for j, exponent in enumerate(exponentList):
        coefM[m, j] = tensorP[m].as_coefficients_dict()[exponent]
    print("P_{}{} = {}".format(m, idx, tensorP[m]))
    print()
np.savetxt("legendre2D_tensor_expos.txt", exponentM, fmt="%d")
np.savetxt("legendre2D_tensor_coefs.txt", coefM, fmt="%d")  # TODO are coefs always integers?


