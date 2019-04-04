import numpy as nm
from sfepy.discrete.dg.dg_field import get_raveler, get_unraveler
from sfepy.base.base import (get_default, output, assert_,
                             Struct, IndexedStruct)

MACHINE_EPS = 1e-40

def minmod(a, b, c):
    """
    Minmod function of three variables, returns:
     _/ 0           , where sign(a) != sign(b) != sign(c)
      \ min(a,b,c)  , elsewhere
    :param a:
    :param b:
    :param c:
    :return:
    """
    seq = (nm.sign(a) == nm.sign(b)) & (nm.sign(b) == nm.sign(c))

    res = nm.zeros(nm.shape(a))
    res[seq] = nm.sign(a[seq]) * nm.minimum.reduce([nm.abs(b[seq]),
                                                    nm.abs(a[seq]),
                                                    nm.abs(c[seq])])

    return res

class DGLimiter:

    name = "abstract DG limiter"

    def __init__(self, n_el_nod, n_cell, verbose=False):
        self.n_el_nod = n_el_nod
        self.n_cell = n_cell
        self.ravel = get_raveler(n_el_nod, n_cell)
        self.unravel = get_unraveler(n_el_nod, n_cell)
        self.verbose = verbose

    def __call__(self, u):
        raise NotImplementedError("Called abstract limiter")


class IdentityLimiter(DGLimiter):

    name = "identity"

    def __call__(self, u):
        if verbose: ouput(self.name)
        return u


class Moment1DLimiter(DGLimiter):
    """
    Krivodonova(2007): Limiters for high-order discontinuous Galerkin methods
    """
    name = "moment_1D_limiter"

    def __call__(self, u):
        """"
        :param u: solution at time step n in shape
        (order, n_space_nod)
        :return: limited solution
        """
        # for convenience do not try to limit FV
        if self.n_el_nod == 1:
            return u
        u = self.unravel(u).swapaxes(0, 1)

        idx = nm.arange(nm.shape(u[0, 1:-1])[0])
        nu = nm.copy(u)
        for l in range(self.n_el_nod - 1, 0, -1):
            tilu = minmod(nu[l, 1:-1][idx],
                          nu[l - 1, 2:][idx] - nu[l - 1, 1:-1][idx],
                          nu[l - 1, 1:-1][idx] - nu[l - 1, :-2][idx])
            idx = abs(tilu - nu[l ,1:-1][idx]) > MACHINE_EPS
            if self.verbose:
                output(self.name + " limiting in {} cells of {} :".format(sum(idx[:,0]), self.n_cell))
                # output(nm.where(idx[:, 0]))
            nu[l, 1:-1][idx] = tilu[idx]


        return self.ravel(nu.swapaxes(0, 1))[:, 0]
