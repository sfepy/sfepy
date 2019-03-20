import numpy as nm
from sfepy.discrete.dg.dg_field import get_raveler, get_unraveler


def minmod(a, b, c):
    seq = (nm.sign(a) == nm.sign(b)) & (nm.sign(b) == nm.sign(c))

    res = nm.zeros(nm.shape(a))
    res[seq] = nm.sign(a[seq]) * nm.minimum.reduce([nm.abs(b[seq]),
                                                    nm.abs(a[seq]),
                                                    nm.abs(c[seq])])

    return res


class DGLimiter:

    name = "abstract DG limiter"

    def __init__(self, n_el_nod, n_cell):
        self.n_el_nod = n_el_nod
        self.n_cell = n_cell
        self.ravel = get_raveler(n_el_nod, n_cell)
        self.unravel = get_unraveler(n_el_nod, n_cell)

    def __call__(self, u):
        raise NotImplementedError("Called abstract limiter")


class IdentityLimiter(DGLimiter):

    name = "identity DG limiter"

    def __call__(self, u):
        return u


class Moment1DLimiter(DGLimiter):
    """
    Krivodonova(2007): Limiters for high-order discontinuous Galerkin methods
    """
    name = "krivodonova moment 1D limiter"

    def __call__(self, u):
        """"
        :param u: solution at time step n in shape
        (order, n_space_nod)
        :return: limited solution
        """
        # for convenience do not try to limit FV
        if u.shape[0] == 1:
            return u

        u = self.unravel(u).swapaxes(0, 1)

        idx = nm.arange(nm.shape(u[0, 1:-1])[0])
        nu = nm.copy(u)
        n_el_nod = u.shape[0]
        for l in range(n_el_nod - 1, 0, -1):
            tilu = minmod(nu[l, 1:-1][idx],
                          nu[l - 1, 2:][idx] - nu[l - 1, 1:-1][idx],
                          nu[l - 1, 1:-1][idx] - nu[l - 1, :-2][idx])
            idx = tilu != nu
            nu[l, 1:-1][idx] = tilu[idx]

        return self.ravel(nu.swapaxes(0, 1))
