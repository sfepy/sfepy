import numpy as nm


# from dg_terms import unravel_sol

# TODO move the into dg_field
def get_unraveler(n_el_nod, n_cell):

    def unravel(u):
        ur = nm.zeros((n_cell, n_el_nod, 1))
        for i in range(n_el_nod):
            ur[:, i] = u[n_cell * i: n_cell * (i + 1), None]
        return ur

    return unravel


def get_raveler(n_el_nod, n_cell):

    def ravel(u):
        ur = nm.zeros((n_cell * n_el_nod, 1))
        for i in range(n_el_nod):
            ur[n_cell * i: n_cell * (i + 1)] = u[:, i]
        return ur

    return ravel


def moment_limiter_1D(u):
    """
    Krivodonova(2007): Limiters for high-order discontinuous Galerkin methods

    :param u: solution at time step n in shape
    (order, n_space_nod)
    :return: limited solution
    """

    def minmod(a, b, c):
        seq = (nm.sign(a) == nm.sign(b)) & (nm.sign(b) == nm.sign(c))

        res = nm.zeros(nm.shape(a))
        res[seq] = nm.sign(a[seq]) * nm.minimum.reduce([nm.abs(b[seq]),
                                                        nm.abs(a[seq]),
                                                        nm.abs(c[seq])])

        return res

    idx = nm.arange(nm.shape(u[0, 1:-1])[0])
    nu = nm.copy(u)
    for l in range(1, 0, -1):
        # tODO is limiter right this way?
        tilu = minmod(nu[l, 1:-1][idx],
                      nu[l - 1, 2:][idx] - nu[l - 1, 1:-1][idx],
                      nu[l - 1, 1:-1][idx] - nu[l - 1, :-2][idx])
        idx = tilu != nu
        nu[l, 1:-1][idx] = tilu[idx]
    return nu