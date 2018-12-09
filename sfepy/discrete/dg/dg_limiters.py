import numpy as nm

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
        tilu = minmod(nu[l, 1:-1][idx],
                      nu[l - 1, 2:][idx] - nu[l - 1, 1:-1][idx],
                      nu[l - 1, 1:-1][idx] - nu[l - 1, :-2][idx])
        idx = tilu != nu
        nu[l, 1:-1][idx] = tilu[idx]
    return nu