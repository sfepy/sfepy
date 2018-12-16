import numpy as nm

from numpy import newaxis as nax


from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.base.base import Struct


class CanonicalPolySPace(PolySpace):

    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):

        if isinstance(coors, (int, float)):
            n = 1
        else:
            n = len(coors)
        values = nm.ones((n, self.order + 1, 1))
        for i in range(1, self.order + 1):
            values[:, i] = coors * values[:, i-1]
        return values


class LegendrePolySpace(PolySpace):
    """
    Legendre hierarchical polynomials basis, over [0, 1] domain
    use transform x = (y + 1)/2 to get basis over [0, 1]
    """

    def __init__(self, name, geometry, order):
        from toolz import map, reduce
        from operator import add, mul
        """
        Does not use init_context
        :param name:
        :param geometry: so far only 1_2 supported
        :param order:
        :param init_context: not used!
        """
        # TODO how is PolySpace supposed to look and work?
        # FIXME - complete LegendrePolySpace

        PolySpace.__init__(self, name, geometry, order)

        self.n_v = geometry.n_vertex,
        self.dim = geometry.dim
        self.n_nod = int(reduce(mul, map(lambda i: order + i + 1, range(self.dim))) /
                         reduce(mul, range(1, self.dim+1)))  # number of DOFs per element


    funs = [lambda x: 1,
            lambda x: x,
            lambda x: (3*x**2 - 1)/2,
            lambda x: (5*x**3 - 3*x)/2,
            lambda x: (35*x**4 - 30*x**2 + 3)/8,
            lambda x: (63*x**5 - 70*x**3 + 15*x)/8
            ]

    def legendreP(self, coors):
        """
        :param coors: coordinates, preferably in interval [0, 1] for which this basisi is intented
        :return: values in coors of all the legendre polynomials up to self.order
        """
        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)
        values = nm.ones(sh + (self.order + 1,))

        if self.order == 0:
            return values

        values[..., 1] = 2 * coors - 1
        for i in range(2, self.order + 1):
            # values[..., i] = ((2*i + 1) * coors * values[..., i-1] - i * values[..., i-2]) / (i + 1)
            # FIXME transform recursive formula
            values[..., i] = self.get_nth_fun(i)(coors)

        return values

    def jacobiP(self, coors, alpha, beta):
        """
        Evals very general jacobi polynomials for all orders from 0 to self.order
        Used formula works only for alpha, beta > -1 and alpha+beta != -1

        TODO move to JacobiPolySpace

        :param coors:
        :param alpha:
        :param beta:
        :return:
        """
        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)
        PL = nm.ones(sh + (self.order + 1,))

        from scipy.special import gamma
        gamma0 = 2 ** (alpha + beta + 1) / (alpha + beta + 1) * gamma(alpha + 1) * \
                 gamma(beta + 1) / gamma(alpha + beta + 1)

        PL[..., 0] = 1 / nm.sqrt(gamma0)
        if self.order == 0:
            return PL

        gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0

        PL[..., 1] = ((alpha + beta + 2) * coors / 2 + (alpha - beta) / 2) / nm.sqrt(gamma1)
        if self.order == 1:
            return PL

        aold = 2 / (2 + alpha + beta) * nm.sqrt((alpha + 1) * (beta + 1) / (alpha + beta + 3))

        for i in range(self.order - 1):
            # FIXME transform recursive formula to interval [0, 1]
            h1 = 2 * i + alpha + beta;
            anew = 2 / (h1 + 2) * nm.sqrt((i + 1) * (i + 1 + alpha + beta) * (i + 1 + alpha) *
                                          (i + 1 + beta) / (h1 + 1) / (h1 + 3))
            bnew = - (alpha ^ 2 - beta ^ 2) / h1 / (h1 + 2)
            PL[..., i + 2] = 1 / anew * (-aold * PL[..., i] + (coors - bnew) * PL[..., i + 1])
            aold = anew

        return PL


    def _eval_base(self, coors, diff=0, ori=None,
                       suppress_errors=False, eps=1e-15):
        return self.legendreP(coors)

    def get_nth_fun(self, n):
        """
        Convenience function for testing
        :param n: 0,1 , 2, 3, ...
        :return: n-th function of the legendre basis
        """

        if n < 6:
            return lambda x: self.funs[n](2*x - 1)
        else:
            from scipy.misc import comb as comb

            def fun(x):
                val = 0
                for k in range(n):
                    val = val + comb(n, k) * comb(n + k, k) * (((2*x-1)-1)/2.)**k

            return fun


class LegendreSimplexPolySpace(LegendrePolySpace):
    name = "legendre_simplex"

    def __init__(self, name, geometry, order):
        LegendrePolySpace.__init__(self, name, geometry, order)

        # self.nodes = nm.array([[1, 0], [0, 1]])
        # self.nts = nm.array([[0, 0], [0, 1]])
        # self.node_coors = nm.array([[0.], [1.]])


    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):
        """
        expects coors to be in shape (..., dim),
        returns values of the simplex basis in shape
            (coors.shape[:-1], n_el_nod)

        :param coors:
        :param diff:
        :param ori:
        :param suppress_errors:
        :param eps:
        :return:
        """
        porder = self.order + 1
        polyvals = nm.zeros(coors.shape + (porder,))
        values = nm.zeros(coors.shape[:-1] + (self.n_nod,))
        for i in range(self.dim):
            polyvals[..., i, :] = self.legendreP(coors[..., i])

        sq2 = nm.sqrt(2)
        # TODO so far only for 1 and 2D
        for i in range(porder):
            for j in range(porder - i):
                m = int(j + porder*i - i/2 * (i - 1))
                # reduce(mul, [polyvals[..., d, ind] for (d, ind) in zip(range(self.dim), (i, j, k))])
                values[..., m] = sq2 * polyvals[..., 0, i] * polyvals[..., 1, j]*(1 - coors[..., 1])**2

        return values


class LegendreTensorProductPolySpace(LegendrePolySpace):
    name = "legendre_tensor_product"

    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):

        # P_i(x) * P_j(y)

        return None


if __name__ == '__main__':
    from matplotlib import pylab as plt

    coors = nm.linspace(0, 1)[:, nax]
    geometry = Struct(n_vertex=2,
                      dim=1,
                      coors=coors.copy())

    # bs = CanonicalPolySPace('primb', geometry, 5)
    # vals = bs.eval_base(coors)

    bs = LegendrePolySpace('legb', geometry, 4)
    Legvals = bs.eval_base(coors)

    # plt.figure("Primitive polyspace")
    # plt.plot(nm.linspace(-1, 1), vals[: ,: ,0])

    plt.figure("Legendre polyspace")
    plt.plot(coors, Legvals[:, 0, :])
    plt.plot(coors, bs.get_nth_fun(2)(coors), "--")
    plt.plot(coors, bs.get_nth_fun(3)(coors), "--")
    plt.plot(coors, bs.get_nth_fun(4)(coors), "--")
    plt.plot([0, 1], [0, 0], 'k')

    plt.show()
    # geometry = Struct(n_vertex = 2,
    #              dim = 1,
    #              coors = self.bbox[:,0:1].copy())
