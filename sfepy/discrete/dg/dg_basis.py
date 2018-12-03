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

    def __init__(self, name, geometry, order, init_context=True):
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

        n_v, dim = geometry.n_vertex, geometry.dim
        self.n_nod = (order + 1) ** dim  # number of DOFs per element

        self.nodes = nm.array([[1, 0], [0, 1]])
        self.nts = nm.array([[0, 0], [0, 1]])
        self.node_coors = nm.array([[-1.], [1.]])

        self.eval_ctx = None

    funs = [lambda x: 1,
            lambda x: x,
            lambda x: (3*x**2 - 1)/2,
            lambda x: (5*x**3 - 3*x)/2,
            lambda x: (35*x**4 - 30*x**2 + 3)/8,
            lambda x: (63*x**5 - 70*x**3 + 15*x)/8
            ]

    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):
        """
        Numpy valuation of basis functions
        :param coors: coordinates, preferably in interval [0, 1] for which this basisi is intented
        :param diff: not supported!
        :param ori: not supported!
        :param suppress_errors:
        :param eps: ???
        :return: values in coors of all the basis function up to order

        """
        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)
        values = nm.ones(sh + (self.order + 1,))

        if self.order == 0:
            return values

        values[..., 1] = 2*coors-1
        for i in range(2, self.order + 1):
            # values[..., i] = ((2*i + 1) * coors * values[..., i-1] - i * values[..., i-2]) / (i + 1)
            # FIXME tranform recursive formula
            values[..., i] = self.get_nth_fun(i)(coors)

        return values

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
