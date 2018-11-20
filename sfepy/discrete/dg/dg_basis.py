import numpy as nm

from numpy import newaxis as nax


from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.base.base import Struct

class PrimitivePolySPace(PolySpace):

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
    Legendre hierarchical polynomials basis, over [-1, 1] domain
    use transform y = 2*x-1 to get basis over [0, 1]
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

        # "C:\Users\Lucia Moa\Python_Projects\sfepy\sfepy\discrete\fem\mappings.py", line 99, in get_mapping
        # cmap.describe(self.coors, self.conn, bf_g, ebf_g, weights)
        #
        # "sfepy\discrete\common\extmods\mappings.pyx", line 104, in sfepy.discrete.common.extmods.mappings.CMapping.describe
        # ValueError: ccore error(seeabove)

        PolySpace.__init__(self, name, geometry, order)

        n_v, dim = geometry.n_vertex, geometry.dim
        self.n_nod = (order + 1) ** dim

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
        :param coors:
        :param diff: not supported!
        :param ori: not supported!
        :param suppress_errors:
        :param eps: ???
        :return: values in coors of all the basis function up to order
        shape = (order + 1, ) + coors.shape() or (order + 1, 1) of coors is scalar
        """
        # TODO this does not correspond to _eval_base of other basis!
        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)
        values = nm.ones((self.order + 1,) + sh)
        values[1, :] = coors
        for i in range(2, self.order + 1):
            values[i, :] = ((2*i + 1) * coors * values[i-1, :] - i * values[i-2, :]) / (i + 1)

        return values

    def get_nth_fun(self, n):
        """
        Convenience function for testing
        :param n:
        :return: n-th function of the legendre basis
        """

        if n < 6:
            return self.funs[n]
        else:
            from scipy.misc import comb as comb

            def fun(x):
                val = 0
                for k in range(n):
                    val = val + comb(n, k) * comb(n + k, k) * ((x-1)/2.)**k

            return fun


if __name__ == '__main__':
    from matplotlib import pylab as plt

    coors = nm.linspace(-1, 1)[:, nax]
    geometry = Struct(n_vertex=2,
                      dim=1,
                      coors=coors.copy())

    bs = PrimitivePolySPace('primb', geometry, 5)
    vals = bs.eval_base(coors)

    bs = LegendrePolySpace('legb', geometry, 2)
    Legvals = bs.eval_base(coors)

    # plt.figure("Primitive polyspace")
    # plt.plot(nm.linspace(-1, 1), vals[: ,: ,0])

    plt.figure("Legendre polyspace")
    plt.plot(nm.linspace(-1, 1), Legvals[:, :, 0].T)
    plt.show()
    # geometry = Struct(n_vertex = 2,
    #              dim = 1,
    #              coors = self.bbox[:,0:1].copy())
