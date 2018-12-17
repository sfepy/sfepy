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

        PolySpace.__init__(self, name, geometry, order)

        self.n_v = geometry.n_vertex,
        self.dim = geometry.dim
        self.n_nod = int(reduce(mul, map(lambda i: order + i + 1, range(self.dim))) /
                         reduce(mul, range(1, self.dim+1)))  # number of DOFs per element

    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):
        """
        Calls combine_polyvals or combine_polyvals_der to build multidimensional basis
        implement these methods in subclasses to get different basis.
        expects coors to be in shape
            (..., dim),
        Returns values of the basis in shape
            (1, coors.shape[:-1], 1, n_el_nod)
        or values of the gradient in shape
            (1, coors.shape[:-1], dim, n_el_nod)

        :param coors:
        :param diff:
        :param ori: TODO what is this for?
        :param suppress_errors:
        :param eps:
        :return:
        """
        porder = self.order + 1
        if diff:
            values = nm.zeros((1,) + coors.shape[:-1] + (self.dim, self.n_nod))  # number of values of derivative is equal to dimension
            polyvals = nm.zeros(coors.shape + (porder, 2))  # 2 is derivative order plus 1, so far we support only first derivative
            polyvals[..., 0] = self.legendreP(coors)
            polyvals[..., 1] = self.gradlegendreP(coors)

            for m, idx in enumerate(self.iter_by_order(self.order, self.dim)):
                for d in range(self.dim):
                    values[..., d, m] = self.combine_polyvals_diff(coors, polyvals, d, idx)

        else:
            values = nm.zeros(coors.shape[:-1] + (1, self.n_nod,))  # 1, because no matter the dimension functions have only one value
            polyvals = self.legendreP(coors)
            for m, idx in enumerate(self.iter_by_order(self.order, self.dim)):
                values[..., 0, m] = self.combine_polyvals(coors, polyvals, idx)

        return values

    def iter_by_order(self, order, dim):
        """
        Iterates over all combinations of basis functions indexes
        needed to create multidimensional basis.
        :param order: desired order of multidimensional basis
        :param dim: dimension of the basis
        :yields: tuple containing indexes, use in combine_polyvals and combine_polyvals_der
        :return: None
        """

        # nth(iter(map(lambda x: x + (order - reduce(add,x),)), range(order)), dim)
        # nth(dim, iterate(map(lambda x: x + (order - reduce(add,x),)), map(tuple, range(order))))
        # nth(2, iterate(map(lambda x: x + (order - reduce(add,x),)), map(lambda x: (x,), range(order))))
        porder = order + 1
        if dim == 1:
            for i in range(porder):
                yield((i,))
            return
        elif dim == 2:
            for i in range(porder):
                for j in range(porder - i):
                    m = int(j + porder * i - i / 2 * (i - 1))
                    yield((i, j))
            return
        elif dim == 3:
            for i in range(porder):
                for j in range(porder - i):
                    for k in range(porder - i - j):
                        m = int(1 + ((11 + 12*porder + 3*porder**2)*i)/6 + ((2*porder + 3)*j)/2
                                + k - (2 + porder)*i**2/2 - i*j - j**2/2 + i**3/6)
                        yield((i, j, k))
            return

    def combine_polyvals(self, coors, polyvals, idx):
        raise NotImplementedError("Called abstract method, use some subclass: LegendreTensorPolySpace or LegendreSimplexPolySpace")

    def combine_polyvals_diff(self, coors, polyvals, der, idx):
        raise NotImplementedError("Called abstract method, use some subclass: LegendreTensorPolySpace or LegendreSimplexPolySpace")


    # --------------------- #
    # 1D legendre polyspace #
    # --------------------- #
    funs = [lambda x: x - x + 1,
            lambda x: 2 * x - 1,
            lambda x: (6 * x ** 2 - 6 * x + 1),
            lambda x: (20 * x ** 3 - 30 * x ** 2 + 12 * x - 1),
            lambda x: (70 * x ** 4 - 140 * x ** 3 + 90 * x ** 2 - 20 * x + 1),
            lambda x: (252 * x ** 5 - 630 * x ** 4 + 560 * x ** 3 - 210 * x ** 2 + 30 * x - 1)
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

    def gradlegendreP(self, coors, d=1):
        """
               :param coors: coordinates, preferably in interval [0, 1] for which this basisi is intented
               :return: values in coors of all the legendre polynomials up to self.order
               """
        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)

        values = nm.zeros(sh + (self.order + 1,))

        if self.order == 0:
            return values

        values[..., 1] = 2
        for i in range(2, self.order + 1):
            # values[..., i] = ((2*i + 1) * coors * values[..., i-1] - i * values[..., i-2]) / (i + 1)
            # FIXME transform recursive formula for derivatives too
            values[..., i] = self.get_nth_fun_der(i, d)(coors)
        return values

    def get_nth_fun(self, n):
        """
        Convenience function for testing
        :param n: 0,1 , 2, 3, ...
        :return: n-th function of the legendre basis
        """

        if n < 6:
            return self.funs[n]
        else:
            from scipy.special import comb as comb

            def fun(x):
                val = 0.0
                for k in range(n + 1):
                    val = val + comb(n, k) * comb(n + k, k) * (-x) ** k
                val *= -1 if n % 2 else 1
                return val

            return fun

    def get_nth_fun_der(self, n, l=1):
        """
        Returns lth derivative of nth function
        :param n:
        :param l:
        :return:
        """

        def dfun(x):
            from scipy.special import comb as comb, factorial
            val = x[:] * 0.0
            for k in range(l, n + 1):
                val += comb(n, k) * comb(n + k, k) * factorial(k) / factorial(k - l) * (-x) ** (k - l)
            val *= -1 if (n - l) % 2 else 1
            return val

        return dfun

    # --------------------- #
    # 1D jacobi polyspace #
    # --------------------- #
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


class LegendreTensorProductPolySpace(LegendrePolySpace):
    name = "legendre_tensor_product"

    def __init__(self, name, geometry, order):
        LegendrePolySpace.__init__(self, name, geometry, order)

    def combine_polyvals(self, coors, polyvals, idx):
        return nm.prod(polyvals[..., range(len(idx)), idx], axis=-1)

    def combine_polyvals_diff(self, coors, polyvals, der, idx):
        dimz = range(len(idx))
        derz = nm.zeros(len(idx), dtype=nm.int32)
        derz[der] = 1
        return nm.prod(polyvals[..., dimz, idx, derz], axis=-1)


class LegendreSimplexPolySpace(LegendrePolySpace):
    name = "legendre_simplex"

    def combine_polyvals(self, coors, polyvals, idx):
        if len(idx) == 1:
            return nm.prod(polyvals[..., range(len(idx)), idx], axis=-1)
        else:

            return nm.prod(polyvals[..., range(len(idx)), idx], axis=-1) * (1 - coors[..., 1])**idx[0]

    def combine_polyvals_diff(self, coors, polyvals, der, idx):
        dimz = range(len(idx))
        derz = nm.zeros(len(idx), dtype=nm.int32)
        derz[der] = 1
        return nm.prod(polyvals[..., dimz, idx, derz], axis=-1)


def plot_2Dtensor_basis_grad():
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter


    gel_coors = nm.array([[0, 0], [0, 1], [1, 1], [0, 1]])
    geometry = Struct(n_vertex=4,
                      dim=2,
                      coors=gel_coors)

    order = 3
    bs = LegendreSimplexPolySpace('legb', geometry, order)

    # Make data.
    X = nm.arange(0, 1, 0.025)
    Y = nm.arange(0, 1, 0.025)
    coors = nm.array(nm.meshgrid(X, Y)).T
    Z = bs.eval_base(coors, diff=False)

    Xgrad = nm.linspace(0, 1, 10)
    Ygrad = nm.linspace(0, 1, 10)
    coorsgrad = nm.array(nm.meshgrid(Xgrad, Ygrad)).T
    # Zgrad = bs.eval_base(coorsgrad, diff=True)

    # Zgrad[:,:,:,1:] = Zgrad[:,:,:,1:] /100  # nm.sum(Zgrad[:,:,:,1:]**2, axis=2)[:,:, None, :]

    for i, idx in enumerate(bs.iter_by_order(order, 2)):
        fig = plt.figure("{}>{}".format(i, idx))
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(coors[:, :, 0], coors[:, :, 1], Z[:, :, 0, i], cmap=cm.coolwarm,
                               linewidth=0, antialiased=False, alpha=.6)
        # vec_field = ax.quiver(coorsgrad[:, :, 0], coorsgrad[:, :, 1], -1*nm.ones((10, 10)),
        #                       Zgrad[:, :, 0, i], Zgrad[:, :, 1, i], nm.zeros((10, 10)))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")


    plt.show()


def plot_1D_basis():
    from matplotlib import pylab as plt

    coors = nm.linspace(0, 1)[:, nax]
    geometry = Struct(n_vertex=2,
                      dim=1,
                      coors=coors.copy())

    # bs = CanonicalPolySPace('primb', geometry, 5)
    # vals = bs.eval_base(coors)

    bs = LegendreTensorProductPolySpace('legb', geometry, 3)
    Legvals = bs.eval_base(coors, diff=False)

    # plt.figure("Primitive polyspace")
    # plt.plot(nm.linspace(-1, 1), vals[: ,: ,0])

    plt.figure("Legendre polyspace")
    plt.plot(coors, Legvals[:, 0, :])

    # plt.plot(coors, bs.get_nth_fun(2)(coors))
    # plt.plot(coors, bs.get_nth_fun_der(2)(coors), "--")
    # plt.plot(coors, bs.get_nth_fun_der(2, 2)(coors), "--")
    # plt.plot(coors, bs.get_nth_fun_der(2, 3)(coors), "--")
    # plt.plot([0, 1], [0, 0], 'k')

    plt.show()
    # geometry = Struct(n_vertex = 2,
    #              dim = 1,
    #              coors = self.bbox[:,0:1].copy())


if __name__ == '__main__':
    plot_2Dtensor_basis_grad()
