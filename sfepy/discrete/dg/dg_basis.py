import numpy as nm

from numpy import newaxis as nax

from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.base.base import Struct


def iter_by_order(order, dim):
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
            yield ((i,))
        return
    elif dim == 2:
        for i in range(porder):
            for j in range(porder - i):
                m = int(j + porder * i - i / 2 * (i - 1))
                yield ((j, i))
        return
    elif dim == 3:
        for i in range(porder):
            for j in range(porder - i):
                for k in range(porder - i - j):
                    m = int(1 + ((11 + 12 * porder + 3 * porder ** 2) * i) / 6 + ((2 * porder + 3) * j) / 2
                            + k - (2 + porder) * i ** 2 / 2 - i * j - j ** 2 / 2 + i ** 3 / 6)
                    yield ((i, j, k))
        return


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
        coors = 2 * coors - 1
        if diff:
            values = nm.zeros((1,) + coors.shape[:-1] + (self.dim, self.n_nod))  # number of values of derivative is equal to dimension
            polyvals = nm.zeros(coors.shape + (porder, 2))  # 2 is derivative order plus 1, so far we support only first derivative
            polyvals[..., 0] = self.legendreP(coors)
            polyvals[..., 1] = 2 * self.gradlegendreP(coors)  # 2 due to transformation of coordinates

            for m, idx in enumerate(iter_by_order(self.order, self.dim)):
                for d in range(self.dim):
                    values[..., d, m] = self.combine_polyvals_diff(coors, polyvals, d, idx)

        else:
            values = nm.zeros(coors.shape[:-1] + (1, self.n_nod,))  # 1, because no matter the dimension functions have only one value
            polyvals = self.legendreP(coors)
            for m, idx in enumerate(iter_by_order(self.order, self.dim)):
                values[..., 0, m] = self.combine_polyvals(coors, polyvals, idx)

        return values

    def combine_polyvals(self, coors, polyvals, idx):
        raise NotImplementedError("Called abstract method, use some subclass: LegendreTensorPolySpace or LegendreSimplexPolySpace")

    def combine_polyvals_diff(self, coors, polyvals, der, idx):
        raise NotImplementedError("Called abstract method, use some subclass: LegendreTensorPolySpace or LegendreSimplexPolySpace")

    def get_interpol_scheme(self):
        raise NotImplementedError("Called abstract method, use some subclass: LegendreTensorPolySpace or LegendreSimplexPolySpace")


    # --------------------- #
    # 1D legendre polyspace #
    # --------------------- #
    funs = [lambda x: x - x + 1,
            lambda x: 2*x - 1,
            lambda x: (6*x ** 2 - 6*x + 1),
            lambda x: (20*x ** 3 - 30*x ** 2 + 12*x - 1),
            lambda x: (70*x ** 4 - 140*x ** 3 + 90*x ** 2 - 20*x + 1),
            lambda x: (252*x ** 5 - 630*x ** 4 + 560*x ** 3 - 210*x ** 2 + 30*x - 1)
            ]

    def legendreP(self, coors):
        """
        :param coors: coordinates, preferably in interval [-1, 1] for which this basis is intented
        :return: values in coors of all the legendre polynomials up to self.order
        """
        return self.jacobiP(coors, alpha=0, beta=0)

    def gradlegendreP(self, coors, diff=1):
        """
        :param coors: coordinates, preferably in interval [-1, 1] for which this basis is intented
        :return: values in coors of all the legendre polynomials up to self.order
        """
        return self.gradjacobiP(coors, 0, 0, diff=diff)

    def get_nth_fun(self, n):
        """
        Uses shifted legendre
        polynomials formula. Useful for testing.
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

    def get_nth_fun_der(self, n, diff=1):
        """
        Returns diff derivative of nth function. Uses shifted legendre
        polynomials formula. Useful for testing.
        :param n:
        :param diff:
        :return:
        """

        def dfun(x):
            from scipy.special import comb as comb, factorial
            val = x[:] * 0.0
            for k in range(diff, n + 1):
                val += comb(n, k) * comb(n + k, k) * factorial(k) / factorial(k - diff) * (-x) ** (k - diff)
            val *= -1 if (n - diff) % 2 else 1
            return val

        return dfun

    # --------------------- #
    #  1D jacobi polyspace  #
    # --------------------- #
    def jacobiP(self, coors, alpha, beta):
        """
        Values of the jacobi polynomials shifted to interval [-1, 1]
        up to self.order + 1 at coors
        :param coors:
        :param alpha:
        :param beta:
        :return: output shape is shape(coor) + (self.order + 1,)
        """
        from scipy.special import eval_jacobi
        if not isinstance(coors, nm.ndarray):
            sh = (1,)
        else:
            sh = nm.shape(coors)
        values = nm.ones(sh + (self.order + 1,))

        for i in range(self.order + 1):
            values[..., i] = eval_jacobi(i, alpha, beta, coors)
            # for some reason eval_jacobi consumes last dimension if it is one, when called
            # with order array

        return values

    def gradjacobiP(self, coors, alpha, beta, diff=1):
        """
         diff derivative of the jacobi polynomials shifted to interval [-1, 1]
        up to self.order + 1 at coors

        Warning
        Computing values of high-order polynomials (around order > 20) using polynomial coefficients is numerically unstable.
        To evaluate polynomial values, the eval_* functions should be used instead. - but how to get derivative?
        :param coors:
        :param alpha:
        :param beta:
        :param diff:
        :return: output shape is shape(coor) + (self.order + 1,)
        """

        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)

        from scipy.special import jacobi
        values = nm.zeros(sh + (self.order + 1,))
        for i in range(self.order + 1):
            jacob_poly = jacobi(i, alpha, beta)
            jacob_poly.deriv(m=diff)
            values[..., i] = jacob_poly.deriv(m=diff)(coors)

        return values


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

    def get_interpol_scheme(self):
        """
        Builds F and P matrices according to gmsh basis specification:

        Let us assume that the approximation of the view's value over an element is written as a linear
        combination of d basis functions f[i], i=0, ..., d-1 (the coefficients being stored in
        list-of-values). Defining

        f[i] = Sum(j=0, ..., d-1) F[i][j]*p[j],

        with p[j] = u^P[j][0]*v^P[j][1]*w^P[j][2] (u, v and w being the coordinates in the element's parameter
        space), then val-coef-matrix denotes the d x d matrix F and val-exp-matrix denotes the d x 3 matrix P.

        :return:
        """

        from scipy.special import jacobi

        def flattten_upper_left_triag(A):
            res = []
            for i, row in enumerate(A):
                res.append(row[:A.shape[1] - i])
            return nm.concatenate(res)

        P = nm.zeros((self.n_nod, 3))
        for m, idx in enumerate(iter_by_order(self.order, 2)):
            P[m, :self.dim] = idx

        F = nm.zeros((self.n_nod, self.n_nod))
        for m, idx in enumerate(iter_by_order(self.order, 2)):
            xcoefs = list(jacobi(idx[0], 0, 0).coef)[::-1]
            xcoefs = nm.array(xcoefs + [0] * (self.order + 1 - len(xcoefs)))
            ycoefs = list(jacobi(idx[1], 0, 0).coef)[::-1]
            ycoefs = nm.array(ycoefs + [0] * (self.order + 1 - len(ycoefs)))
            coef_mat = nm.outer(ycoefs, xcoefs)
            F[m, :] = flattten_upper_left_triag(coef_mat)

        interp_scheme_struct = Struct(name=self.name + "_interpol_scheme",
                                    F=F, P=P, desc=self.geometry.name)
        return interp_scheme_struct


class LegendreSimplexPolySpace(LegendrePolySpace):
    name = "legendre_simplex"

    def combine_polyvals(self, coors, polyvals, idx):

        from scipy.special import jacobi, eval_jacobi
        if len(idx) == 1:
            nm.prod(polyvals[..., range(len(idx)), idx], axis=-1)
        elif len(idx) == 2:
            r = coors[..., 0]
            s = coors[..., 1]
            a = 2 * (1 + r) / (1 - s) - 1
            b = s
            a[nm.isnan(a)] = 1.

            # a = (a + 1) / 2
            # b = (b + 1) / 2
            # TODO maybe, just maybe somehow transform this to be able to use with jacobi polys on interval [0, 1]
            return nm.sqrt(2) * eval_jacobi(idx[0], 0, 0, a) * eval_jacobi(idx[1], 2*idx[0] + 1, 0, b)*(1 - b)**idx[0]
        elif len(idx) == 3:
            r = coors[..., 0]
            s = coors[..., 1]
            t = coors[..., 2]
            a = -2 * (1 + r) / (s + t) - 1
            b = 2 * (1 + s) / (1 - t) - t
            c = t
            a[nm.isnan(a)] = 1.
            b[nm.isnan(b)] = 1.
            # a = (a + 1) / 2
            # b = (b + 1) / 2
            # TODO maybe, just maybe somehowtransform this to be able to use with jacobi polys on interval [0, 1]
            return nm.sqrt(8) * eval_jacobi(idx[0], 0, 0, a) * \
                   eval_jacobi(idx[1], 2*idx[0] + 1, 0, 0, b) * \
                   eval_jacobi(idx[2], 2*idx[0] + 2*idx[1] + 2, 0, c) * (1 - c)**(idx[0] + idx[1])

    def combine_polyvals_diff(self, coors, polyvals, di, idx):
        dimz = range(len(idx))
        derz = nm.zeros(len(idx), dtype=nm.int32)
        derz[di] = 1
        if len(idx) == 1:
            nm.prod(polyvals[..., range(len(idx)), idx], axis=-1)
        elif len(idx) == 2:
            r = coors[..., 0]
            s = coors[..., 1]
            a = 2 * (1 + r) / (1 - s) - 1
            b = s
            a[nm.isnan(a)] = 1.

            if di == 0:
                return polyvals[..., 0, idx[0], 1] * \
                       self.jacobiP(b, 2*idx[0] + 1, 0)[..., idx[1]]*(1 - b)**idx[0]
            elif di == 1:
                return 2 * polyvals[..., 0, idx[0], 0] * \
                       (self.gradjacobiP(b, 2*idx[0] + 1, 0)[..., idx[1]]*(1 - b)**idx[0] -
                        self.jacobiP(b, 2 * idx[0] + 1, 0,)[..., idx[1]]*(1 - b) ** (idx[0] - 1))
                        # 2 due to transformation of coor in eval basis, TODO refactor!
        elif len(idx) == 3:
            r = 2 * coors[..., 0] - 1
            s = 2 * coors[..., 1] - 1
            t = 2 * coors[..., 2] - 1
            a = -2 * (1 + r) / (s + t) - 1
            b = 2 * (1 + s) / (1 - t) - t
            c = t
            a[nm.isnan(a)] = 1.
            b[nm.isnan(b)] = 1.
            # a = (a + 1) / 2
            # b = (b + 1) / 2
            raise NotImplementedError("Not implemented yet, tough luck :-|")
            # TODO maybe, just maybe somehow transform this to be able to use jacobi polys on interval [0, 1]
            return nm.sqrt(8) * eval_jacobi(idx[0], 0, 0, a) * \
                   eval_jacobi(idx[1], 2 * idx[0] + 1, 0, 0, b) * \
                   eval_jacobi(idx[2], 2 * idx[0] + 2 * idx[1] + 2, 0, c) * (1 - c) ** (idx[0] + idx[1])

    def get_interpol_scheme(self):
        """
        Builds F and P matrices according to gmsh basis specification:

        Let us assume that the approximation of the view's value over an element is written as a linear
        combination of d basis functions f[i], i=0, ..., d-1 (the coefficients being stored in
        list-of-values). Defining

        f[i] = Sum(j=0, ..., d-1) F[i][j]*p[j],

        with p[j] = u^P[j][0]*v^P[j][1]*w^P[j][2] (u, v and w being the coordinates in the element's parameter
        space), then val-coef-matrix denotes the d x d matrix F and val-exp-matrix denotes the d x 3 matrix P.

        :return: struct with name of the scheme, geometry desc and P and F
        """
        F = nm.array([[ 1,   0,  0,  0,  0,  0],
                      [-2,   4,  0,  2,  0,  0],
                      [ 4, -24, 24, -8, 24,  4],
                      [-1,   0,  0,  3,  0,  0],
                      [ 2,  -4,  0,-12, 20, 10],
                      [ 1,   0,  0, -8,  0, 10]])
        P = nm.array([[0, 0, 0],
                      [1, 0, 0],
                      [2, 0, 0],
                      [0, 1, 0],
                      [1, 1, 0],
                      [0, 2, 0]])
        interp_scheme_struct = Struct(name=self.name + "_interpol_scheme",
                                      F=F, P=P, desc=self.geometry.name)
        return interp_scheme_struct


def plot_2Dtensor_basis_grad():
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm

    gel_coors = nm.array([[0, 0], [0, 1], [1, 1], [0, 1]])
    geometry = Struct(n_vertex=4,
                      dim=2,
                      coors=gel_coors)

    order = 1
    bs = LegendreTensorProductPolySpace('legb', geometry, order)

    # Make data.
    X = nm.arange(0, 1, 0.025)
    Y = nm.arange(0, 1, 0.025)
    coors = nm.array(nm.meshgrid(X, Y)).T
    Z = bs.eval_base(coors, diff=False)

    Xgrad = nm.linspace(0, 1, 10)
    Ygrad = nm.linspace(0, 1, 10)
    coorsgrad = nm.array(nm.meshgrid(Xgrad, Ygrad)).T
    Zgrad = bs.eval_base(coorsgrad, diff=True)

    Zgrad[:,:,:,1:] = Zgrad[:,:,:,1:] /100  # nm.sum(Zgrad[:,:,:,1:]**2, axis=2)[:,:, None, :]

    for i, idx in enumerate(iter_by_order(order, 2)):
        fig = plt.figure("{}>{}".format(i, idx))
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(coors[:, :, 0], coors[:, :, 1], Z[:, :, 0, i], cmap=cm.coolwarm,
                               linewidth=0, antialiased=False, alpha=.6)
        vec_field = ax.quiver(coorsgrad[:, :, 0], coorsgrad[:, :, 1], -1*nm.ones((10, 10)),
                               Zgrad[:, :, 0, i], Zgrad[:, :, 1, i], nm.zeros((10, 10)), color='r')
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

    plt.show()


def plot_2Dsimplex_basis_grad():
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

    def simplex_nodes2D(N):
        # function [x,y] = EquiNodes2D(N);
        # Purpose  : Compute (x,y) nodes in equilateral triangle for polynomial of order N

        # total number of nodes
        Np = int((N + 1) * (N + 2) / 2)

        # Create equidistributed nodes on equilateral triangle
        x = nm.zeros((Np,))
        y = nm.zeros((Np,))
        start = 0
        end = 0
        xs = nm.linspace(0, 1, N + 1)
        for n in range(0, N + 1):
            end = end + N - n + 1
            x[start: end] = xs[n]
            y[start: end] = nm.linspace(0, float(N - n) / N, (N - n + 1))
            start = start + N - n + 1
        return x, y, Np

    x, y, Np = simplex_nodes2D(20)

    coors = nm.zeros((Np, 2))
    coors[:, 0] = x
    coors[:, 1] = y
    z = bs.eval_base(coors, diff=False)

    gx, gy, Np = simplex_nodes2D(5)
    coorsgrad = nm.zeros((Np, 2))
    coorsgrad[:, 0] = gx
    coorsgrad[:, 0] = gy
    zgrad = bs.eval_base(coorsgrad, diff=True)

    for i, idx in enumerate(iter_by_order(order, 2)):
        fig = plt.figure("{}>{}".format(i, idx))
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(coors[:, 0], coors[:, 1], z[:, 0, i], linewidth=0.2, antialiased=True)
        vec_field = ax.quiver(coorsgrad[:,  0], coorsgrad[:, 1], -1*nm.ones((Np,)),
                              zgrad[0, :,  0, i], zgrad[0, :,  1, i], nm.zeros((Np,)), color='r')
        #  ax.scatter(x, y, z[:, 0, i])
        ax.scatter(gx, gy, nm.sqrt(zgrad[0, :,  0, i]**2 + zgrad[0, :,  1, i]**2), 'k')
        ax.plot([0, 0, 1, 0], [0, 1, 0, 0], 'k')
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

    order = 4
    bs = LegendreTensorProductPolySpace('legb', geometry, order)
    vals = bs.eval_base(coors)
    grad_vals = bs.eval_base(coors, diff=1)

    plt.figure("Jacobi polyspace")
    plt.plot(coors, vals[:, 0, :])

    plt.plot(coors, bs.get_nth_fun(1)(coors), '.-')
    plt.plot(coors, bs.get_nth_fun(2)(coors), '.-')
    plt.plot(coors, bs.get_nth_fun(3)(coors), '.-')
    plt.plot(coors, bs.get_nth_fun(4)(coors), '.-')
    plt.plot([0, 1], [0, 0], 'k')

    plt.figure("Jacobi polyspace diff")
    plt.plot(coors, grad_vals[0, :,0, :])

    plt.plot(coors, bs.get_nth_fun_der(2)(coors), ".-")
    plt.plot(coors, bs.get_nth_fun_der(3)(coors), ".-")
    plt.plot(coors, bs.get_nth_fun_der(4)(coors), ".-")
    plt.plot([0, 1], [0, 0], 'k')

    plt.show()


if __name__ == '__main__':
    plot_2Dtensor_basis_grad()
    # plot_2Dsimplex_basis_grad()
    # plot_1D_basis()