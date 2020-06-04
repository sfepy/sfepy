# -*- coding: utf-8 -*-
from abc import abstractmethod
from functools import reduce
from operator import mul

import numpy as nm

from sfepy.base.ioutils import InDir
from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.base.base import Struct
from scipy.special import jacobi as scp_jacobi
from scipy.special import eval_jacobi as scp_eval_jacobi


def iter_by_order(order, dim, extended=False):
    """Iterates over all combinations of basis functions indexes
    needed to create multidimensional basis in a way that creates hierarchical
    basis

    Parameters
    ----------
    order : int
        desired order of multidimensional basis
    dim : int
        dimension of the basis
    extended : bool
        iterate over extended tensor product basis
        (Default value = False)

    Yields
    ------
    idx : tuple
        containing basis function indexes, used in
        ``_combine_polyvals`` and ``_combine_polyvals_der``
    """

    # nth(iter(map(lambda x: x + (order - reduce(add,x),)), range(order)), dim)
    # nth(dim, iterate(map(lambda x: x + (order - reduce(add,x),)),
    #                  map(tuple, range(order))))
    # nth(2, iterate(map(lambda x: x + (order - reduce(add,x),)),
    #                map(lambda x: (x,), range(order))))
    porder = order + 1
    if dim == 1:
        for i in range(porder):
            yield i,
        return
    elif dim == 2:
        for k in range(porder):
            for i in range(k + 1):
                yield i, k - i
        if not extended: return
        for k in range(1, porder):
            for i in range(1, porder):
                if i + k <= porder - 1:
                    continue
                yield i, k

    elif dim == 3:
        for k in range(porder):
            for j in range(k + 1):
                for i in range(j + 1):
                    yield i, j - i, k - j
        return


def get_n_el_nod(order, dim, extended=False):
    r"""Number of nodes per element for discontinuous legendre basis, i.e.
    number of iterations yielded by iter_by_order

    When extended is False

    .. math::
        N_p =  \frac{(n + 1) \cdot (n + 2) \cdot ... \cdot (n + d)}{d!}
    

    where `n` is the order and `d` the dimension.
    When extended is True

    .. math::
        N_p = (n + 1) ^ d

    where `n` is the order and `d` the dimension.

    Parameters
    ----------
    order : int
        desired order of multidimensional basis
    dim : int
        dimension of the basis
    extended : bool
        iterate over extended tensor product basis
        (Default value = False)


    Returns
    -------
    n_el_nod : int
        number of basis functions in basis
    """
    if extended:
        return (order + 1) ** dim
    return int(reduce(mul, map(lambda i: order + i + 1, range(dim))) /
               reduce(mul, range(1, dim + 1)))


class LegendrePolySpace(PolySpace):
    """Legendre hierarchical polynomials basis, over [0, 1] domain."""

    def __init__(self, name, geometry, order, extended):
        """

        Parameters
        ----------
        name
        geometry
        order : int
            approximation order, 0 for constant functions, 1 for linear etc.
        extended : bool
            for extended tensor product space
        """
        PolySpace.__init__(self, name, geometry, order)
        self.extended = extended  # only tensor product polyspace is extended
        self.n_v = geometry.n_vertex,
        self.dim = geometry.dim
        self.n_nod = get_n_el_nod(self.order, self.dim, self.extended)

        self.coefM = None
        self.expoM = None

    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):
        """Calls _combine_polyvals or _combine_polyvals_diff to build
        multidimensional basis implemented in subclasses to get basis for
        different geometries expects coors to be in shape
            (..., dim),
        Returns values of the basis in shape
            (coors.shape[:-1], 1, n_el_nod)
        or values of the gradient in shape
            (1, coors.shape[:-1], dim, n_el_nod)

        Parameters
        ----------
        coors : array_like
        ori : unused
            we do not need ori, because the basis is discontinous across
            elements this is treated in DGField(Default value = None)
        suppress_errors : has no effect
        diff : int or truthy
             (Default value = 0)
        eps : unused
             (Default value = 1e-15)

        Returns
        -------
        values : ndarray
        """
        coors = 2 * coors - 1  # transofrm from [0, 1] to [-1, 1]
        porder = self.order + 1
        if diff:
            diff = int(diff)
            values = nm.zeros((1,) + coors.shape[:-1] +
                              # (1,) for dummy axis used throughout sfepy
                              (self.dim,) * diff +
                              # (dim,) * diff order is shape of derivation tensor,
                              # we support only first derivative
                              (self.n_nod,))
            polyvals = nm.zeros(coors.shape + (porder,) + (diff + 1,))
            # diff + 1 is number of values of one dimensional base
            polyvals[..., 0] = self.legendreP(coors)
            polyvals[..., 1] = self.gradlegendreP(coors)

            for m, idx in enumerate(
                    iter_by_order(self.order, self.dim, self.extended)):
                for d in range(self.dim):
                    values[..., d, m] = 2 * self._combine_polyvals_diff(coors,
                                                                        polyvals,
                                                                        d, idx)
                    # 2 is due to transformation from [0,1] to [-1,1]
        else:
            values = nm.zeros(coors.shape[:-1] + (1, self.n_nod,))
            # 1, because no matter the dimension functions have only one value
            polyvals = self.legendreP(coors)
            for m, idx in enumerate(
                    iter_by_order(self.order, self.dim, self.extended)):
                values[..., 0, m] = self._combine_polyvals(coors, polyvals, idx)

        return values

    def get_interpol_scheme(self):
        r"""For dim > 1 returns F and P matrices according to gmsh basis
        specification [1]_:
        Let us assume that the approximation of the view's value over an element
        is written as a linear  combination of d basis functions
        :math:`f_i, i=0, ..., n-1` (the coefficients being stored
        in list-of-values).
        
        Defining

        .. math ::
            f_i = \sum\limits_{j=0}^{d-1} F_{ij}\cdot p_j,

        with

        :math:`p_j = u^{P_j^{(0)}}\cdot v^{P_(j)^1}\cdot w^{P_j^{(2)}}`
        (`u`, `v` and `w` being the
        coordinates in the element's parameter space), then val-coef-matrix
        denotes the n x n matrix F and val-exp-matrix denotes the n x 3 matrix P
        where n is number of basis functions as calculated by ``get_n_el_nod``.
        
        Expects matrices to be saved in atributes coefM and expoM!

        .. [1] Remacle, J.-F., Chevaugeon, N., Marchandise, E., & Geuzaine, C.
           (2007). Efficient visualization of high-order finite elements.
           International Journal for Numerical Methods in Engineering, 69(4),
           750-771. https://doi.org/10.1002/nme.1787

        Returns
        -------
        interp_scheme_struct : Struct
            Struct with name of the scheme, geometry desc and P and F
        """
        interp_scheme_struct = None
        if self.dim > 1:
            interp_scheme_struct = Struct(name=self.name + "_interpol_scheme",
                                          F=self.coefM, P=self.expoM,
                                          desc=self.geometry.name)
        return interp_scheme_struct

    @abstractmethod
    def _combine_polyvals(self, coors, polyvals, idx):
        """Combines Legendre or Jacobi polynomials to get muldidimensional
        basis values according to element topolgy

        Parameters
        ----------
        coors : array_like
            coordinates of evaluation points
        polyvals : array_like
            values of legendre polynomials precomputed in _eval_base
        idx : tuple
            function index, for tensor-product correspond to orders of polynomials in variables

        Returns
        -------
        values : ndarray
        """

    @abstractmethod
    def _combine_polyvals_diff(self, coors, polyvals, der, idx):
        """Combines Legendre or Jacobi polynomials to get muldidimensional
        basis derivative values according to element topolgy.

        Parameters
        ----------
        coors : array_like
            coordinates of evaluation points
        polyvals : array_like
            values of legendre polynomials precomputed in _eval_base
        der : int
            derivative variable
        idx : tuple
            function index, for tensor-product correspond to orders of polynomials in variables

        Returns
        -------
        values : ndarray
        """

    # --------------------- #
    # 1D legendre polyspace #
    # --------------------- #
    def legendreP(self, coors):
        """
        Parameters
        ----------
        coors : array_like
            coordinates, preferably in interval [-1, 1] for which this basis is
            intented

        Returns
        -------
        values : ndarray
            values at coors of all the legendre polynomials up to self.order

        """
        return self.jacobiP(coors, alpha=0, beta=0)

    def gradlegendreP(self, coors, diff=1):
        """
        Parameters
        ----------
        diff : int
            default 1
        coors : array_like
            coordinates, preferably in interval [-1, 1] for which this basis is
            intented

        Returns
        -------
        values : ndarray
            values at coors of all the legendre polynomials up to self.order

        """
        return self.gradjacobiP(coors, 0, 0, diff=diff)

    """
    Explict legendre polynomials up to order 5
    """
    legendre_funs = [lambda x: 0 * x + 1,
                     # we need constant preserving shape and type of x
                     lambda x: 2 * x - 1,
                     lambda x: (6 * x ** 2 - 6 * x + 1),
                     lambda x: (20 * x ** 3 - 30 * x ** 2 + 12 * x - 1),
                     lambda x: (70 * x ** 4 - 140 * x ** 3
                                + 90 * x ** 2 - 20 * x + 1),
                     lambda x: (252 * x ** 5 - 630 * x ** 4 + 560 * x ** 3
                                - 210 * x ** 2 + 30 * x - 1)
                     ]

    def get_nth_fun(self, n):
        """Uses shifted Legendre polynomials formula on interval [0, 1].

        Convenience function for testing

        Parameters
        ----------
        n : int

        Returns
        -------
        fun : callable
            n-th function of the legendre basis
        """

        if n < 6:
            return self.legendre_funs[n]
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
        """Returns diff derivative of nth function. Uses shifted legendre
        polynomials formula on interval [0, 1].

        Useful for testing.

        Parameters
        ----------
        n : int
        diff : int
             (Default value = 1)

        Returns
        -------
        fun : callable
            derivative of n-th function of the 1D legendre basis

        """

        def dfun(x):
            """

            Parameters
            ----------
            x :
                

            Returns
            -------

            """
            from scipy.special import comb as comb, factorial
            val = x[:] * 0.0
            for k in range(diff, n + 1):
                val += comb(n, k) * comb(n + k, k) * factorial(k) / \
                       factorial(k - diff) * (-x) ** (k - diff)
            val *= -1 if (n - diff) % 2 else 1
            return val

        return dfun
    # --------------------- #

    # --------------------- #
    #  1D jacobi polyspace  #
    # --------------------- #
    def jacobiP(self, coors, alpha, beta):
        """Values of the jacobi polynomials on interval [-1, 1]
        up to self.order + 1 at coors

        Parameters
        ----------
        coors : array_like
        beta : float
        alpha : float
            

        Returns
        -------
        values : ndarray
            output shape is shape(coor) + (self.order + 1,)

        """
        if not isinstance(coors, nm.ndarray):
            sh = (1,)
        else:
            sh = nm.shape(coors)
        values = nm.ones(sh + (self.order + 1,))

        for i in range(self.order + 1):
            values[..., i] = scp_eval_jacobi(i, alpha, beta, coors)
            # for some reason eval_jacobi consumes last dimension if it is one,
            # when called with order array

        return values

    def gradjacobiP(self, coors, alpha, beta, diff=1):
        """diff derivative of the jacobi polynomials on interval [-1, 1]
        up to self.order + 1 at coors

        Parameters
        ----------
        coors :

        alpha : float

        beta : float

        diff : int
             (Default value = 1)

        Returns
        -------
        values : ndarray
            output shape is shape(coor) + (self.order + 1,)

        """
        if isinstance(coors, (int, float)):
            sh = (1,)
        else:
            sh = nm.shape(coors)

        values = nm.zeros(sh + (self.order + 1,))
        for i in range(self.order + 1):
            jacob_poly = scp_jacobi(i, alpha, beta)
            values[..., i] = jacob_poly.deriv(m=diff)(coors)
            # Warning
            # Computing values of high-order polynomials (around order > 20)
            # using polynomial coefficients is numerically unstable. To evaluate
            # polynomial values, the eval_* functions should be used instead.
            # On that note jacob_poly.deriv seems to use stable version.

        return values
    # --------------------- #


class LegendreTensorProductPolySpace(LegendrePolySpace):
    name = "legendre_tensor_product"

    def __init__(self, name, geometry, order, ):
        super().__init__(name, geometry, order, extended=True)
        self.n_nod = get_n_el_nod(self.order, self.dim, self.extended)
        if self.dim > 1:
            self.coefM, self.expoM = self.build_interpol_scheme()

    def _combine_polyvals(self, coors, polyvals, idx):
        return nm.prod(polyvals[..., range(len(idx)), idx], axis=-1)

    def _combine_polyvals_diff(self, coors, polyvals, dvar, idx):
        dimz = range(len(idx))
        derz = nm.zeros(len(idx), dtype=nm.int32)
        derz[dvar] = 1
        return nm.prod(polyvals[..., dimz, idx, derz], axis=-1)

    def build_interpol_scheme(self):
        """Builds F and P matrices returned by self.get_interpol_scheme.

        Note that this function returns coeficients according to gmsh
        parametrization of Quadrangle i.e. [-1, 1] x [-1, 1] and hence the form
        of basis function is not the same as exhibited by the
        LegendreTensorProductPolySpace object which acts on parametrization
        [0, 1] x [0, 1].


        Returns
        -------
        F : ndarray
            coefficient matrix
        P : ndarray
            exponent matrix
        """
        P = nm.zeros((self.n_nod, 3), dtype=nm.int32)
        for m, idx in enumerate(
                iter_by_order(self.order, self.dim, self.extended)):
            P[m, :self.dim] = idx

        F = nm.zeros((self.n_nod, self.n_nod))
        for m, idx in enumerate(
                iter_by_order(self.order, self.dim, self.extended)):
            xcoefs = list(scp_jacobi(idx[0], 0, 0).coef)[::-1]
            xcoefs = nm.array(xcoefs + [0] * (self.order + 1 - len(xcoefs)))
            ycoefs = list(scp_jacobi(idx[1], 0, 0).coef)[::-1]
            ycoefs = nm.array(ycoefs + [0] * (self.order + 1 - len(ycoefs)))
            coef_mat = nm.outer(xcoefs, ycoefs)
            F[m, :] = [coef_mat[idx] for idx in
                       iter_by_order(self.order, self.dim, self.extended)]
        return F, P


class LegendreSimplexPolySpace(LegendrePolySpace):
    name = "legendre_simplex"

    def __init__(self, name, geometry, order, extended=False):
        super().__init__(name, geometry, order, extended)
        if self.dim > 1:
            indir = InDir(__file__)
            try:
                self.coefM = nm.loadtxt(
                        indir("legendre2D_simplex_coefs.txt")
                )[:self.n_nod, :self.n_nod]
                self.expoM = nm.loadtxt(
                        indir("legendre2D_simplex_expos.txt")
                )[:self.n_nod, :]
            except IOError as e:
                raise IOError(
                    ("File {} not found, run gen_legendre_simplex_base.py"
                     + " to generate it.")
                        .format(e.args[0]))

    def _combine_polyvals(self, coors, polyvals, idx):

        from scipy.special import eval_jacobi
        if len(idx) == 1:  # 1D
            nm.prod(polyvals[..., range(len(idx)), idx], axis=-1)
        elif len(idx) == 2:  # 2D
            r = coors[..., 0]
            s = coors[..., 1]
            a = 2 * (1 + r) / (1 - s) - 1
            a[nm.isnan(a)] = -1.
            b = s
            return eval_jacobi(idx[0], 0, 0, a) \
                   * eval_jacobi(idx[1],2 * idx[0] + 1, 0,b) * (1 - b) ** idx[0]
        elif len(idx) == 3:  # 3D
            r = coors[..., 0]
            s = coors[..., 1]
            t = coors[..., 2]
            a = -2 * (1 + r) / (s + t) - 1
            b = 2 * (1 + s) / (1 - t) - t
            c = t
            a[nm.isnan(a)] = -1.
            b[nm.isnan(b)] = -1.
            return eval_jacobi(idx[0], 0, 0, a) * \
                   eval_jacobi(idx[1], 2 * idx[0] + 1, 0, 0, b) * \
                   eval_jacobi(idx[2], 2 * idx[0] + 2 * idx[1] + 2, 0, c) * \
                                                    (1 - c) ** (idx[0] + idx[1])

    def _combine_polyvals_diff(self, coors, polyvals, dvar, idx):

        if len(idx) == 1:  # 1D
            dimz = range(len(idx))
            derz = nm.zeros(len(idx), dtype=nm.int32)
            derz[dvar] = 1
            return nm.prod(polyvals[..., dimz, idx, derz], axis=-1)
        elif len(idx) == 2:  # 2D
            r = coors[..., 0]
            s = coors[..., 1]
            a = 2 * (1 + r) / (1 - s) - 1
            b = s
            a[nm.isnan(a)] = -1.
            di = idx[0]
            dj = idx[1]

            fa = self.jacobiP(a, 0, 0)[..., di]
            dfa = self.gradjacobiP(a, 0, 0)[..., di]
            gb = self.jacobiP(b, 2 * di + 1, 0)[..., dj]
            dgb = self.gradjacobiP(b, 2 * di + 1, 0)[..., dj]

            if dvar == 0:  # d/dr
                # r - derivative
                # d / dr = da / dr * d / da + db / dr * d / db
                #        = (2 / (1 - s)) d / da = (2 / (1 - b)) d / da
                dmodedr = dfa * gb
                dmodedr = dmodedr * (
                            (0.5 * (1 - b)) ** (di - 1)) if di > 0 else dmodedr
                return 2 ** di * dmodedr

            elif dvar == 1:  # d/ds
                # s - derivative
                # d / ds = ((1 + a) / 2) / ((1 - b) / 2) d / da + d / db
                dmodeds = dfa * (gb * (0.5 * (1 + a)))
                dmodeds = dmodeds * (
                            (0.5 * (1 - b)) ** (di - 1)) if di > 0 else dmodeds

                tmp = dgb * ((0.5 * (1 - b)) ** di)
                tmp = tmp - 0.5 * di * gb * (
                            (0.5 * (1 - b)) ** (di - 1)) if di > 0 else tmp
                dmodeds = dmodeds + fa * tmp
                return 2 ** di * dmodeds
        elif len(idx) == 3:  # 3D
            #  UNTESTED!
            r = coors[..., 0]
            s = coors[..., 1]
            t = coors[..., 2]
            a = -2 * (1 + r) / (s + t) - 1
            b = 2 * (1 + s) / (1 - t) - t
            c = t
            a[nm.isnan(a)] = -1.
            b[nm.isnan(b)] = -1.
            di = idx[0]
            dj = idx[1]
            dk = idx[2]
            fa = self.jacobiP(a, 0, 0)[..., di]
            dfa = self.gradjacobiP(a, 0, 0)[..., di]
            gb = self.jacobiP(b, 2 * di + 1, 0)[..., dj]
            dgb = self.gradjacobiP(b, 2 * di + 1, 0)[..., dj]
            hc = self.jacobiP(c, 2 * (di + dj) + 2, 0, dk)
            dhc = self.gradjacobiP(c, 2 * (di + dj) + 2, 0, dk)

            V3Dr = dfa * (gb * hc)
            V3Dr = V3Dr * ((0.5 * (1 - b)) ** (di - 1)) if di > 0 else V3Dr
            V3Dr = V3Dr * ((0.5 * (1 - c)) ** (
                        di + dj - 1)) if di + dj > 0 else V3Dr

            tmp = dgb * ((0.5 * (1 - b)) ** di)
            tmp = tmp + (-0.5 * di) * (
                        gb * (0.5 * (1 - b)) ** (di - 1)) if di > 0 else tmp
            tmp = tmp * ((0.5 * (1 - c)) ** (
                        di + dj - 1)) if di + dj > 0 else tmp
            tmp = fa * (tmp * hc)

            if dvar == 0:
                return V3Dr * (2 ** (2 * di + dj + 1))
            elif dvar == 1:
                V3Ds = 0.5 * (1 + a) * V3Dr

                V3Ds = V3Ds + tmp
                return V3Ds * (2 ** (2 * di + dj + 1))
            elif dvar == 2:
                V3Dt = 0.5 * (1 + a) * V3Dr + 0.5 * (1 + b) * tmp
                tmp = dhc * ((0.5 * (1 - c)) ** (di + dj))
                tmp = tmp - 0.5 * (di + dj) * (hc * ((0.5 * (1 - c)) ** (
                            di + dj - 1))) if di + dj > 0 else tmp
                tmp = fa * (gb * tmp)
                tmp = tmp * ((0.5 * (1 - b)) ** di)
                V3Dt = V3Dt + tmp
                return V3Dt * (2 ** (2 * di + dj + 1))
