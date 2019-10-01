from __future__ import absolute_import
import numpy as nm
from six.moves import range

try:
    import sympy as sm
except ImportError:
    sm = None

from sfepy.base.testing import TestCommon
from sfepy.base.base import assert_, ordered_iteritems

def symarray(prefix, shape):
    """
    Copied from SymPy so that the tests pass for its different versions.
    """
    arr = nm.empty(shape, dtype=object)
    for index in nm.ndindex(shape):
        arr[index] = sm.Symbol('%s_%s' % (prefix, '_'.join(map(str, index))))
    return arr

def get_poly(order, dim, is_simplex=False):
    """
    Construct a polynomial of given `order` in space dimension `dim`,
    and integrate it symbolically over a rectangular or simplex domain
    for coordinates in [0, 1].
    """
    xs = symarray('x', dim)

    opd = max(1, int((order + 1) / dim))

    poly = 1.0
    oo = 0
    for ii, x in enumerate(xs):
        if ((oo + opd) > order) or (ii == (len(xs) - 1)):
            opd = max(order - oo, 0)

        poly *= (x**opd + 1)
        oo += opd

    assert_(oo == order)

    limits = [[xs[ii], 0, 1] for ii in range(dim)]
    if is_simplex:
        for ii in range(1, dim):
            for ip in range(0, ii):
                limits[ii][2] -= xs[ip]

    integral = sm.integrate(poly, *reversed(limits))

    return xs, poly, limits, integral

class Test(TestCommon):

    @staticmethod
    def from_conf( conf, options ):

        test = Test(conf=conf, options=options)
        return test

    def test_weight_consistency(self):
        """
        Test if integral of 1 (= sum of weights) gives the domain volume.
        """
        from sfepy.discrete.quadratures import quadrature_tables

        ok = True
        for geometry, qps in ordered_iteritems(quadrature_tables):
            self.report('geometry:', geometry)
            for order, qp in ordered_iteritems(qps):
                diff = nm.abs(qp.weights.sum() - qp.volume)
                _ok = diff < 1e-14
                self.report('order %d: %s (%.2e)' % (order, _ok, diff))

                ok = ok and _ok

        return ok

    def test_quadratures(self):
        """
        Test if the quadratures have orders they claim to have, using
        symbolic integration by sympy.
        """
        from sfepy.discrete.quadratures import quadrature_tables, tp_geometries
        from sfepy.discrete.integrals import Integral

        if sm is None:
            self.report('cannot import sympy, skipping')
            return True

        quad = Integral('test_integral')

        ok = True
        failed = []
        for geometry, qps in ordered_iteritems(quadrature_tables):
            self.report('geometry:', geometry)

            if geometry == '0_1':
                continue

            elif geometry in tp_geometries:
                iter_qp = range(1, 11)

            elif geometry == '2_3':
                iter_qp = range(1, 25)

            elif geometry == '3_4':
                iter_qp = range(1, 12)

            else:
                iter_qp = sorted(qps.keys())

            for order in iter_qp:
                self.report('order:', order)

                aux = self._test_quadratures(quad, geometry, order)
                _ok, integral, val = aux

                if not _ok:
                    failed.append((geometry, order, float(integral), val))

                ok = ok and _ok

        if not ok:
            self.report('failed:')
            for aux in failed:
                self.report(aux, '%.1e' % abs(aux[2] - aux[3]))

        return ok

    def _test_quadratures(self, quad, geometry, order):
        dim = int(geometry[0])
        n_v = int(geometry[2])
        is_simplex = n_v == (dim + 1)

        porder = order if is_simplex else dim * order

        xs, poly, limits, integral = get_poly(porder, dim,
                                              is_simplex=is_simplex)

        self.report('  polynomial:', poly)
        self.report('  limits:', limits)
        self.report('  integral:', integral)

        def fun(coors):
            vals = nm.empty(coors.shape[0], dtype=nm.float64)

            subs = {}
            for ir, cc in enumerate(coors):
                for ic, x in enumerate(xs):
                    subs[x] = coors[ir,ic]
                vals[ir] = float(poly.subs(subs))

            return vals

        val = quad.integrate(fun, order=order, geometry=geometry)
        ok = nm.allclose(val, float(integral), rtol=0.0, atol=1e-11)

        self.report('  sym. == num.: %f == %f -> %s' %
                    (float(integral), val, ok))

        return ok, integral, val
