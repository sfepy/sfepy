import numpy as nm

try:
    import sympy as sm
except ImportError:
    sm = None

import pytest

from sfepy.base.base import assert_, ordered_iteritems
import sfepy.base.testing as tst

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

def test_weight_consistency():
    """
    Test if integral of 1 (= sum of weights) gives the domain volume.
    """
    from sfepy.discrete.quadratures import quadrature_tables

    ok = True
    for geometry, qps in ordered_iteritems(quadrature_tables):
        tst.report('geometry:', geometry)
        for order, qp in ordered_iteritems(qps):
            diff = nm.abs(qp.weights.sum() - qp.volume)
            _ok = diff < 1e-14
            tst.report('order %d: %s (%.2e)' % (order, _ok, diff))

            ok = ok and _ok

    assert_(ok)

def _test_quadratures(quad, geometry, order):
    dim = int(geometry[0])
    n_v = int(geometry[2])
    is_simplex = n_v == (dim + 1)

    porder = order if is_simplex else dim * order

    xs, poly, limits, integral = get_poly(porder, dim,
                                          is_simplex=is_simplex)

    tst.report('  polynomial:', poly)
    tst.report('  limits:', limits)
    tst.report('  integral:', integral)

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

    tst.report('  sym. == num.: %f == %f -> %s' %
               (float(integral), val, ok))

    return ok, integral, val

@pytest.mark.slow
def test_quadratures():
    """
    Test if the quadratures have orders they claim to have, using
    symbolic integration by sympy.
    """
    from sfepy.discrete.quadratures import quadrature_tables, tp_geometries
    from sfepy.discrete.integrals import Integral

    if sm is None:
        tst.report('cannot import sympy, skipping')
        return

    quad = Integral('test_integral')

    ok = True
    failed = []
    for geometry, qps in ordered_iteritems(quadrature_tables):
        tst.report('geometry:', geometry)

        if geometry == '0_1':
            continue

        elif geometry in tp_geometries:
            iter_qp = range(1, 11)

        elif geometry == '2_3':
            iter_qp = range(1, 25)

        elif geometry == '3_4':
            iter_qp = range(1, 12)

        elif geometry == '3_6':
            iter_qp = range(4, 9)

        else:
            iter_qp = sorted(qps.keys())

        for order in iter_qp:
            tst.report('order:', order)

            _ok, integral, val = _test_quadratures(quad, geometry, order)

            if not _ok:
                failed.append((geometry, order, float(integral), val))

            ok = ok and _ok

    if not ok:
        tst.report('failed:')
        for aux in failed:
            tst.report(aux, '%.1e' % abs(aux[2] - aux[3]))

    assert_(ok)
