"""
Isogeometric analysis utilities.

Notes
-----
The function :func:`compute_bezier_extraction_1d()` implements the algorithm
described in [1].

[1] Michael J. Borden, Michael A. Scott, John A. Evans, Thomas J. R. Hughes:
    Isogeometric finite element data structures based on Bezier extraction of
    NURBS, Institute for Computational Engineering and Sciences, The University
    of Texas at Austin, Austin, Texas, March 2010.
"""
import numpy as nm

from sfepy.base.base import assert_

def _get_knots_tuple(knots):
    if isinstance(knots, nm.ndarray) and (knots.ndim == 1):
        knots = (knots,)

    elif not isinstance(knots, tuple):
        raise ValueError('knots must be 1D array or a tuple of 1D arrays!')

    return knots

def get_raveled_index(indices, shape):
    """
    Get a global raveled index corresponding to nD indices into an array of the
    given shape.
    """
    return nm.ravel_multi_index(indices, shape)

def compute_bezier_extraction_1d(knots, degree):
    """
    Compute local (element) Bezier extraction operators for a 1D B-spline
    parametric domain.

    Parameters
    ----------
    knots : array
        The knot vector.
    degree : int
        The curve degree.

    Returns
    -------
    cs : list of 2D arrays
        The element extraction operators.
    """
    knots = nm.asarray(knots, dtype=nm.float64)
    n_knots = knots.shape[0]

    a = degree
    b = a + 1

    # The first element extraction operator.
    cs = [nm.eye(degree + 1, degree + 1, dtype=nm.float64)]
    while (b + 1) < n_knots:
        # The current extraction operator.
        cc = cs[-1]

        # Multiplicity of the knot at location b.
        b0 = b
        while ((b + 1) < n_knots) and (knots[b] == knots[b + 1]):
            b += 1
        mult = b - b0 + 1

        if mult >= degree:
            continue

        # The next extraction operator.
        cn = nm.eye(degree + 1, degree + 1, dtype=nm.float64)
        cs.append(cn)

        alphas = nm.zeros(degree - mult, dtype=nm.float64)
        numer = knots[b] - knots[a]
        for ij in xrange(degree, mult, -1):
            alphas[ij - mult - 1] = numer / (knots[a + ij] - knots[a])

        r = degree - mult
        for ij in xrange(0, r):
            save = r - ij - 1
            s = mult + ij

            for ik in xrange(degree, s, -1):
                alpha = alphas[ik - s - 1]
                cc[:, ik] = alpha * cc[:, ik] + (1.0 - alpha) * cc[:, ik - 1]

            if (b + 1) < n_knots:
                # Update overlapping coefficients for the next operator.
                cn[save : ij + save + 2,
                   save] = cc[degree - ij - 1: degree + 1, degree]

        if (b + 1) < n_knots:
            # The next knot vector interval.
            a = b
            b = b + 1

    return cs

def compute_bezier_extraction(knots, degrees):
    """
    Compute local (element) Bezier extraction operators for a nD B-spline
    parametric domain.

    Parameters
    ----------
    knots : sequence of array or array
        The knot vectors.
    degrees : tuple of ints or int
        Polynomial degrees in each parametric dimension.

    Returns
    -------
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    """
    if isinstance(degrees, int): degrees = [degrees]

    knots = _get_knots_tuple(knots)
    dim = len(knots)
    assert_(dim == len(degrees))

    cs = []
    for ii, knots1d in enumerate(knots):
        cs1d = compute_bezier_extraction_1d(knots1d, degrees[ii])
        cs.append(cs1d)

    return cs

def eval_bernstein_basis(x, degree):
    """
    Evaluate the Bernstein polynomial basis of the given `degree`, and its
    derivatives, in a point `x` in [0, 1].

    Parameters
    ----------
    x : float
        The point in [0, 1].
    degree : int
        The basis degree.

    Returns
    -------
    funs : array
        The `degree + 1` values of the Bernstein polynomial basis.
    ders : array
        The `degree + 1` values of the Bernstein polynomial basis derivatives.
    """
    n_fun = degree + 1

    funs = nm.zeros(n_fun, dtype=nm.float64)
    ders = nm.zeros(n_fun, dtype=nm.float64)

    funs[0] = 1.0

    if degree == 0: return funs, ders

    for ip in xrange(1, n_fun - 1):
        prev = 0.0
        for ifun in xrange(ip + 1):
            tmp = x * funs[ifun]
            funs[ifun] = (1.0 - x) * funs[ifun] + prev
            prev = tmp

    for ifun in xrange(n_fun):
        ders[ifun] = degree * (funs[ifun - 1] - funs[ifun])

    prev = 0.0
    for ifun in xrange(n_fun):
        tmp = x * funs[ifun]
        funs[ifun] = (1.0 - x) * funs[ifun] + prev
        prev = tmp

    return funs, ders
