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

def compute_bezier_extraction_1d(knots, degree):
    """
    Compute local (element) Bezier extraction operators for a 1D B-spline
    parameteric domain.

    Parameters
    ----------
    knots : array
        The knot vector
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
