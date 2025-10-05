"""
Generate simplex quadrature points. Code taken and adapted from pytools/hedge
by Andreas Kloeckner.
"""
import numpy as nm
from functools import reduce

def generate_decreasing_nonnegative_tuples_summing_to(n, length, min=0,
                                                      max=None):
    if length == 0:
        yield ()
    elif length == 1:
        if n <= max:
            #print "MX", n, max
            yield (n,)
        else:
            return
    else:
        if max is None or n < max:
            max = n

        for i in range(min, max+1):
            #print "SIG", sig, i
            for remainder in generate_decreasing_nonnegative_tuples_summing_to(
                    n-i, length-1, min, i):
                yield (i,) + remainder

def generate_permutations(original):
    """
    Generate all permutations of the list `original'.

    Nicked from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/252178
    """
    if len(original) <=1:
        yield original
    else:
        for perm in generate_permutations(original[1:]):
            for i in range(len(perm)+1):
                #nb str[0:1] works in both string and list contexts
                yield perm[:i] + original[0:1] + perm[i:]

def generate_unique_permutations(original):
    """
    Generate all unique permutations of the list `original'.
    """

    had_those = set()

    for perm in generate_permutations(original):
        if perm not in had_those:
            had_those.add(perm)
            yield perm

def wandering_element(length, wanderer=1, landscape=0):
    for i in range(length):
        yield i*(landscape,) + (wanderer,) + (length-1-i)*(landscape,)

def factorial(n):
    from operator import mul
    assert n == int(n)
    return reduce(mul, (i for i in range(1,n+1)), 1)

def _extended_euclidean(q, r):
    """
    Return a tuple (p, a, b) such that p = aq + br,
    where p is the greatest common divisor.
    """

    # see [Davenport], Appendix, p. 214

    if abs(q) < abs(r):
        p, a, b = _extended_euclidean(r, q)
        return p, b, a

    Q = 1, 0
    R = 0, 1

    while r:
        quot, t = divmod(q, r)
        T = Q[0] - quot*R[0], Q[1] - quot*R[1]
        q, r = r, t
        Q, R = R, T

    return q, Q[0], Q[1]

def _gcd(q, r):
    return _extended_euclidean(q, r)[0]

def _simplify_fraction(a_b):
    (a, b) = a_b
    gcd = _gcd(a,b)
    return (a//gcd, b//gcd)

def get_simplex_cubature(order, dimension):
    r"""
    Cubature on an M{n}-simplex.

    cf.
    A. Grundmann and H.M. Moeller,
    Invariant integration formulas for the n-simplex by combinatorial methods,
    SIAM J. Numer. Anal.  15 (1978), 282--290.

    This cubature rule has both negative and positive weights.  It is exact for
    polynomials up to order :math:`2s+1`, where :math:`s` is given as *order*.
    The integration domain is the unit simplex

    .. math::

        T_n := \{(x_1, \dots, x_n): x_i \ge -1, \sum_i x_i \le -1\}
    """
    s = order
    n = dimension
    d = exact_to = 2*s+1

    points_to_weights = {}

    for i in range(s+1):
        weight = (-1)**i * 2**(-2*s) \
                * (d + n-2*i)**d \
                / factorial(i) \
                / factorial(d+n-i)

        for t in generate_decreasing_nonnegative_tuples_summing_to(s-i, n+1):
            for beta in generate_unique_permutations(t):
                denominator = d+n-2*i
                point = tuple(
                        _simplify_fraction((2*beta_i+1, denominator))
                        for beta_i in beta)

                points_to_weights[point] = points_to_weights.get(point, 0) \
                                           + weight
    from operator import add

    vertices = [-1 * nm.ones((n,))] \
               + [nm.array(x)
                  for x in wandering_element(n, landscape=-1, wanderer=1)]

    pos_points = []
    pos_weights = []
    neg_points = []
    neg_weights = []

    dim_factor = 2**n
    for p, w in points_to_weights.items():
        real_p = reduce(add, (a/b*v for (a,b),v in zip(p, vertices)))
        if w > 0:
            pos_points.append(real_p)
            pos_weights.append(dim_factor*w)
        else:
            neg_points.append(real_p)
            neg_weights.append(dim_factor*w)

    points = nm.array(pos_points + neg_points)
    weights = nm.array(pos_weights + neg_weights)

    return points, weights, exact_to
