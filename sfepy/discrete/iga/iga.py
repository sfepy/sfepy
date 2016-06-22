"""
Isogeometric analysis utilities.

Notes
-----
The functions :func:`compute_bezier_extraction_1d()` and
:func:`eval_nurbs_basis_tp()` implement the algorithms described in [1].

[1] Michael J. Borden, Michael A. Scott, John A. Evans, Thomas J. R. Hughes:
    Isogeometric finite element data structures based on Bezier extraction of
    NURBS, Institute for Computational Engineering and Sciences, The University
    of Texas at Austin, Austin, Texas, March 2010.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_
from six.moves import range

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

def get_unraveled_indices(index, shape):
    """
    Get nD indices into an array of the given shape corresponding to a global
    raveled index.
    """
    return nm.unravel_index(index, shape)

def tensor_product(a, b):
    """
    Compute tensor product of two 2D arrays with possibly different shapes. The
    result has the form::

       c = [[a00 b, a01 b, ...],
            [a10 b, a11 b, ...],
             ...
             ...               ]
    """
    c = nm.empty((a.shape[0] * b.shape[0],
                  a.shape[1] * b.shape[1]), dtype=b.dtype)

    n0 = b.shape[0]
    n1 = b.shape[1]
    for ir in range(a.shape[0]):
        for ic in range(a.shape[1]):
            c[n1 * ir : n1 * (ir + 1),
              n0 * ic : n0 * (ic + 1)] = a[ir, ic] * b

    return c

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
    cs : array of 2D arrays (3D array)
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

        # The next extraction operator.
        if (b + 1) < n_knots:
            cn = nm.eye(degree + 1, degree + 1, dtype=nm.float64)
            cs.append(cn)

        if mult < degree:
            alphas = nm.zeros(degree - mult, dtype=nm.float64)
            numer = knots[b] - knots[a]
            for ij in range(degree, mult, -1):
                alphas[ij - mult - 1] = numer / (knots[a + ij] - knots[a])

            r = degree - mult
            for ij in range(0, r):
                save = r - ij - 1
                s = mult + ij

                for ik in range(degree, s, -1):
                    alpha = alphas[ik - s - 1]
                    cc[:, ik] = (alpha * cc[:, ik]
                                 + (1.0 - alpha) * cc[:, ik - 1])

                if (b + 1) < n_knots:
                    # Update overlapping coefficients for the next operator.
                    cn[save : ij + save + 2,
                       save] = cc[degree - ij - 1: degree + 1, degree]

        if (b + 1) < n_knots:
            # The next knot vector interval.
            a = b
            b = b + 1

    return nm.asarray(cs, dtype=nm.float64)

def compute_bezier_extraction(knots, degrees):
    """
    Compute local (element) Bezier extraction operators for a nD B-spline
    parametric domain.

    Parameters
    ----------
    knots : sequence of array or array
        The knot vectors.
    degrees : sequence of ints or int
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

def combine_bezier_extraction(cs):
    """
    For a nD B-spline parametric domain, combine the 1D element extraction
    operators in each parametric dimension into a single operator for each nD
    element.

    Parameters
    ----------
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.

    Returns
    -------
    ccs : list of 2D arrays
        The combined element extraction operators.
    """
    dim = len(cs)

    if dim == 3:
        c0, c1, c2 = cs[0], cs[1], cs[2]

        ncc = (len(c0), len(c1), len(c2))
        ccs = [None] * nm.prod(ncc)
        for i0 in range(len(c0)):
            for i1 in range(len(c1)):
                for i2 in range(len(c2)):
                    cc = tensor_product(c0[i0], tensor_product(c1[i1], c2[i2]))
                    ii = get_raveled_index([i0, i1, i2], ncc)
                    ccs[ii] = cc

    elif dim == 2:
        c0, c1 = cs[0], cs[1]

        ncc = (len(c0), len(c1))
        ccs = [None] * nm.prod(ncc)
        for i0 in range(len(c0)):
            for i1 in range(len(c1)):
                cc = tensor_product(c0[i0], c1[i1])
                ii = get_raveled_index([i0, i1], ncc)
                ccs[ii] = cc

    else:
        ccs = cs[0]

    return ccs

def create_connectivity_1d(n_el, knots, degree):
    """
    Create connectivity arrays of 1D Bezier elements.

    Parameters
    ----------
    n_el : int
        The number of elements.
    knots : array
        The knot vector.
    degree : int
        The basis degree.

    Returns
    -------
    conn : array
        The connectivity of the global NURBS basis.
    bconn : array
        The connectivity of the Bezier basis.
    """
    # Get multiplicities of NURBS knots.
    n_knots = len(knots)
    mul = [0]
    ii = degree + 1
    while ii < (n_knots - degree - 1):
        i0 = ii
        while (ii < (n_knots - degree - 2)) and (knots[ii] == knots[ii + 1]):
            ii += 1
        mul.append(ii - i0 + 1)
        ii += 1

    mul = nm.array(mul)[:, None]

    aux1 = nm.arange(degree + 1)[None, :]
    conn = aux1 + nm.cumsum(mul, 0)

    # Bezier basis knots have multiplicity equal to degree.
    aux2 = nm.arange(n_el)[:, None]
    bconn = aux1 + degree * aux2

    return conn.astype(nm.int32), bconn.astype(nm.int32)

def create_connectivity(n_els, knots, degrees):
    """
    Create connectivity arrays of nD Bezier elements.

    Parameters
    ----------
    n_els : sequence of ints
        The number of elements in each parametric dimension.
    knots : sequence of array or array
        The knot vectors.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.

    Returns
    -------
    conn : array
        The connectivity of the global NURBS basis.
    bconn : array
        The connectivity of the Bezier basis.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    knots = _get_knots_tuple(knots)
    dim = len(n_els)
    assert_(dim == len(degrees) == len(knots))

    conns = []
    bconns = []
    n_gfuns = []
    n_gbfuns = []
    for ii, n_el in enumerate(n_els):
        conn1d, bconn1d = create_connectivity_1d(n_el, knots[ii], degrees[ii])
        conns.append(conn1d)
        bconns.append(bconn1d)

        n_gfuns.append(conn1d.max() + 1)
        n_gbfuns.append(bconn1d.max() + 1)

    n_el = nm.prod(n_els)
    n_efuns = degrees + 1
    n_efun = nm.prod(n_efuns)

    if dim == 3:
        def make_conn_3d(conns, n_gfuns):
            conn = nm.empty((n_el, n_efun), dtype=nm.int32)
            for ie0 in range(n_els[0]):
                c0 = conns[0][ie0]
                for ie1 in range(n_els[1]):
                    c1 = conns[1][ie1]
                    for ie2 in range(n_els[2]):
                        c2 = conns[2][ie2]
                        ie = get_raveled_index([ie0, ie1, ie2], n_els)

                        for il0 in range(n_efuns[0]):
                            cl0 = c0[il0]
                            for il1 in range(n_efuns[1]):
                                cl1 = c1[il1]
                                for il2 in range(n_efuns[2]):
                                    cl2 = c2[il2]

                                    iloc = get_raveled_index([il0, il1, il2],
                                                             n_efuns)
                                    ig = get_raveled_index([cl0, cl1, cl2],
                                                           n_gfuns)
                                    conn[ie, iloc] = ig

            return conn

        conn = make_conn_3d(conns, n_gfuns)
        bconn = make_conn_3d(bconns, n_gbfuns)

    elif dim == 2:
        def make_conn_2d(conns, n_gfuns):
            conn = nm.empty((n_el, n_efun), dtype=nm.int32)
            for ie0 in range(n_els[0]):
                c0 = conns[0][ie0]
                for ie1 in range(n_els[1]):
                    c1 = conns[1][ie1]
                    ie = get_raveled_index([ie0, ie1], n_els)

                    for il0 in range(n_efuns[0]):
                        cl0 = c0[il0]
                        for il1 in range(n_efuns[1]):
                            cl1 = c1[il1]

                            iloc = get_raveled_index([il0, il1], n_efuns)
                            ig = get_raveled_index([cl0, cl1], n_gfuns)
                            conn[ie, iloc] = ig

            return conn

        conn = make_conn_2d(conns, n_gfuns)
        bconn = make_conn_2d(bconns, n_gbfuns)

    else:
        conn = conns[0]
        bconn = bconns[0]

    return conn, bconn

def compute_bezier_control(control_points, weights, ccs, conn, bconn):
    """
    Compute the control points and weights of the Bezier mesh.

    Parameters
    ----------
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    ccs : list of 2D arrays
        The combined element extraction operators.
    conn : array
        The connectivity of the global NURBS basis.
    bconn : array
        The connectivity of the Bezier basis.

    Returns
    -------
    bezier_control_points : array
        The control points of the Bezier mesh.
    bezier_weights : array
        The weights of the Bezier mesh.
    """
    n_bpoints = bconn.max() + 1
    dim = control_points.shape[1]

    bezier_control_points = nm.zeros((n_bpoints, dim), dtype=nm.float64)
    bezier_weights = nm.zeros(n_bpoints, dtype=nm.float64)

    for ie, ec in enumerate(conn):
        cc = ccs[ie]
        bec = bconn[ie]

        ew = weights[ec]
        ecp = control_points[ec]

        bew = nm.dot(cc.T, ew)

        becp = (1.0 / bew[:, None]) * nm.dot(cc.T, ew[:, None] * ecp)
        bezier_control_points[bec] = becp
        bezier_weights[bec] = bew

    return bezier_control_points, bezier_weights

def get_bezier_topology(bconn, degrees):
    """
    Get a topology connectivity corresponding to the Bezier mesh connectivity.

    In the referenced Bezier control points the Bezier mesh is interpolatory.

    Parameters
    ----------
    bconn : array
        The connectivity of the Bezier basis.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.

    Returns
    -------
    tconn : array
        The topology connectivity (corner nodes, or vertices, of Bezier
        elements) with vertex ordering suitable for a FE mesh.
    """
    shape = nm.asarray(degrees) + 1
    dim = len(shape)

    ii = nm.arange(bconn.shape[1]).reshape(shape)

    if dim == 3:
        corners = [ii[0, 0, 0], ii[-1, 0, 0], ii[-1, -1, 0], ii[0, -1, 0],
                   ii[0, 0, -1], ii[-1, 0, -1], ii[-1, -1, -1], ii[0, -1, -1]]

    elif dim == 2:
        corners = [ii[0, 0], ii[-1, 0], ii[-1, -1], ii[0, -1]]

    else:
        corners = [ii[0], ii[-1]]

    tconn = bconn[:, corners]

    return tconn

def get_patch_box_regions(n_els, degrees):
    """
    Get box regions of Bezier topological mesh in terms of element corner
    vertices of Bezier mesh.

    Parameters
    ----------
    n_els : sequence of ints
        The number of elements in each parametric dimension.
    degrees : sequence of ints or int
        Polynomial degrees in each parametric dimension.

    Returns
    -------
    regions : dict
        The Bezier mesh vertices of box regions.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    n_els = nm.asarray(n_els)
    dim = len(n_els)

    shape = n_els * degrees + 1

    regions = {}
    if dim == 3:
        aux0 = nm.arange(0, shape[2], degrees[2], dtype=nm.uint32)
        aux1 = nm.arange(0, shape[2] * shape[1], shape[2] * degrees[1],
                         dtype=nm.uint32)
        aux2 = nm.arange(0, shape[2] * shape[1] * shape[0],
                         shape[2] * shape[1] * degrees[0], dtype=nm.uint32)

        aux01 = (aux0[None, :] + aux1[:, None]).ravel()
        aux02 = (aux0[None, :] + aux2[:, None]).ravel()
        aux12 = (aux1[None, :] + aux2[:, None]).ravel()

        regions.update({
            'xi00' : aux01,
            'xi01' : aux01 + shape[2] * shape[1] * (shape[0] - 1),
            'xi10' : aux02,
            'xi11' : aux02 + shape[2] * (shape[1] - 1),
            'xi20' : aux12,
            'xi21' : aux12 + shape[2] - 1,
        })

    elif dim == 2:
        aux0 = nm.arange(0, shape[1], degrees[1], dtype=nm.uint32)
        aux1 = nm.arange(0, shape[1] * shape[0], shape[1] * degrees[0],
                         dtype=nm.uint32)

        regions.update({
            'xi00' : aux0,
            'xi01' : aux0 + shape[1] * (shape[0] - 1),
            'xi10' : aux1,
            'xi11' : aux1 + shape[1] - 1,
        })

    else:
        regions.update({
            'xi00' : nm.array([0], dtype=nm.uint32),
            'xi01' : nm.array([shape[0] - 1], dtype=nm.uint32),
        })

    return regions

def get_facet_axes(dim):
    """
    For each reference Bezier element facet return the facet axes followed by
    the remaining (perpendicular) axis, as well as the remaining axis
    coordinate of the facet.

    Parameters
    ----------
    dim : int
        The topological dimension.

    Returns
    -------
    axes : array
        The axes of the reference element facets.
    coors : array
        The remaining coordinate of the reference element facets.
    """
    if dim == 3:
        axes = [[1, 0, 2], [2, 1, 0], [0, 2, 1],
                [0, 1, 2], [1, 2, 0], [2, 0, 1]]
        coors = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    elif dim == 2:
        axes = [[0, 1], [1, 0], [0, 1], [1, 0]]
        coors = [0.0, 1.0, 1.0, 0.0]

    else:
        axes = [[0]]
        coors = None

    return nm.array(axes, dtype=nm.uint32), nm.array(coors, dtype=nm.float64)

def get_surface_degrees(degrees):
    """
    Get degrees of the NURBS patch surfaces.

    Parameters
    ----------
    degrees : sequence of ints or int
        Polynomial degrees in each parametric dimension.

    Returns
    -------
    sdegrees : list of arrays
        The degrees of the patch surfaces, in the order of the reference Bezier
        element facets.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    dim = len(degrees)

    if dim == 3:
        sdegrees = [(degrees[0], degrees[1]),
                    (degrees[1], degrees[2]),
                    (degrees[0], degrees[2]),
                    (degrees[0], degrees[1]),
                    (degrees[1], degrees[2]),
                    (degrees[0], degrees[2])]
        sdegrees = nm.array(sdegrees, dtype=nm.uint32)

    elif dim == 2:
        sdegrees = degrees[[0, 1, 0, 1]]

    else:
        sdegrees = None

    return sdegrees

def create_boundary_qp(coors, dim):
    """
    Create boundary quadrature points from the surface quadrature points.

    Uses the Bezier element tensor product structure.

    Parameters
    ----------
    coors : array, shape (n_qp, d)
        The coordinates of the surface quadrature points.
    dim : int
        The topological dimension.

    Returns
    -------
    bcoors : array, shape (n_qp, d + 1)
        The coordinates of the boundary quadrature points.
    """
    # Boundary QP - use tensor product structure.
    axes, acoors = get_facet_axes(dim)
    n_f = len(axes)

    bcoors = nm.empty((n_f, coors.shape[0], coors.shape[1] + 1),
                      dtype=nm.float64)
    ii = nm.arange(bcoors.shape[1], dtype=nm.uint32)
    for ik in range(n_f):
        for ic in range(bcoors.shape[2] - 1):
            bcoors[ik, :, axes[ik, ic]] = coors[:, ic]
        bcoors[ik, ii, axes[ik, -1]] = acoors[ik]

    return bcoors

def get_bezier_element_entities(degrees):
    """
    Get faces and edges of a Bezier mesh element in terms of indices into the
    element's connectivity (reference Bezier element entities).

    Parameters
    ----------
    degrees : sequence of ints or int
        Polynomial degrees in each parametric dimension.

    Returns
    -------
    faces : list of arrays
        The indices for each face or None if not 3D.
    edges : list of arrays
        The indices for each edge or None if not at least 2D.
    vertices : list of arrays
        The indices for each vertex.

    Notes
    -----
    The ordering of faces and edges has to be the same as in
    :data:`sfepy.discrete.fem.geometry_element.geometry_data`.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    dim = len(degrees)
    shape = degrees + 1

    n_dof = nm.prod(shape)

    aux = nm.arange(n_dof, dtype=nm.uint32).reshape(shape)
    if dim == 3:
        faces = [aux[:, :, 0],
                 aux[0, :, :],
                 aux[:, 0, :],
                 aux[:, :, -1],
                 aux[-1, :, :],
                 aux[:, -1, :]]
        faces = [ii.ravel() for ii in faces]
        edges = [aux[:, 0, 0],
                 aux[-1, :, 0],
                 aux[:, -1, 0],
                 aux[0, :, 0],
                 aux[:, 0, -1],
                 aux[-1, :, -1],
                 aux[:, -1, -1],
                 aux[0, :, -1],
                 aux[0, 0, :],
                 aux[0, -1, :],
                 aux[-1, -1, :],
                 aux[-1, 0, :]]
        vertices = [aux[0, 0, 0],
                    aux[-1, 0, 0],
                    aux[-1, -1, 0],
                    aux[0, -1, 0],
                    aux[0, 0, -1],
                    aux[-1, 0, -1],
                    aux[-1, -1, -1],
                    aux[0, -1, -1]]
        vertices = [ii[None] for ii in vertices]

    elif dim == 2:
        faces = None
        edges = [aux[:, 0],
                 aux[-1, :],
                 aux[:, -1],
                 aux[0, :]]
        vertices = [aux[0, 0],
                    aux[-1, 0],
                    aux[-1, -1],
                    aux[0, -1]]
        vertices = [ii[None] for ii in vertices]

    else:
        faces, edges = None, None
        vertices = [aux[:1], aux[-1:]]

    return faces, edges, vertices

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

    for ip in range(1, n_fun - 1):
        prev = 0.0
        for ifun in range(ip + 1):
            tmp = x * funs[ifun]
            funs[ifun] = (1.0 - x) * funs[ifun] + prev
            prev = tmp

    for ifun in range(n_fun):
        ders[ifun] = degree * (funs[ifun - 1] - funs[ifun])

    prev = 0.0
    for ifun in range(n_fun):
        tmp = x * funs[ifun]
        funs[ifun] = (1.0 - x) * funs[ifun] + prev
        prev = tmp

    return funs, ders

def eval_nurbs_basis_tp(qp, ie, control_points, weights, degrees, cs, conn):
    """
    Evaluate the tensor-product NURBS shape functions in a quadrature point for
    a given Bezier element.

    Parameters
    ----------
    qp : array
        The quadrature point coordinates with components in [0, 1] reference
        element domain.
    ie : int
        The Bezier element index.
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    conn : array
        The connectivity of the global NURBS basis.

    Returns
    -------
    R : array
        The NURBS shape functions.
    dR_dx : array
        The NURBS shape functions derivatives w.r.t. the physical coordinates.
    det : array
        The Jacobian of the mapping to the unit reference element.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    dim = len(degrees)
    assert_(dim == len(qp) == len(cs))

    n_efuns = degrees + 1
    n_efun = nm.prod(n_efuns)
    n_efuns_max = n_efuns.max()

    assert_(n_efun == conn.shape[1])

    # Element connectivity.
    ec = conn[ie]

    # Element control points and weights.
    W = weights[ec]
    P = control_points[ec]

    # 1D Bernstein basis B, dB/dxi.
    B = nm.empty((dim, n_efuns_max), dtype=nm.float64)
    dB_dxi = nm.empty((dim, n_efuns_max), dtype=nm.float64)
    for ii in range(dim):
        (B[ii, :n_efuns[ii]],
         dB_dxi[ii, :n_efuns[ii]]) = eval_bernstein_basis(qp[ii], degrees[ii])

    # 1D B-spline basis N = CB, dN/dxi = C dB/dxi.
    N = nm.empty((dim, n_efuns_max), dtype=nm.float64)
    dN_dxi = nm.empty((dim, n_efuns_max), dtype=nm.float64)
    n_els = [len(ii) for ii in cs]
    ic = get_unraveled_indices(ie, n_els)
    for ii in range(dim):
        C = cs[ii][ic[ii]]

        N[ii, :n_efuns[ii]] = nm.dot(C, B[ii, :n_efuns[ii]])
        dN_dxi[ii, :n_efuns[ii]] = nm.dot(C, dB_dxi[ii, :n_efuns[ii]])

    # Numerators and denominator for tensor-product NURBS basis R, dR/dxi.
    R = nm.empty(n_efun, dtype=nm.float64)
    dR_dxi = nm.empty((n_efun, dim), dtype=nm.float64)

    w = 0 # w_b
    dw_dxi = nm.zeros(dim, dtype=nm.float64) # dw_b/dxi
    a = 0 # Basis function index.
    if dim == 3:
        for i0 in range(n_efuns[0]):
            for i1 in range(n_efuns[1]):
                for i2 in range(n_efuns[2]):
                    R[a] = N[0, i0] * N[1, i1] * N[2, i2] * W[a]
                    w += R[a]

                    dR_dxi[a, 0] = dN_dxi[0, i0] * N[1, i1] * N[2, i2] * W[a]
                    dw_dxi[0] += dR_dxi[a, 0]

                    dR_dxi[a, 1] = N[0, i0] * dN_dxi[1, i1] * N[2, i2] * W[a]
                    dw_dxi[1] += dR_dxi[a, 1]

                    dR_dxi[a, 2] = N[0, i0] * N[1, i1] * dN_dxi[2, i2] * W[a]
                    dw_dxi[2] += dR_dxi[a, 2]

                    a += 1

    elif dim == 2:
        for i0 in range(n_efuns[0]):
            for i1 in range(n_efuns[1]):
                R[a] = N[0, i0] * N[1, i1] * W[a]
                w += R[a]

                dR_dxi[a, 0] = dN_dxi[0, i0] * N[1, i1] * W[a]
                dw_dxi[0] += dR_dxi[a, 0]

                dR_dxi[a, 1] = N[0, i0] * dN_dxi[1, i1] * W[a]
                dw_dxi[1] += dR_dxi[a, 1]

                a += 1

    else:
        for i0 in range(n_efuns[0]):
            R[a] = N[0, i0] * W[a]
            w += R[a]

            dR_dxi[a, 0] = dN_dxi[0, i0] * W[a]
            dw_dxi[0] += dR_dxi[a, 0]

            a += 1

    # Finish R <- R / w_b.
    R /= w

    # Finish dR/dxi. D == W C dB/dxi, dR/dxi = (D - R dw_b/dxi) / w_b.
    dR_dxi = (dR_dxi - R[:, None] * dw_dxi) / w

    # Mapping reference -> physical domain dxi/dx.
    # x = sum P_a R_a, dx/dxi = sum P_a dR_a/dxi, invert.
    dx_dxi = nm.dot(P.T, dR_dxi)
    det = nm.linalg.det(dx_dxi)

    dxi_dx = nm.linalg.inv(dx_dxi)

    # dR/dx.
    dR_dx = nm.dot(dR_dxi, dxi_dx)

    return R, dR_dx, det

def eval_mapping_data_in_qp(qps, control_points, weights, degrees, cs, conn,
                            cells=None):
    """
    Evaluate data required for the isogeometric domain reference mapping in the
    given quadrature points. The quadrature points are the same for all Bezier
    elements and should correspond to the Bernstein basis degree.

    Parameters
    ----------
    qps : array
        The quadrature points coordinates with components in [0, 1] reference
        element domain.
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    conn : array
        The connectivity of the global NURBS basis.
    cells : array, optional
        If given, use only the given Bezier elements.

    Returns
    -------
    bfs : array
        The NURBS shape functions in the physical quadrature points of all
        elements.
    bfgs : array
        The NURBS shape functions derivatives w.r.t. the physical coordinates
        in the physical quadrature points of all elements.
    dets : array
        The Jacobians of the mapping to the unit reference element in the
        physical quadrature points of all elements.
    """
    if cells is None:
        cells = nm.arange(conn.shape[0])

    n_el = len(cells)
    n_qp = qps.shape[0]
    dim = control_points.shape[1]
    n_efuns = degrees + 1
    n_efun = nm.prod(n_efuns)

    # Output Jacobians.
    dets = nm.empty((n_el, n_qp, 1, 1), dtype=nm.float64)

    # Output shape functions.
    bfs = nm.empty((n_el, n_qp, 1, n_efun), dtype=nm.float64)

    # Output gradients of shape functions.
    bfgs = nm.empty((n_el, n_qp, dim, n_efun), dtype=nm.float64)

    # Loop over elements.
    for iseq, ie in enumerate(cells):
        # Loop over quadrature points.
        for iqp, qp in enumerate(qps):
            bf, bfg, det = eval_nurbs_basis_tp(qp, ie,
                                               control_points, weights,
                                               degrees, cs, conn)
            bfs[iseq, iqp] = bf
            bfgs[iseq, iqp] = bfg.T
            dets[iseq, iqp] = det

    return bfs, bfgs, dets

def eval_variable_in_qp(variable, qps,
                        control_points, weights, degrees, cs, conn,
                        cells=None):
    """
    Evaluate a field variable in the given quadrature points. The quadrature
    points are the same for all Bezier elements and should correspond to the
    Bernstein basis degree. The field variable is defined by its DOFs - the
    coefficients of the NURBS basis.

    Parameters
    ----------
    variable : array
        The DOF values of the variable with n_c components, shape (:, n_c).
    qps : array
        The quadrature points coordinates with components in [0, 1] reference
        element domain.
    control_points : array
        The NURBS control points.
    weights : array
        The NURBS weights.
    degrees : sequence of ints or int
        The basis degrees in each parametric dimension.
    cs : list of lists of 2D arrays
        The element extraction operators in each parametric dimension.
    conn : array
        The connectivity of the global NURBS basis.
    cells : array, optional
        If given, use only the given Bezier elements.

    Returns
    -------
    coors : array
        The physical coordinates of the quadrature points of all elements.
    vals : array
        The field variable values in the physical quadrature points.
    dets : array
        The Jacobians of the mapping to the unit reference element in the
        physical quadrature points.
    """
    if cells is None:
        cells = nm.arange(conn.shape[0])

    n_el = len(cells)
    n_qp = qps.shape[0]
    dim = control_points.shape[1]
    nc = variable.shape[1]

    # Output values of the variable.
    vals = nm.empty((n_el * n_qp, nc), dtype=nm.float64)

    # Output physical coordinates of QPs.
    coors = nm.empty((n_el * n_qp, dim), dtype=nm.float64)

    # Output Jacobians.
    dets = nm.empty((n_el * n_qp, 1), dtype=nm.float64)

    # Loop over elements.
    for iseq, ie in enumerate(cells):
        ec = conn[ie]
        vals_e = variable[ec]
        cps_e = control_points[ec]

        # Loop over quadrature points.
        for iqp, qp in enumerate(qps):
            ii = n_qp * iseq + iqp
            bf, bfg, det = eval_nurbs_basis_tp(qp, ie,
                                               control_points, weights,
                                               degrees, cs, conn)
            vals_qp = nm.dot(bf, vals_e)
            vals[ii, :] = vals_qp

            coors_qp = nm.dot(bf, cps_e)
            coors[ii, :] = coors_qp

            dets[ii] = det

    return coors, vals, dets
