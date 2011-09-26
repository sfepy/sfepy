# -*- Mode: Python -*-
"""
Polynomial base functions and related utilities.
"""
cimport cython

cimport numpy as np
import numpy as np

from types cimport int32, float64, complex128

@cython.boundscheck(False)
cdef _get_barycentric_coors(np.ndarray[float64, ndim=2] bc,
                            np.ndarray[float64, ndim=2] coors,
                            np.ndarray[float64, ndim=2] mtx_i,
                            float64 eps=1e-8,
                            int check_errors=False):
    """
    Get barycentric (area in 2D, volume in 3D) coordinates of points.

    Parameters
    ----------
    bc : array
        The barycentric coordinates, shape `(n_coor, dim + 1)`. Then
        reference element coordinates `xi = dot(bc, ref_coors)`.
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.
    """
    cdef int error
    cdef int ir, ic, ii
    cdef int n_coor = coors.shape[0]
    cdef int n_v = mtx_i.shape[0]
    cdef int dim = n_v - 1
    cdef float64 val

    for ir in range(0, n_coor):
        for ic in range(0, n_v):
            val = 0.0;
            for ii in range(0, dim):
                val += mtx_i[ic, ii] * coors[ir, ii]
            val += mtx_i[ic, dim]

            error = False
            if val < 0.0:
                if val > -eps:
                    val = 0.0

                else:
                    error = True

            if val > 1.0:
                if val < (1.0 + eps):
                    val = 1.0

                else:
                    error = True

            if check_errors and error:
                msg = 'point %d outside of element! (%s)' % (ic, coors[ic])
                raise ValueError(msg)

            bc[ir, ic] = val

@cython.boundscheck(False)
def get_barycentric_coors(np.ndarray[float64, ndim=2] coors not None,
                          np.ndarray[float64, ndim=2] mtx_i not None,
                          float64 eps=1e-8,
                          int check_errors=False):
    """
    Get barycentric (area in 2D, volume in 3D) coordinates of points.

    Parameters
    ----------
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.

    Returns
    -------
    bc : array
        The barycentric coordinates, shape `(n_coor, dim + 1)`. Then
        reference element coordinates `xi = dot(bc, ref_coors)`.
    """
    cdef int n_coor = coors.shape[0]
    cdef int n_v = mtx_i.shape[0]
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, n_v),
                                                   dtype=np.float64)

    _get_barycentric_coors(bc, coors, mtx_i, eps, check_errors)
    return bc

@cython.boundscheck(False)
@cython.cdivision(True)
cdef _eval_lagrange_simplex(np.ndarray[float64, ndim=3] out,
                            np.ndarray[float64, ndim=2] bc,
                            np.ndarray[float64, ndim=2] mtx_i,
                            np.ndarray[int32, ndim=2] nodes,
                            int order, int diff=False):
    """
    Evaluate Lagrange base polynomials in given points on simplex domain.
    """
    cdef int n_coor = bc.shape[0]
    cdef int n_v = bc.shape[1]
    cdef int dim = n_v - 1
    cdef int n_nod = nodes.shape[0]
    cdef int ii, ir, ic, i1, i2, inod, n_i1, n_ii
    cdef float64 dval, dd, vv
    cdef float64 *pout

    assert n_coor == out.shape[0]
    assert n_v == nodes.shape[1]

    if not diff:
        for ic in range(0, n_coor):
            pout = &out[ic, 0, 0]

            for inod in range(0, n_nod):
                pout[inod] = 1.0
                for i1 in range(0, n_v):
                    n_i1 = nodes[inod, i1]
                    for i2 in range(0, n_i1):
                        pout[inod] *= (order * bc[ic, i1] - i2) / (i2 + 1.0)

    else:
        for ic in range(0, n_coor):
            pout = &out[ic, 0, 0]

            for inod in range(0, n_nod):
                pout[inod] = 0.0

                for ii in range(0, n_v):
                    vv = 1.0

                    for i1 in range(0, n_v):
                        if i1 == ii: continue
                        n_i1 = nodes[inod, i1]

                        for i2 in range(0, n_i1):
                            vv *= (order * bc[ic, i1] - i2) / (i2 + 1.0)

                    dval = 0.0
                    n_ii = nodes[inod, ii]
                    for i1 in range(0, n_ii):
                        dd = 1.0

                        for i2 in range(0, n_ii):
                            if i1 == i2: continue

                            dd *= (order * bc[ic, ii] - i2) / (i2 + 1.0)

                        dval += dd * order / (i1 + 1.0)

                    for ir in range(0, dim):
                        pout[n_nod*ir+inod] += vv * dval * mtx_i[ii, ir]

@cython.boundscheck(False)
def eval_lagrange_simplex(np.ndarray[float64, ndim=2] coors not None,
                          np.ndarray[float64, ndim=2] mtx_i not None,
                          np.ndarray[int32, ndim=2] nodes not None,
                          int order, int diff=False,
                          float64 eps=1e-15,
                          int check_errors=True):
    """
    Evaluate Lagrange base polynomials in given points on simplex domain.

    Parameters
    ----------
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of simplex coordinates matrix, shape `(dim + 1, dim + 1)`.
    nodes : array
        The description of finite element nodes, shape `(n_nod, dim + 1)`.
    order : int
        The polynomial order.
    diff : bool
        If True, return base function derivatives.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.

    Returns
    -------
    out : array
        The evaluated base functions, shape `(n_coor, 1 or dim, n_nod)`.
    """
    cdef int bdim
    cdef int n_coor = coors.shape[0]
    cdef int dim = mtx_i.shape[0] - 1
    cdef int n_nod = nodes.shape[0]
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, dim + 1),
                                                   dtype=np.float64)

    if diff:
        bdim = dim

    else:
        bdim = 1

    cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_nod),
                                                    dtype=np.float64)

    _get_barycentric_coors(bc, coors, mtx_i, eps, check_errors)
    _eval_lagrange_simplex(out, bc, mtx_i, nodes, order, diff)

    return out

@cython.boundscheck(False)
cdef _eval_lagrange_tensor_product(np.ndarray[float64, ndim=3] out,
                                   np.ndarray[float64, ndim=3] bc,
                                   np.ndarray[float64, ndim=2] mtx_i,
                                   np.ndarray[float64, ndim=3] base1d,
                                   np.ndarray[int32, ndim=2] nodes,
                                   int order, int diff=False):
    """
    Evaluate Lagrange base polynomials in given points on tensor product
    domain.
    """
    cdef int ii, idim, im, ic
    cdef int n_coor = out.shape[0]
    cdef int nr = out.shape[1]
    cdef int n_nod = out.shape[2]
    cdef int dim = bc.shape[0]
    cdef int out_size = n_coor * nr * n_nod
    cdef float64 *val, *val1d
    cdef np.ndarray[int32, ndim=2] pnodes
    cdef np.ndarray[float64, ndim=2] bc1

    val = &out[0, 0, 0]
    for im in range(0, out_size):
        val[im] = 1.0

    val1d = &base1d[0, 0, 0]

    if not diff:
        for ii in range(0, dim):
            pnodes = nodes[:, 2 * ii : 2 * ii + 2]
            bc1 = bc[ii, :, :]

            _eval_lagrange_simplex(base1d, bc1, mtx_i, pnodes, order, diff)

            for im in range(0, n_coor):
                for ic in range(0, n_nod):
                    out[im, 0, ic] *= base1d[im, 0, ic]

    else:
        for ii in range(0, dim):
            pnodes = nodes[:, 2 * ii : 2 * ii + 2]
            bc1 = bc[ii, :, :]

            for idim in range(0, dim):
                if ii == idim:
                    _eval_lagrange_simplex(base1d, bc1, mtx_i, pnodes, order,
                                           diff)

                else:
                    _eval_lagrange_simplex(base1d, bc1, mtx_i, pnodes, order,
                                           False)

                for im in range(0, n_coor):
                    for ic in range(0, n_nod):
                        out[im, idim, ic] *= base1d[im, 0, ic]

@cython.boundscheck(False)
def eval_lagrange_tensor_product(np.ndarray[float64, ndim=2] coors not None,
                                 np.ndarray[float64, ndim=2] mtx_i not None,
                                 np.ndarray[int32, ndim=2] nodes not None,
                                 int order, int diff=False,
                                 float64 eps=1e-15,
                                 int check_errors=True):
    """
    Evaluate Lagrange base polynomials in given points on tensor product
    domain.

    Parameters
    ----------
    coors : array
        The coordinates of the points, shape `(n_coor, dim)`.
    mtx_i : array
        The inverse of 1D simplex coordinates matrix, shape `(2, 2)`.
    nodes : array
        The description of finite element nodes, shape `(n_nod, 2 * dim)`.
    order : int
        The polynomial order.
    diff : bool
        If True, return base function derivatives.
    eps : float
        The tolerance for snapping out-of-simplex point back to the simplex.
    check_errors : bool
        If True, raise ValueError if a barycentric coordinate is outside
        the snap interval `[-eps, 1 + eps]`.

    Returns
    -------
    out : array
        The evaluated base functions, shape `(n_coor, 1 or dim, n_nod)`.
    """
    cdef int ii, idim, im, ic
    cdef int n_coor = coors.shape[0]
    cdef int n_nod = nodes.shape[0]
    cdef int dim = coors.shape[1]
    cdef np.ndarray[float64, ndim=3] bc = np.zeros((dim, n_coor, 2),
                                                   dtype=np.float64)
    cdef np.ndarray[float64, ndim=3] base1d = np.zeros((n_coor, 1, n_nod),
                                                       dtype=np.float64)
    if diff:
        bdim = dim

    else:
        bdim = 1

    cdef np.ndarray[float64, ndim=3] out = np.zeros((n_coor, bdim, n_nod),
                                                    dtype=np.float64)

    for ii in range(0, dim):
        _get_barycentric_coors(bc[ii, :, :], coors[:, ii : ii + 1],
                               mtx_i, eps, check_errors)

    _eval_lagrange_tensor_product(out, bc, mtx_i, base1d,
                                  nodes, order, diff)

    return out
