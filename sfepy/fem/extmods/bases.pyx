# -*- Mode: Python -*-
"""
Polynomial base functions and related utilities.
"""
cimport cython

import numpy as np
cimport numpy as np

ctypedef np.complex128_t complex128
ctypedef np.float64_t float64
ctypedef np.int32_t int32

@cython.boundscheck(False)
cpdef get_barycentric_coors(np.ndarray[float64, ndim=2]
                            coors,
                            np.ndarray[float64, ndim=2]
                            mtx_i,
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
    cdef int error
    cdef int ir, ic, ii
    cdef int n_coor = coors.shape[0]
    cdef int dim = coors.shape[1]
    cdef int n_v = dim + 1
    cdef float64 val
    cdef np.ndarray[float64, ndim=2] bc = np.zeros((n_coor, n_v),
                                                   dtype=np.float64)

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

    return bc
