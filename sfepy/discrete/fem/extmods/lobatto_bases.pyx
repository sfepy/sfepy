# -*- Mode: Python -*-
# cython: language_level=3
"""
Interface to Lobatto bases.
"""
cimport cython

import numpy as np
cimport numpy as np

cimport sfepy.discrete.common.extmods._fmfield as _f
from sfepy.discrete.common.extmods._fmfield cimport FMField

from sfepy.discrete.common.extmods.types cimport int32, float64, complex128

cdef extern from 'lobatto.h':
    int32 _eval_lobatto1d \
          'eval_lobatto1d'(FMField *out, FMField *coors, int32 order)

    int32 _eval_lobatto_tensor_product \
          'eval_lobatto_tensor_product'(FMField *out, FMField *coors,
                                        int32 *nodes,
                                        float64 cmin, float64 cmax,
                                        int32 diff)

    int32 max_order

@cython.boundscheck(False)
def eval_lobatto1d(float64[::1] coors not None,
                   int32 order):
    """
    Evaluate 1D Lobatto functions of the given order in given points.
    """
    cdef int32 n_coor = coors.shape[0]
    cdef np.ndarray[float64, ndim=1] out = np.zeros(n_coor, dtype=np.float64)
    cdef FMField[1] _coors, _out

    _f.array2fmfield1(_coors, coors)
    _f.array2fmfield1(_out, out)

    _eval_lobatto1d(_out, _coors, order)

    return out

@cython.boundscheck(False)
def eval_lobatto_tensor_product(float64[:, ::1]
                                coors not None,
                                int32[:, ::1]
                                nodes not None,
                                float64 cmin, float64 cmax, int32 order,
                                int32 diff=False):
    """
    Evaluate tensor product Lobatto functions of the given order in given
    points.

    Base functions are addressed using the `nodes` array with rows
    corresponding to individual functions and columns to 1D indices (= orders
    when >= 1) into lobatto[] and d_lobatto[] lists for each axis.
    """
    cdef np.ndarray[float64, ndim=3] out
    cdef FMField[1] _coors, _out
    cdef int32 *_nodes = &nodes[0, 0]
    cdef int32 n_coor = coors.shape[0]
    cdef int32 dim = coors.shape[1]
    cdef int32 n_fun = nodes.shape[0]
    cdef int32 nr

    if (order < 1) or (order > max_order):
        raise ValueError('order must be in [1, %d]! (was %d)'
                         % (max_order, order))

    nr = 1 if not diff else dim
    out = np.zeros((n_coor, nr, n_fun), dtype=np.float64)

    _f.array2fmfield2(_coors, coors)
    _f.array2fmfield3(_out, out)

    _eval_lobatto_tensor_product(_out, _coors, _nodes, cmin, cmax, diff)

    return out
