# -*- Mode: Python -*-
"""
Template file for 'lobatto.pyx'.
"""
cimport cython

import numpy as np
cimport numpy as np

from types cimport int32, float64, complex128

cdef extern from 'math.h':
    cdef float64 sqrt(float64 x)
    cdef float64 pow(float64 x, float64 y)

ctypedef float64 (*fun)(float64 x)

# Start of generated code.
REPLACE_TEXT
# End of generated code.

@cython.boundscheck(False)
def eval_lobatto(np.ndarray[float64, mode='c', ndim=1] coors not None,
                 int32 order):
    """
    Evaluate Lobatto function of the given order in given points.
    """
    cdef int32 ii
    cdef fun eval_fun
    cdef int32 n_coor = coors.shape[0]
    cdef np.ndarray[float64, ndim=1] out = np.zeros(n_coor, dtype=np.float64)
    cdef float64 *_coors = &coors[0]
    cdef float64 *_out = &out[0]

    if (order < 0) or (order > max_order):
        raise ValueError('order must be in [0, %d]! (was %d)'
                         % (max_order, order))

    eval_fun = lobatto[order]
    for ii in range(0, n_coor):
        _out[ii] = eval_fun(_coors[ii])

    return out

@cython.boundscheck(False)
def eval_lobatto_tensor_product(np.ndarray[float64, mode='c', ndim=2]
                                coors not None,
                                np.ndarray[int32, mode='c', ndim=2]
                                nodes not None,
                                float64 cmin, float64 cmax, int32 order,
                                int32 diff=False):
    cdef np.ndarray[float64, ndim=3] out
    cdef np.ndarray[float64, ndim=2] lambdas
    cdef np.ndarray[float64, ndim=2] xis
    cdef float64 dx = cmax - cmin
    cdef int32 ii, ifun, ic, ir, io, nr
    cdef int32 n_coor = coors.shape[0]
    cdef int32 dim = coors.shape[1]
    cdef int32 n_fun = nodes.shape[0]
    cdef float64 *_xis, *_out
    cdef int32 *_nodes = &nodes[0, 0]
    cdef fun eval_fun

    if (order < 1) or (order > max_order):
        raise ValueError('order must be in [1, %d]! (was %d)'
                         % (max_order, order))

    nr = 1 if not diff else dim
    out = np.ones((n_coor, nr, n_fun), dtype=np.float64)

    # Transform coordinates via affine coordinates lambda to be in [-1, 1].
    lambdas = (coors - cmin) / dx
    xis = 2.0 * lambdas - 1.0
    _xis = &xis[0, 0]
    _out = &out[0, 0, 0]
    if not diff:
        for ii in range(0, dim):
            for ifun in range(0, n_fun):
                eval_fun = lobatto[_nodes[dim * ifun + ii]]
                for ic in range(0, n_coor):
                    _out[n_fun * ic + ifun] *= eval_fun(_xis[dim * ic + ii])

    else:
        for ii in range(0, dim):
            for ifun in range(0, n_fun):
                for ir in range(0, dim):
                    if ir == ii:
                        eval_fun = d_lobatto[_nodes[dim * ifun + ii]]

                    else:
                        eval_fun = lobatto[_nodes[dim * ifun + ii]]

                    for ic in range(0, n_coor):
                        io = n_fun * (nr * ic + ir) + ifun
                        _out[io] *= eval_fun(_xis[dim * ic + ii])

        # Multiply by 2 due to the transformation of coordinates.
        for ii in range(0, n_coor * nr * n_fun):
            _out[ii] *= 2.0

    return out
