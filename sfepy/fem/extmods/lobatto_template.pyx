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
