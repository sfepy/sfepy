# -*- Mode: Python -*-
"""
Low level finite element assembling functions.
"""
cimport cython

import numpy as np
cimport numpy as np

ctypedef np.complex128_t complex128
ctypedef np.float64_t float64
ctypedef np.int32_t int32

@cython.boundscheck(False)
def assemble_vector(np.ndarray[float64, mode='c', ndim=1] vec not None,
                    np.ndarray[float64, mode='c', ndim=4] vec_in_els not None,
                    np.ndarray[int32, mode='c', ndim=1] iels not None,
                    float64 sign,
                    np.ndarray[int32, mode='c', ndim=2] conn not None):
    cdef int ii, iel, ir, irg
    cdef int num = iels.shape[0]
    cdef int n_ep = conn.shape[1]
    cdef int32 *pconn
    cdef float64 *val = &vec[0]
    cdef float64 *vec_in_el

    assert num == vec_in_els.shape[0]

    for ii in range(0, num):
        iel = iels[ii]

        pconn = &conn[iel, 0]
        vec_in_el = &vec_in_els[ii, 0, 0, 0]

        for ir in range(0, n_ep):
            irg = pconn[ir]
            if irg < 0: continue

            val[irg] += sign * vec_in_el[ir]

@cython.boundscheck(False)
def assemble_vector_complex(np.ndarray[complex128, mode='c', ndim=1]
                            vec not None,
                            np.ndarray[complex128, mode='c', ndim=4]
                            vec_in_els not None,
                            np.ndarray[int32, mode='c', ndim=1] iels not None,
                            complex128 sign,
                            np.ndarray[int32, mode='c', ndim=2] conn not None):
    cdef int ii, iel, ir, irg
    cdef int num = iels.shape[0]
    cdef int n_ep = conn.shape[1]
    cdef int32 *pconn
    cdef complex128 *val = &vec[0]
    cdef complex128 *vec_in_el

    assert num == vec_in_els.shape[0]

    for ii in range(0, num):
        iel = iels[ii]

        pconn = &conn[iel, 0]
        vec_in_el = &vec_in_els[ii, 0, 0, 0]

        for ir in range(0, n_ep):
            irg = pconn[ir]
            if irg < 0: continue

            val[irg] += sign * vec_in_el[ir]
