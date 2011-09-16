# -*- Mode: Python -*-
cimport cython

import numpy as np
cimport numpy as np

ctypedef np.float64_t float64

cdef extern from "dft.h":
    double vxc(double n, int mode, int relat)

@cython.boundscheck(False)
def get_vxc(np.ndarray[float64, mode='c', ndim=4] n_qp not None,
            int mode=0, int relat=0):
    cdef unsigned int ii
    cdef int num = n_qp.shape[0]
    cdef int nq = n_qp.shape[1]
    cdef np.ndarray[float64, ndim=4] out = np.zeros((num, nq, 1, 1),
                                                    dtype=np.float64)
    cdef float64 *pout = &out[0, 0, 0, 0]
    cdef float64 *pn_qp = &n_qp[0, 0, 0, 0]

    assert (n_qp.shape[2] == 1) and (n_qp.shape[3] == 1)

    for ii in range(0, num * nq):
        pout[ii] = vxc(pn_qp[ii], mode, relat)

    return out
