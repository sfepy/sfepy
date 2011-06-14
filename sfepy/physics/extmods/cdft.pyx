# -*- Mode: Python -*-
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef extern from "dft.h":
    double vxc(double n, int relat)

cpdef get_vxc(np.ndarray[DTYPE_t, ndim=1] n_qp, int relat=0):
    cdef unsigned int ii
    cdef int num = n_qp.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] out = np.zeros(num, dtype=DTYPE)

    for ii in range(0, num):
        out[ii] = vxc(n_qp[ii], relat)

    return out
