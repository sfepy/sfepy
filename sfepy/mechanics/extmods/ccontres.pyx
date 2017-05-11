# -*- Mode: Python -*-
cimport cython

import numpy as np
cimport numpy as np

from sfepy.discrete.common.extmods.types cimport int32, uint32, float64

cdef extern from 'contres.h':
     cdef void _getLongestEdgeAndGPs \
     'getLongestEdgeAndGPs'(float64* longestEdge, float64* GPs,
                            int n, int nsd, int ngp, int neq, int nsn,
                            int nes, int nen,
                            uint32* elementID, uint32* segmentID,
                            int32* ISN, int32* IEN, float64* H, float64* X)

def get_longest_edge_and_gps(
    np.ndarray[float64, mode='fortran', ndim=2] GPS not None,
    int neq,
    np.ndarray[uint32, mode='c', ndim=1] elementID not None,
    np.ndarray[uint32, mode='c', ndim=1] segmentID not None,
    np.ndarray[int32, mode='c', ndim=2] ISN not None,
    np.ndarray[int32, mode='c', ndim=2] IEN not None,
    np.ndarray[float64, mode='fortran', ndim=2] H not None,
    np.ndarray[float64, mode='fortran', ndim=2] X not None,
    ):

    cdef float64 tmp[1]
    cdef int n, nsd, ngp, nsn, nes, nen

    n = elementID.shape[0]
    nsd = X.shape[1]
    ngp = H.shape[0]
    nsn = ISN.shape[0]
    nes = ISN.shape[1]
    nen = IEN.shape[1]
    _getLongestEdgeAndGPs(<float64 *> &tmp, &GPS[0, 0],
                          n, nsd, ngp, neq, nsn, nes, nen,
                          &elementID[0], &segmentID[0], &ISN[0, 0], &IEN[0, 0],
                          &H[0, 0], &X[0, 0])

    longestEdge = tmp[0]

    return longestEdge, GPS
