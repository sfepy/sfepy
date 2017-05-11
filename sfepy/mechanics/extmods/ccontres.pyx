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
     cdef void _getAABB \
     'getAABB'(float64* AABBmin, float64* AABBmax, int nsd, int nnod,
               float64* X, float64 longestEdge, int32* IEN, int32* ISN,
               uint32* elementID, uint32* segmentID,
               int n, int nsn, int nes, int nen, int neq)

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

def get_AABB(
    np.ndarray[float64, mode='fortran', ndim=2] X not None,
    int longestEdge,
    np.ndarray[int32, mode='c', ndim=2] IEN not None,
    np.ndarray[int32, mode='c', ndim=2] ISN not None,
    np.ndarray[uint32, mode='c', ndim=1] elementID not None,
    np.ndarray[uint32, mode='c', ndim=1] segmentID not None,
    int neq,
    ):

    cdef np.ndarray[float64, mode='c', ndim=1] AABBmin
    cdef np.ndarray[float64, mode='c', ndim=1] AABBmax
    cdef int n, nnod, nsd, nsn, nes, nen

    n = elementID.shape[0]
    nnod = X.shape[0]
    nsd = X.shape[1]
    nsn = ISN.shape[0]
    nes = ISN.shape[1]
    nen = IEN.shape[1]

    AABBmin = np.empty((nsd,), dtype=np.float64)
    AABBmax = np.empty((nsd,), dtype=np.float64)

    _getAABB(&AABBmin[0], &AABBmax[0], nsd, nnod,
             &X[0, 0], longestEdge, &IEN[0, 0], &ISN[0, 0],
             &elementID[0], &segmentID[0],
             n, nsn, nes, nen, neq)

    return AABBmin, AABBmax

def init_global_search(
    np.ndarray[int32, mode='c', ndim=1] N not None,
    np.ndarray[float64, mode='c', ndim=1] AABBmin not None,
    np.ndarray[float64, mode='c', ndim=1] AABBmax not None,
    np.ndarray[float64, mode='fortran', ndim=2] X not None,
    ):
    """
    The linked list initialization. The head array contains, at the position
    Ic, the index of the first point that belongs to the cell Ic, the second
    point index is then next[head[Ic]], the third point index is
    next[next[head[Ic]]] etc. - the next array points from the i-th point in
    each cell to the (i+1)-th point, until -1 is reached.
    """
    cdef int i, lower, upper, Ic
    cdef np.ndarray[int32, mode='c', ndim=1] I
    cdef np.ndarray j
    cdef np.ndarray[int32, mode='c', ndim=1] head
    cdef np.ndarray[int32, mode='c', ndim=1] next

    dim = X.shape[1]
    nnod = X.shape[0]

    head = np.empty((np.prod(N),), dtype=np.int32);
    head[:] = -1

    next = np.empty((nnod,), dtype=np.int32);
    next[:] = -1

    # Loop over points to find a cell each point is in.
    for i in range(nnod):
        lower = np.sum(X[i] < AABBmin)
        higher = np.sum(X[i] > AABBmax)

        if lower + higher: continue

        I = np.floor(N * (X[i] - AABBmin) / (AABBmax-AABBmin)).astype(np.int32)

        j = I == N
        I[j] = N[j] - 1

        Ic = np.ravel_multi_index(I, N)

        next[i] = head[Ic]
        head[Ic] = i

    return head, next
