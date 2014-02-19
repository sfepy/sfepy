# -*- Mode: Python -*-
cimport cython

import numpy as np
cimport numpy as np

from sfepy.discrete.fem.extmods.types cimport int32

cdef extern from 'common.h':
    cdef void errclear()

cdef extern from 'rcm.h':
    cdef void rcm_genrcm(int32 *perm, int32 neqns, int32 *xadj, int32 n_ptr,
                         int32 *adjncy, int32 n_indx)

    cdef int32 gr_permuteInPlace(int32 *row, int32 n_row,
                                 int32 *col, int32 n_col,
                                 int32 *perm, int32 n_perm,
                                 int32 *permI, int32 n_permI)

@cython.boundscheck(False)
def rcm(mtx not None):
    """
    Generate the reversed Cuthil-McKee permutation for a CSR matrix.
    """
    cdef np.ndarray[int32, mode='c', ndim=1] perm
    cdef np.ndarray[int32, mode='c', ndim=1] indices
    cdef np.ndarray[int32, mode='c', ndim=1] indptr
    cdef int32 *_perm, *_indptr, *_indices
    cdef int32 n_row, n_ptr, n_indx

    n_row = mtx.shape[0]
    assert n_row == mtx.shape[1]
    perm = np.zeros((n_row,), dtype=np.int32)
    indptr = mtx.indptr
    indices = mtx.indices

    _perm = &perm[0]
    _indptr = &indptr[0]
    n_ptr = indptr.shape[0]
    _indices = &indices[0]
    n_indx = indices.shape[0]

    rcm_genrcm(_perm, n_row, _indptr, n_ptr, _indices, n_indx)

    return perm

@cython.boundscheck(False)
def permute_in_place(mtx not None,
                     np.ndarray[int32, mode='c', ndim=1] perm not None,
                     np.ndarray[int32, mode='c', ndim=1] perm_i=None,
                     inverse=False):
    """
    Permute a graph (= CSR sparse matrix with boolean values) in place,
    given a permuation vector.
    """
    cdef int32 ret = 0
    cdef np.ndarray[int32, mode='c', ndim=1] indices
    cdef np.ndarray[int32, mode='c', ndim=1] indptr
    cdef int32 *_perm, *_perm_i, *_indptr, *_indices
    cdef int32 n_row, n_ptr, n_indx

    if perm_i is None:
        perm_i = np.empty_like(perm)
        perm_i[perm] = np.arange(perm.shape[0], dtype=perm_i.dtype)

    if inverse:
        perm, perm_i = perm_i, perm.copy()
    else:
        perm_i = perm_i.copy()

    indptr = mtx.indptr
    indices = mtx.indices

    _perm = &perm[0]
    n_row = perm.shape[0]
    _perm_i = &perm_i[0]
    _indptr = &indptr[0]
    n_ptr = indptr.shape[0]
    _indices = &indices[0]
    n_indx = indices.shape[0]

    assert n_row == perm_i.shape[0]
    assert n_row == mtx.shape[0]
    assert n_row == mtx.shape[1]

    # Destroys perm_i!
    ret = gr_permuteInPlace(_indptr, n_ptr, _indices, n_indx,
                            _perm, n_row, _perm_i, n_row)
    if ret:
        errclear()
        raise ValueError('ccore error (see above)')
