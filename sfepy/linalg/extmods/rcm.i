/* -*- C -*- */
#ifdef SWIGPYTHON

%module rcm
%{
#include "rcm.h"
%}

%include "types.h"

%include common.i
%include array.i

%apply (int32 *array, int32 len) {
  (int32 *perm, int32 neqns),
  (int32 *xadj, int32 n_ptr),
  (int32 *adjncy, int32 n_indx),
  (int32 *row, int32 n_row),
  (int32 *col, int32 n_col),
  (int32 *perm, int32 n_perm),
  (int32 *permI, int32 n_permI)
};

void rcm_genrcm(int32 *perm, int32 neqns, int32 *xadj, int32 n_ptr,
		int32 *adjncy, int32 n_indx);

int32 gr_permuteInPlace(int32 *row, int32 n_row,
			int32 *col, int32 n_col,
			int32 *perm, int32 n_perm,
			int32 *permI, int32 n_permI);

%pythoncode %{
import numpy as nm

def rcm(mtx):
    """
    Generate the reversed Cuthil-McKee permutation for a CSR matrix.
    """
    n_row = mtx.shape[0]
    perm = nm.zeros((n_row,), dtype=nm.int32)
    rcm_genrcm(perm, mtx.indptr, mtx.indices)

    return perm

def permute_in_place(mtx, perm, perm_i=None, inverse=False):
    """Permute a graph (= CSR sparse matrix with boolean values) in place,
    given a permuation vector."""
    if perm_i is None:
        perm_i = nm.empty_like(perm)
        perm_i[perm] = nm.arange(perm.shape[0], dtype=perm_i.dtype)

    if inverse:
        perm, perm_i = perm_i, perm.copy()
    else:
        perm_i = perm_i.copy()

    # Destroys perm_i!
    gr_permuteInPlace(mtx.indptr, mtx.indices, perm, perm_i)

%}

#endif
