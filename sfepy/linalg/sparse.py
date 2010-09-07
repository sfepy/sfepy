"""Some sparse matrix utilities missing in scipy."""
from sfepy.base.base import *

def save_sparse_txt(filename, mtx, fmt='%d %d %f\n'):
    """Save a CSR/CSC sparse matrix into a text file"""
    fd = open(filename, 'w')

    fd.write('%d %d\n' % mtx.shape)
    fd.write('%d\n' % mtx.size)

    if mtx.format == 'csr':
        indptr, indices, data = mtx.indptr, mtx.indices, mtx.data
        for ir in xrange(mtx.shape[0]):
            for ii in xrange(indptr[ir], indptr[ir+1]):
                fd.write(fmt % (ir, indices[ii], data[ii]))

    elif mtx.format == 'csc':
        indptr, indices, data = mtx.indptr, mtx.indices, mtx.data
        for ic in xrange(mtx.shape[0]):
            for ii in xrange(indptr[ir], indptr[ir+1]):
                fd.write(fmt % (indices[ii], ic, data[ii]))

    else:
        raise ValueError('matrix format not supported! (%s)' % mtx.format)

def insert_sparse_to_csr(mtx1, mtx2, irs, ics):
    """
    Insert a sparse matrix `mtx2` into a CSR sparse matrix `mtx1` at
    rows `irs` and columns `ics`. The submatrix `mtx1[irs,ics]` must
    already be preallocated and have the same structure as `mtx2`.
    """
    import sfepy.fem.extmods.fem as fem

    if isinstance(irs, slice):
        irs = nm.arange(irs.start, irs.stop, irs.step, dtype=nm.int32)

    if isinstance(ics, slice):
        ics = nm.arange(ics.start, ics.stop, ics.step, dtype=nm.int32)

    n_row, n_col = mtx1.shape

    assert_((irs.min() >= 0) and (irs.max() < n_row))
    assert_((ics.min() >= 0) and (ics.max() < n_col))

    aux = mtx2.tocoo()
    data = nm.ascontiguousarray(aux.data[:,None,None,None])
    rows = irs[aux.row[:,None]]
    cols = ics[aux.col[:,None]]

    iels = nm.arange(rows.shape[0], dtype=nm.int32)
    fem.assemble_matrix(mtx1.data, mtx1.indptr, mtx1.indices, data,
                        iels, 1.0, rows, cols)
