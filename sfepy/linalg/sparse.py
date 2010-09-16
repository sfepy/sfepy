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

def compose_sparse(blocks):
    """
    Compose sparse matrices into a global sparse matrix.

    Parameters
    ----------
    blocks : sequence of sequences
        The sequence of sequences of equal lengths - the individual
        sparse matrix blocks. The integer 0 can be used to mark an all-zero
        block, if its size can be determined from the other blocks.

    Returns
    -------
    mtx : coo_matrix
        The sparse matrix (COO format) composed from the given blocks.

    Examples
    --------
    Stokes-like problem matrix.

    >>> import scipy.sparse as sp
    >>> A = sp.csr_matrix([[1, 0], [0, 1]])
    >>> B = sp.coo_matrix([[1, 1]])
    >>> K = compose_sparse_to_csr([[A, B.T], [B, 0]])
    >>> print K.todense()
    [[1 0 1]
     [0 1 1]
     [1 1 0]]
    """
    if not len(blocks):
        raise ValueError('no matrix blocks!')

    row_sizes = nm.array([-1] * len(blocks))
    col_sizes = nm.array([-1] * len(blocks[0]))

    for ir, row in enumerate(blocks):
        for ic, mtx in enumerate(row):
            if mtx == 0:
                continue

            if ic >= len(col_sizes):
                raise ValueError('invalid row size at (%d, %d)!' % (ir, ic))

            if row_sizes[ir] == -1:
                row_sizes[ir] = mtx.shape[0]

            elif row_sizes[ir] != mtx.shape[0]:
                msg = 'incompatible matrix block row size at (%d, %d)!' \
                      % (ir, ic)
                raise ValueError(msg)

            if col_sizes[ic] == -1:
                col_sizes[ic] = mtx.shape[1]

            elif col_sizes[ic] != mtx.shape[1]:
                msg = 'incompatible matrix block column size at (%d, %d)!' \
                      % (ic, ic)
                raise ValueError(msg)

    if nm.any(row_sizes == -1):
        raise ValueError('incomplete row block sizes! (%s)' % row_sizes)

    if nm.any(col_sizes == -1):
        raise ValueError('incomplete column block sizes! (%s)' % row_sizes)

    row_offsets = nm.cumsum(nm.r_[0, row_sizes])
    col_offsets = nm.cumsum(nm.r_[0, col_sizes])

    rows = []
    cols = []
    datas = []
    for ir, row in enumerate(blocks):
        for ic, mtx in enumerate(row):
            if mtx == 0:
                continue

            aux = mtx.tocoo()

            rows.append(aux.row + row_offsets[ir])
            cols.append(aux.col + col_offsets[ic])
            datas.append(aux.data)

    rows = nm.concatenate(rows)
    cols = nm.concatenate(cols)
    datas = nm.concatenate(datas)

    mtx = sp.coo_matrix((datas, (rows, cols)),
                        shape=(row_offsets[-1], col_offsets[-1]))

    return mtx
