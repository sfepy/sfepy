"""Some sparse matrix utilities missing in scipy."""
from __future__ import absolute_import
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import assert_
from six.moves import range

def save_sparse_txt(filename, mtx, fmt='%d %d %f\n'):
    """Save a CSR/CSC sparse matrix into a text file"""
    fd = open(filename, 'w')

    fd.write('%d %d\n' % mtx.shape)
    fd.write('%d\n' % mtx.size)

    if mtx.format == 'csr':
        indptr, indices, data = mtx.indptr, mtx.indices, mtx.data
        for ir in range(mtx.shape[0]):
            for ii in range(indptr[ir], indptr[ir+1]):
                fd.write(fmt % (ir, indices[ii], data[ii]))

    elif mtx.format == 'csc':
        indptr, indices, data = mtx.indptr, mtx.indices, mtx.data
        for ic in range(mtx.shape[0]):
            for ii in range(indptr[ir], indptr[ir+1]):
                fd.write(fmt % (indices[ii], ic, data[ii]))

    else:
        raise ValueError('matrix format not supported! (%s)' % mtx.format)

def insert_sparse_to_csr(mtx1, mtx2, irs, ics):
    """
    Insert a sparse matrix `mtx2` into a CSR sparse matrix `mtx1` at
    rows `irs` and columns `ics`. The submatrix `mtx1[irs,ics]` must
    already be preallocated and have the same structure as `mtx2`.
    """
    import sfepy.discrete.common.extmods.assemble as asm

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
    asm.assemble_matrix(mtx1.data, mtx1.indptr, mtx1.indices, data,
                        iels, 1.0, rows, cols)

def _normalize_sizes(sizes):
    """
    Checks whether all the sizes are either slices or not. Transforms
    slices into their sizes.
    """
    out = []
    ns = 0
    for size in sizes:
        if isinstance(size, slice):
            size = size.stop - size.start
            ns += 1

        else:
            size = int(size)

        out.append(size)

    if ns:
        if ns != len(sizes):
            raise ValueError('cannot mix sizes with slices! (%s)' % (sizes,))

        is_slice = True

    else:
        is_slice = False

    return out, is_slice

def compose_sparse(blocks, row_sizes=None, col_sizes=None):
    """
    Compose sparse matrices into a global sparse matrix.

    Parameters
    ----------
    blocks : sequence of sequences
        The sequence of sequences of equal lengths - the individual
        sparse matrix blocks. The integer 0 can be used to mark an all-zero
        block, if its size can be determined from the other blocks.
    row_sizes : sequence, optional
        The required row sizes of the blocks. It can be either a
        sequence of non-negative integers, or a sequence of slices with
        non-negative limits. In any case the sizes have to be compatible
        with the true block sizes. This allows to extend the matrix
        shape as needed and to specify sizes of all-zero blocks.
    col_sizes : sequence, optional
        The required column sizes of the blocks. See `row_sizes`.

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
    >>> K = compose_sparse([[A, B.T], [B, 0]])
    >>> print K.todense()
    [[1 0 1]
     [0 1 1]
     [1 1 0]]
    """
    if not len(blocks):
        raise ValueError('no matrix blocks!')

    if row_sizes is None:
        row_sizes = nm.array([-1] * len(blocks))

    else:
        assert_(len(row_sizes) == len(blocks))

    if col_sizes is None:
        col_sizes = nm.array([-1] * len(blocks[0]))

    else:
        assert_(len(col_sizes) == len(blocks[0]))

    rs, is_slice_r = _normalize_sizes(row_sizes)
    cs, is_slice_c = _normalize_sizes(col_sizes)

    for ir, row in enumerate(blocks):
        for ic, mtx in enumerate(row):
            if isinstance(mtx, int) and (mtx == 0):
                continue

            if ic >= len(col_sizes):
                raise ValueError('invalid row size at (%d, %d)!' % (ir, ic))

            if rs[ir] == -1:
                rs[ir] = mtx.shape[0]

            elif rs[ir] != mtx.shape[0]:
                msg = 'incompatible matrix block row size at (%d, %d)!' \
                      % (ir, ic)
                raise ValueError(msg)

            if cs[ic] == -1:
                cs[ic] = mtx.shape[1]

            elif cs[ic] != mtx.shape[1]:
                msg = 'incompatible matrix block column size at (%d, %d)!' \
                      % (ic, ic)
                raise ValueError(msg)

    if nm.any(rs == -1):
        raise ValueError('incomplete row block sizes! (%s)' % row_sizes)

    if nm.any(cs == -1):
        raise ValueError('incomplete column block sizes! (%s)' % col_sizes)

    if is_slice_r:
        n_row = row_sizes[-1].stop
        row_offsets = nm.r_[[ii.start for ii in row_sizes], n_row]

    else:
        row_offsets = nm.cumsum(nm.r_[0, rs])
        n_row = row_offsets[-1]

    if is_slice_c:
        n_col = col_sizes[-1].stop
        col_offsets = nm.r_[[ii.start for ii in col_sizes], n_col]

    else:
        col_offsets = nm.cumsum(nm.r_[0, cs])
        n_col = col_offsets[-1]

    rows = []
    cols = []
    datas = []
    for ir, row in enumerate(blocks):
        for ic, mtx in enumerate(row):
            if isinstance(mtx, int) and (mtx == 0):
                continue

            aux = sp.coo_matrix(mtx)

            rows.append(aux.row + row_offsets[ir])
            cols.append(aux.col + col_offsets[ic])
            datas.append(aux.data)

    rows = nm.concatenate(rows)
    cols = nm.concatenate(cols)
    datas = nm.concatenate(datas)

    mtx = sp.coo_matrix((datas, (rows, cols)), shape=(n_row, n_col))

    return mtx

def infinity_norm(mtx):
    """
    Infinity norm of a sparse matrix (maximum absolute row sum).  

    Parameters
    ----------
    mtx : spmatrix or array
        The sparse matrix.
    
    Returns
    -------
    norm : float
        Infinity norm of the matrix.
    
    Notes
    -----
    - This serves as an upper bound on spectral radius.
    - CSR and CSC avoid copying `indices` and `indptr` arrays.
    - inspired by PyAMG

    See Also
    --------
    scipy.linalg.norm : dense matrix norms
    """
    ones = nm.ones(mtx.shape[1], dtype=mtx.dtype)

    if sp.isspmatrix_csr(mtx) or sp.isspmatrix_csc(mtx):
        # Avoid copying index and ptr arrays.
        abs_mtx = mtx.__class__((nm.abs(mtx.data), mtx.indices ,mtx.indptr),
                                shape=mtx.shape)
        norm = (abs_mtx * ones).max()

    elif sp.isspmatrix(mtx):
        norm = (abs(mtx) * ones).max()

    else:
        norm = nm.dot(nm.abs(mtx), ones).max()

    return norm
