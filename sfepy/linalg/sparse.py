"""Some sparse matrix utilities missing in scipy."""

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
