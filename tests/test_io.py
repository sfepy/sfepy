import os.path as op
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import assert_
import sfepy.base.testing as tst

def test_sparse_matrix_hdf5(output_dir):
    from sfepy.base.ioutils import (write_sparse_matrix_hdf5,
                                    read_sparse_matrix_hdf5, pt)
    if pt is None:
        tst.report('skipped (no pytables)')
        return

    filename = op.join(output_dir, 'mtx.h5')

    aux = nm.random.rand(5, 5)
    aux[1,:] = aux[:,2] = aux[3,:] = 0.0

    mtx = sp.csr_matrix(aux, dtype = nm.float64)
    tst.report('saving matrix into %s...' % filename)
    write_sparse_matrix_hdf5(filename, mtx)
    tst.report('reading...')
    mtx2 = read_sparse_matrix_hdf5(filename)
    tst.report('difference:\n%s' % (mtx2 - mtx).__repr__())

    assert_(mtx.shape == mtx2.shape)
    assert_(mtx.dtype == mtx2.dtype)
    assert_(mtx.format == mtx2.format)
    assert_(nm.allclose(mtx.data, mtx2.data))
    assert_(nm.allclose(mtx.indices, mtx2.indices))
    assert_(nm.allclose(mtx.indptr, mtx2.indptr))

def test_recursive_dict_hdf5(output_dir):
    from sfepy.base.ioutils import write_dict_hdf5, read_dict_hdf5, pt
    if pt is None:
        tst.report('skipped (no pytables)')
        return

    filename = op.join(output_dir, 'dict.h5')

    test = {'A' : 0, 'B' : {'C' : [0, 1],
                            'D' : {'E' : {'F' : {'G' : 2.0}}}}}

    tst.report('%s' % test)
    tst.report('saving into %s...' % filename)
    write_dict_hdf5(filename, test)
    tst.report('reading...')
    test2 = read_dict_hdf5(filename)
    tst.report('%s' % test2)

    assert_(test == test2)
