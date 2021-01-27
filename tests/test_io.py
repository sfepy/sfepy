from __future__ import absolute_import
from sfepy.base.base import assert_
from sfepy.base.testing import TestCommon
import numpy as nm
import scipy.sparse as sp
import os.path as op

##
# 02.07.2007, c
class Test( TestCommon ):

    ##
    # 02.07.2007, c
    def from_conf( conf, options ):
        return Test( conf = conf, options = options )
    from_conf = staticmethod( from_conf )

    ##
    # c: 02.07.2007, r: 12.06.2008
    def test_sparse_matrix_hdf5( self ):
        from sfepy.base.ioutils import write_sparse_matrix_hdf5, read_sparse_matrix_hdf5
        from sfepy.base.ioutils import pt
        if pt is None:
            self.report( 'skipped (no pytables)' )
            return True
        filename = op.join( self.options.out_dir, 'mtx.h5' )

        aux = nm.random.rand( 5, 5 )
        aux[1,:] = aux[:,2] = aux[3,:] = 0.0

        mtx = sp.csr_matrix( aux, dtype = nm.float64 )
#        self.report( 'sparse matrix:\n%s' % mtx )
        self.report( 'saving matrix into %s...' % filename )
        write_sparse_matrix_hdf5( filename, mtx )
        self.report( 'reading...' )
        mtx2 = read_sparse_matrix_hdf5( filename )
#        self.report( 'read matrix:\n%s' % mtx2 )
        self.report( 'difference:\n%s' % (mtx2 - mtx).__repr__() )

        assert_( mtx.shape == mtx2.shape )
        assert_( mtx.dtype == mtx2.dtype )
        assert_( mtx.format == mtx2.format )
        assert_( nm.allclose( mtx.data, mtx2.data ) )
        assert_( nm.allclose( mtx.indices, mtx2.indices ) )
        assert_( nm.allclose( mtx.indptr, mtx2.indptr ) )

        return True

    ##
    # c: 09.07.2007, r: 12.06.2008
    def test_recursive_dict_hdf5( self ):
        from sfepy.base.ioutils import write_dict_hdf5, read_dict_hdf5
        from sfepy.base.ioutils import pt
        if pt is None:
            self.report( 'skipped (no pytables)' )
            return True
        filename = op.join( self.options.out_dir, 'dict.h5' )

        test = {'A' : 0, 'B' : {'C' : [0, 1],
                                'D' : {'E' : {'F' : {'G' : 2.0}}}}}

        self.report( '%s' % test )
        self.report( 'saving into %s...' % filename )
        write_dict_hdf5( filename, test )
        self.report( 'reading...' )
        test2 = read_dict_hdf5( filename )
        self.report( '%s' % test2 )

        assert_( test == test2 )

        return True
