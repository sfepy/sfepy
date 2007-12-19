from sfe.base.testing import TestCommon
import numpy as nm
import scipy.sparse as sp
import os.path as op

##
# 02.07.2007, c
class Test( TestCommon ):

    ##
    # 02.07.2007, c
    def fromConf( conf, options ):
        return Test( conf = conf, options = options )
    fromConf = staticmethod( fromConf )

    ##
    # 02.07.2007, c
    # 03.07.2007
    def test_sparseMatrixHDF5( self ):
        from sfe.base.ioutils import writeSparseMatrixHDF5, readSparseMatrixHDF5
        fileName = op.join( self.options.outDir, 'mtx.h5' )

        aux = nm.random.rand( 5, 5 )
        aux[1,:] = aux[:,2] = aux[3,:] = 0.0

        mtx = sp.csr_matrix( aux, dtype = nm.float64 )
#        self.report( 'sparse matrix:\n%s' % mtx )
        self.report( 'saving matrix into %s...' % fileName )
        writeSparseMatrixHDF5( fileName, mtx )
        self.report( 'reading...' )
        mtx2 = readSparseMatrixHDF5( fileName )
#        self.report( 'read matrix:\n%s' % mtx2 )
        self.report( 'difference:\n%s' % (mtx2 - mtx).__repr__() )

        assert mtx.shape == mtx2.shape
        assert mtx.dtype == mtx2.dtype
        assert mtx.format == mtx2.format
        assert nm.allclose( mtx.data, mtx2.data )
        assert nm.allclose( mtx.indices, mtx2.indices )
        assert nm.allclose( mtx.indptr, mtx2.indptr )

        return True

    ##
    # 09.07.2007, c
    def test_recursiveDictHDF5( self ):
        from sfe.base.ioutils import writeDictHDF5, readDictHDF5
        fileName = op.join( self.options.outDir, 'dict.h5' )

        test = {'A' : 0, 'B' : {'C' : [0, 1],
                                'D' : {'E' : {'F' : {'G' : 2.0}}}}}

        self.report( '%s' % test )
        self.report( 'saving into %s...' % fileName )
        writeDictHDF5( fileName, test )
        self.report( 'reading...' )
        test2 = readDictHDF5( fileName )
        self.report( '%s' % test2 )

        assert test == test2

        return True
