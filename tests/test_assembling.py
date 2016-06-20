from __future__ import absolute_import
import numpy as nm
import scipy.sparse as sps

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        conn = nm.array([[0, 1, 2],
                         [2, 3, 4]], dtype=nm.int32)

        num = conn.max() + 1

        iels = nm.array([0, 1], dtype=nm.int32)

        vec_in_els = nm.array([[1, 1, 1], [2, 2, 2]], dtype=nm.float64)
        vec_in_els.shape = (2, 1, 1, 3)

        mtx_in_els = nm.array([nm.ones((3, 3)), 2 * nm.ones((3, 3))],
                              dtype=nm.float64)
        mtx_in_els.shape = (2, 1, 3, 3)

        return Test(conn=conn, num=num, iels=iels,
                    vec_in_els=vec_in_els, mtx_in_els=mtx_in_els,
                    conf=conf, options=options)

    def test_assemble_vector(self):
        from sfepy.discrete.common.extmods.assemble import assemble_vector

        vec = nm.zeros(self.num, dtype=nm.float64)

        assemble_vector(vec, self.vec_in_els, self.iels, 1, self.conn)

        aux = nm.array([1, 1, 3, 2, 2], dtype=nm.float64)

        self.report('assembled: %s' % vec)
        self.report('expected: %s' % aux)
        ok = self.compare_vectors(vec, aux,
                                  label1='assembled',
                                  label2='expected')
        return ok

    def test_assemble_vector_complex(self):
        from sfepy.discrete.common.extmods.assemble import assemble_vector_complex

        vec = nm.zeros(self.num, dtype=nm.complex128)
        vec_in_els = self.vec_in_els.astype(nm.complex128) * (2 - 3j)

        assemble_vector_complex(vec, vec_in_els, self.iels, 1, self.conn)

        aux = nm.array([2-3j, 2-3j, 6-9j, 4-6j, 4-6j],
                       dtype=nm.complex128)

        self.report('assembled: %s' % vec)
        self.report('expected: %s' % aux)
        ok = self.compare_vectors(vec, aux,
                                  label1='assembled',
                                  label2='expected')
        return ok

    def test_assemble_matrix(self):
        from sfepy.discrete.common.extmods.assemble import assemble_matrix

        mtx = sps.csr_matrix(nm.ones((self.num, self.num),
                                     dtype=nm.float64))
        mtx.data[:] = 0.0

        assemble_matrix(mtx.data, mtx.indptr, mtx.indices, self.mtx_in_els,
                        self.iels, 1, self.conn, self.conn)

        aux = nm.array([[1, 1, 1, 0, 0],
                        [1, 1, 1, 0, 0],
                        [1, 1, 3, 2, 2],
                        [0, 0, 2, 2, 2],
                        [0, 0, 2, 2, 2]], dtype=nm.float64)

        self.report('assembled:\n%s' % mtx.toarray())
        self.report('expected:\n%s' % aux)
        ok = self.compare_vectors(mtx, aux,
                                  label1='assembled',
                                  label2='expected')
        return ok

    def test_assemble_matrix_complex(self):
        from sfepy.discrete.common.extmods.assemble import assemble_matrix_complex

        mtx = sps.csr_matrix(nm.ones((self.num, self.num),
                                     dtype=nm.complex128))
        mtx.data[:] = 0.0
        mtx_in_els = self.mtx_in_els.astype(nm.complex128) * (2 - 3j)

        assemble_matrix_complex(mtx.data, mtx.indptr, mtx.indices, mtx_in_els,
                                self.iels, 1, self.conn, self.conn)

        aux = nm.array([[2-3j, 2-3j, 2-3j, 0+0j, 0+0j],
                        [2-3j, 2-3j, 2-3j, 0+0j, 0+0j],
                        [2-3j, 2-3j, 6-9j, 4-6j, 4-6j],
                        [0+0j, 0+0j, 4-6j, 4-6j, 4-6j],
                        [0+0j, 0+0j, 4-6j, 4-6j, 4-6j]], dtype=nm.complex128)

        self.report('assembled:\n%s' % mtx.toarray())
        self.report('expected:\n%s' % aux)
        ok = self.compare_vectors(mtx, aux,
                                  label1='assembled',
                                  label2='expected')
        return ok
