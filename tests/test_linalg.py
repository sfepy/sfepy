from __future__ import absolute_import
from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_tensors(self):
        import numpy as nm
        from sfepy.linalg import dot_sequences, insert_strided_axis

        ok = True

        a = nm.arange(1, 10).reshape(3, 3)
        b = nm.arange(9, 0, -1).reshape(3, 3)

        dab = nm.dot(a, b)
        dabt = nm.dot(a, b.T)
        datb = nm.dot(a.T, b)
        datbt = nm.dot(a.T, b.T)

        sa = insert_strided_axis(a, 0, 10)
        sb = insert_strided_axis(b, 0, 10)

        dsab = dot_sequences(sa, sb, mode='AB')
        _ok = nm.allclose(dab[None, ...], dsab, rtol=0.0, atol=1e-14)
        self.report('dot_sequences AB: %s' % _ok)
        ok = ok and _ok

        dsabt = dot_sequences(sa, sb, mode='ABT')
        _ok = nm.allclose(dabt[None, ...], dsabt, rtol=0.0, atol=1e-14)
        self.report('dot_sequences ABT: %s' % _ok)
        ok = ok and _ok

        dsatb = dot_sequences(sa, sb, mode='ATB')
        _ok = nm.allclose(datb[None, ...], dsatb, rtol=0.0, atol=1e-14)
        self.report('dot_sequences ATB: %s' % _ok)
        ok = ok and _ok

        dsatbt = dot_sequences(sa, sb, mode='ATBT')
        _ok = nm.allclose(datbt[None, ...], dsatbt, rtol=0.0, atol=1e-14)
        self.report('dot_sequences ATBT: %s' % _ok)
        ok = ok and _ok

        return ok

    def test_unique_rows(self):
        import numpy as nm
        from sfepy.linalg import unique_rows

        a = nm.arange(1, 10).reshape(3, 3)

        b = nm.r_[a, a]
        c = unique_rows(b)

        ok = (a == c).all()

        return ok

    def test_assemble1d(self):
        import numpy as nm
        from sfepy.linalg import assemble1d

        a = nm.arange(5)
        b = nm.arange(2)

        assemble1d(b, [1, 1, 1, 1, 0, 0], a[[0, 2, 3, 4, 1, 1]])

        ok = (b == [2, 10]).all()

        return ok

    def test_geometry(self):
        import numpy as nm
        from sfepy.linalg import get_face_areas

        a1 = get_face_areas([[0, 1, 2, 3]],
                            [[0, 0], [1, 0], [1, 1], [0, 1]])

        a2 = get_face_areas([[0, 1, 2, 3]],
                            [[0, 0, 2], [1, 0, 2], [1, 1, 2], [0, 1, 2]])
        ok = nm.allclose([a1, a2], [1, 1], rtol=0, atol=1e-15)

        return ok
