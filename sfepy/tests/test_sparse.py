import sfepy.base.testing as tst

def test_compose_sparse():
    import numpy as nm
    import scipy.sparse as sps
    from sfepy.linalg import compose_sparse

    ok = True

    # basic
    ma = sps.csr_matrix([[1, 0], [0, 1]])
    mb = sps.coo_matrix([[1, 1]])
    mk = compose_sparse([[ma, mb.T], [mb, 0]])
    expected = nm.array([[1, 0, 1],
                         [0, 1, 1],
                         [1, 1, 0]])

    _ok = nm.all(mk.toarray() == expected)
    tst.report('basic: %s' % _ok)
    ok = ok and _ok

    # sizes and slices
    ma = sps.csr_matrix([[2, 3]])
    mb = sps.coo_matrix([[4, 5, 6]])

    mk = compose_sparse([[ma, mb]], col_sizes=[2, 3])
    expected = nm.array([[2, 3, 4, 5, 6]])

    _ok = nm.all(mk.toarray() == expected)
    tst.report('sizes: %s' % _ok)
    ok = ok and _ok

    i1 = slice(1, 3)
    i2 = slice(8, 11)
    mk = compose_sparse([[ma, mb]], col_sizes=[i1, i2])
    expected = nm.array([[0, 2, 3, 0, 0, 0, 0, 0, 4, 5, 6]])

    _ok = nm.all(mk.toarray() == expected)
    tst.report('slices: %s' % _ok)
    ok = ok and _ok

    # zero block sizes and slices
    mk = compose_sparse([[0, ma, 0, mb, 0]], col_sizes=[1, 2, 5, 3, 1])
    expected = nm.array([[0, 2, 3, 0, 0, 0, 0, 0, 4, 5, 6, 0]])

    _ok = nm.all(mk.toarray() == expected)
    tst.report('zero block sizes: %s' % _ok)
    ok = ok and _ok

    expected = nm.array([[0, 2, 3, 0, 4, 5, 6, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 6]])

    i0 = slice(0, 1)
    i1 = slice(1, 3)
    i2 = slice(4, 7)
    i3 = slice(8, 11)
    mk = compose_sparse([[0, ma, mb, 0],
                         [0, 0, 0, mb]], col_sizes=[i0, i1, i2, i3])

    _ok = nm.all(mk.toarray() == expected)
    tst.report('zero block slices: %s' % _ok)
    ok = ok and _ok

    assert ok
