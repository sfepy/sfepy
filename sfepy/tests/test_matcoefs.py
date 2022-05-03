import sfepy.base.testing as tst

def test_elastic_constants():
    import numpy as nm
    from sfepy.mechanics.matcoefs import ElasticConstants

    ok = True

    names = ['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave']

    ec = ElasticConstants(lam=1.0, mu=1.5)
    vals = ec.get(names)

    tst.report('using values:', vals)

    for i1 in range(len(names)):
        for i2 in range(i1+1, len(names)):
            kwargs = {names[i1] : vals[i1], names[i2] : vals[i2]}

            try:
                ec.init(**kwargs)

            except:
                _ok = False

            else:
                _ok = True

            ec_vals = ec.get(names)
            _ok = _ok and nm.allclose(ec_vals, vals)

            tst.report(names[i1], names[i2], '->', _ok)
            if not _ok:
                tst.report('correct:', vals)
                tst.report('    got:', ec_vals)

            ok = ok and _ok

    assert ok

def test_conversion_functions():
    import numpy as nm

    import sfepy.mechanics.matcoefs as mc

    ok = True

    lam = 1.0
    mu = 1.5

    ec = mc.ElasticConstants(lam=lam, mu=mu)
    young, poisson, bulk = ec.get(['young', 'poisson', 'bulk'])

    lam = nm.array([lam] * 3)
    mu = nm.array([mu] * 3)
    young = nm.array([young] * 3)
    poisson = nm.array([poisson] * 3)

    _lam, _mu = mc.lame_from_youngpoisson(young, poisson)
    _ok = (nm.allclose(lam, _lam, rtol=0.0, atol=1e-14) and
           nm.allclose(mu, _mu, rtol=0.0, atol=1e-14))
    tst.report('lame_from_youngpoisson():', _ok)
    if not _ok:
        tst.report('correct:', lam, mu)
        tst.report('    got:', _lam, _mu)
    ok = ok and _ok

    _bulk = mc.bulk_from_youngpoisson(young, poisson)
    _ok = nm.allclose(bulk, _bulk, rtol=0.0, atol=1e-14)
    tst.report('bulk_from_youngpoisson():', _ok)
    if not _ok:
        tst.report('correct:', bulk)
        tst.report('    got:', _bulk)
    ok = ok and _ok

    _bulk = mc.bulk_from_lame(lam, mu)
    _ok = nm.allclose(bulk, _bulk, rtol=0.0, atol=1e-14)
    tst.report('bulk_from_lame():', _ok)
    if not _ok:
        tst.report('correct:', bulk)
        tst.report('    got:', _bulk)
    ok = ok and _ok

    assert ok

def test_stiffness_tensors():
    import numpy as nm

    from sfepy.base.base import assert_
    import sfepy.mechanics.matcoefs as mc

    ok = True

    lam = 1.0
    mu = 4.0

    lam = nm.array([lam] * 3)
    mu = nm.array([mu] * 3)

    d = nm.array([[ 9.,  1.,  1.,  0.,  0.,  0.],
                  [ 1.,  9.,  1.,  0.,  0.,  0.],
                  [ 1.,  1.,  9.,  0.,  0.,  0.],
                  [ 0.,  0.,  0.,  4.,  0.,  0.],
                  [ 0.,  0.,  0.,  0.,  4.,  0.],
                  [ 0.,  0.,  0.,  0.,  0.,  4.]])

    _ds = mc.stiffness_from_lame(3, lam, mu)
    assert_(_ds.shape == (3, 6, 6))

    _ok = True
    for _d in _ds:
        __ok = nm.allclose(_d, d, rtol=0.0, atol=1e-14)
        _ok = _ok and __ok
    tst.report('stiffness_from_lame():', _ok)
    ok = ok and _ok

    _lam, _mu = mc.lame_from_stiffness(d)
    _ok = (_lam == 1) and (_mu == 4)
    tst.report('lame_from_stiffness():', _ok)
    ok = ok and _ok

    young = 1.0
    poisson = 0.25

    d = mc.stiffness_from_youngpoisson(3, young, poisson, plane='strain')
    _young, _poisson = mc.youngpoisson_from_stiffness(d, plane='strain')
    _ok = nm.allclose([young, poisson], [_young, _poisson],
                      rtol=0.0, atol=1e-14)
    tst.report('youngpoisson_from_stiffness(plane="strain"):', _ok)
    ok = ok and _ok

    d = mc.stiffness_from_youngpoisson(2, young, poisson, plane='stress')
    _young, _poisson = mc.youngpoisson_from_stiffness(d, plane='stress')
    _ok = nm.allclose([young, poisson], [_young, _poisson],
                      rtol=0.0, atol=1e-14)
    tst.report('youngpoisson_from_stiffness(plane="stress"):', _ok)
    ok = ok and _ok

    d = 4.0 / 3.0 * nm.array([[ 4., -2., -2.,  0.,  0.,  0.],
                              [-2.,  4., -2.,  0.,  0.,  0.],
                              [-2., -2.,  4.,  0.,  0.,  0.],
                              [ 0.,  0.,  0.,  3.,  0.,  0.],
                              [ 0.,  0.,  0.,  0.,  3.,  0.],
                              [ 0.,  0.,  0.,  0.,  0.,  3.]])

    _ds = mc.stiffness_from_lame_mixed(3, lam, mu)
    assert_(_ds.shape == (3, 6, 6))

    _ok = True
    for _d in _ds:
        __ok = nm.allclose(_d, d, rtol=0.0, atol=1e-14)
        _ok = _ok and __ok
    tst.report('stiffness_from_lame_mixed():', _ok)
    ok = ok and _ok

    blam = - mu * 2.0 / 3.0
    _ds = mc.stiffness_from_lame(3, blam, mu)
    assert_(_ds.shape == (3, 6, 6))

    _ok = True
    for _d in _ds:
        __ok = nm.allclose(_d, d, rtol=0.0, atol=1e-14)
        _ok = _ok and __ok
    tst.report('stiffness_from_lame() with modified lambda:', _ok)
    ok = ok and _ok

    assert ok
