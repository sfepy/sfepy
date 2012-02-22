from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_elastic_constants(self):
        import numpy as nm
        from sfepy.mechanics.matcoefs import ElasticConstants

        ok = True

        names = ['bulk', 'lam', 'mu', 'young', 'poisson', 'p_wave']

        ec = ElasticConstants(lam=1.0, mu=1.5)
        vals = ec.get(names)

        self.report('using values:', vals)

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

                self.report(names[i1], names[i2], '->', _ok)
                if not _ok:
                    self.report('correct:', vals)
                    self.report('    got:', ec_vals)

                ok = ok and _ok

        return ok

    def test_conversion_functions(self):
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

        _lam, _mu = mc.youngpoisson_to_lame(young, poisson)
        _ok = (nm.allclose(lam, _lam, rtol=0.0, atol=1e-14) and
               nm.allclose(mu, _mu, rtol=0.0, atol=1e-14))
        self.report('youngpoisson_to_lame():', _ok)
        if not _ok:
            self.report('correct:', lam, mu)
            self.report('    got:', _lam, _mu)
        ok = ok and _ok

        _bulk = mc.bulk_modulus_youngpoisson(young, poisson)
        _ok = nm.allclose(bulk, _bulk, rtol=0.0, atol=1e-14)
        self.report('bulk_modulus_youngpoisson():', _ok)
        if not _ok:
            self.report('correct:', bulk)
            self.report('    got:', _bulk)
        ok = ok and _ok

        _bulk = mc.bulk_modulus_lame(lam, mu)
        _ok = nm.allclose(bulk, _bulk, rtol=0.0, atol=1e-14)
        self.report('bulk_modulus_lame():', _ok)
        if not _ok:
            self.report('correct:', bulk)
            self.report('    got:', _bulk)
        ok = ok and _ok

        return ok

    def test_stiffness_tensors(self):
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

        _ds = mc.stiffness_tensor_lame(3, lam, mu)
        assert_(_ds.shape == (3, 6, 6))

        _ok = True
        for _d in _ds:
            __ok = nm.allclose(_d, d, rtol=0.0, atol=1e-14)
            _ok = _ok and __ok
        self.report('stiffness_tensor_lame():', _ok)
        ok = ok and _ok

        d = nm.array([[ 6., -2., -2.,  0.,  0.,  0.],
                      [-2.,  6., -2.,  0.,  0.,  0.],
                      [-2., -2.,  6.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  4.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  4.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  4.]])

        _ds = mc.stiffness_tensor_lame_mixed(3, lam, mu)
        assert_(_ds.shape == (3, 6, 6))

        _ok = True
        for _d in _ds:
            __ok = nm.allclose(_d, d, rtol=0.0, atol=1e-14)
            _ok = _ok and __ok
        self.report('stiffness_tensor_lame_mixed():', _ok)
        ok = ok and _ok

        return ok
