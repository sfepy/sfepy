from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_tensors(self):
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
