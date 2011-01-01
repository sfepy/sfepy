from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_tensors(self):
        import numpy as nm
        import sfepy.mechanics.tensors as tn

        ok = True

        a_full = 2.0 * nm.ones((5,3,3), dtype=nm.float64)
        a_sym = 2.0 * nm.ones((5,6), dtype=nm.float64)

        _tr = nm.array([6.0] * 5, dtype=nm.float64)
        _vt_full = 2.0 * nm.tile(nm.eye(3, dtype=nm.float64), (5,1,1))
        _vt_sym = nm.tile(nm.array([2, 2, 2, 0, 0, 0], dtype=nm.float64),
                          (5,1,1))
        _dev_full = a_full - _vt_full
        _dev_sym = a_sym - _vt_sym
        _vms = 6.0 * nm.ones((5,1), dtype=nm.float64)

        tr = tn.get_trace(a_full, sym_storage=False)
        _ok = nm.allclose(tr, _tr, rtol=0.0, atol=1e-14)
        self.report('trace full: %s' % _ok)
        ok = ok and _ok
        
        tr = tn.get_trace(a_sym, sym_storage=True)
        ok = ok and nm.allclose(tr, _tr, rtol=0.0, atol=1e-14)
        self.report('trace sym: %s' % _ok)
        ok = ok and _ok

        vt = tn.get_volumetric_tensor(a_full, sym_storage=False)
        _ok = nm.allclose(vt, _vt_full, rtol=0.0, atol=1e-14)
        self.report('volumetric tensor full: %s' % _ok)
        ok = ok and _ok
        
        vt = tn.get_volumetric_tensor(a_sym, sym_storage=True)
        _ok = nm.allclose(vt, _vt_sym, rtol=0.0, atol=1e-14)
        self.report('volumetric tensor sym: %s' % _ok)
        ok = ok and _ok
        
        dev = tn.get_deviator(a_full, sym_storage=False)
        _ok = nm.allclose(dev, _dev_full, rtol=0.0, atol=1e-14)
        self.report('deviator full: %s' % _ok)
        ok = ok and _ok

        aux = (dev * nm.transpose(dev, (0, 2, 1))).sum(axis=1).sum(axis=1)
        vms2 = nm.sqrt((3.0/2.0) * aux)[:,None]
        
        dev = tn.get_deviator(a_sym, sym_storage=True)
        _ok = nm.allclose(dev, _dev_sym, rtol=0.0, atol=1e-14)
        self.report('deviator sym: %s' % _ok)
        ok = ok and _ok
        
        vms = tn.get_von_mises_stress(a_full, sym_storage=False)
        _ok = nm.allclose(vms, _vms, rtol=0.0, atol=1e-14)
        self.report('von Mises stress full: %s' % _ok)
        ok = ok and _ok

        vms = tn.get_von_mises_stress(a_sym, sym_storage=True)
        _ok = nm.allclose(vms, _vms, rtol=0.0, atol=1e-14)
        self.report('von Mises stress sym: %s' % _ok)
        ok = ok and _ok

        _ok = nm.allclose(vms2, _vms, rtol=0.0, atol=1e-14)
        self.report('von Mises stress via deviator: %s' % _ok)
        ok = ok and _ok

        return ok
