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

    def test_transform_data(self):
        import numpy as nm
        from sfepy.mechanics.tensors import transform_data

        ok = True

        coors = nm.eye(3)

        data = nm.eye(3)
        expected = nm.zeros((3, 3))
        expected[[0, 1, 2], [0, 0, 2]] = 1.0

        out = transform_data(data, coors)

        _ok = nm.allclose(out, expected, rtol=0.0, atol=1e-14)
        self.report('vectors in cylindrical coordinates: %s' % _ok)
        ok = ok and _ok

        data = nm.zeros((3, 6))
        data[:, :3] = [[1, 2, 3]]
        expected = data.copy()
        expected[1, [0, 1]] = expected[1, [1, 0]]

        out = transform_data(data, coors)

        _ok = nm.allclose(out, expected, rtol=0.0, atol=1e-14)
        self.report('sym. tensors in cylindrical coordinates: %s' % _ok)
        ok = ok and _ok

        return ok

    def test_stress_transform(self):
        import numpy as nm
        from sfepy.mechanics.tensors import StressTransform

        stress_2pk = nm.arange(6) + 1

        def_grad = nm.array([[0.5047051 , 0.71142596, 0.10180901],
                             [0.13427707, 0.87156371, 0.42612244],
                             [0.27509466, 0.6262605 , 0.87659051]])
        det = nm.linalg.det(def_grad)

        aux = stress_2pk[[0, 3, 4, 3, 1, 5, 4, 5, 2]].reshape(3, 3)
        expected = nm.dot(nm.dot(def_grad, aux), def_grad.T) / det
        expected = expected.ravel()[[0, 4, 8, 1, 2, 5]][:, None]
        expected = nm.tile(expected, (5, 1, 1, 1))

        transform = StressTransform(nm.tile(def_grad, (5, 1, 1, 1)))

        stress_2pk.shape = (6, 1)
        ts = nm.tile(stress_2pk.reshape((6, 1)), (5, 1, 1, 1))
        stress_cauchy = transform.get_cauchy_from_2pk(ts)

        ok = nm.allclose(stress_cauchy, expected, rtol=0.0, atol=1e-12)
        self.report('stress: Cauchy from second Piola-Kirchhoff: %s' % ok)

        return ok
