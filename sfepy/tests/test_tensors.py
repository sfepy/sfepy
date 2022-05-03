import sfepy.base.testing as tst

def get_ortho_d(phi1, phi2):
    import numpy as nm
    import sfepy.mechanics.tensors as tn

    v1 = nm.array([nm.cos(phi1), nm.sin(phi1), 0])
    v2 = nm.array([nm.cos(phi2), nm.sin(phi2), 0])
    om1 = nm.outer(v1, v1)
    om2 = nm.outer(v2, v2)

    ii = tn.get_sym_indices(3)

    o1 = om1.flat[ii]
    o2 = om2.flat[ii]

    dr = nm.outer(o1, o1) + nm.outer(o2, o2)
    return dr, v1, v2, om1, om2

def test_tensors():
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
    tst.report('trace full: %s' % _ok)
    ok = ok and _ok

    tr = tn.get_trace(a_sym, sym_storage=True)
    ok = ok and nm.allclose(tr, _tr, rtol=0.0, atol=1e-14)
    tst.report('trace sym: %s' % _ok)
    ok = ok and _ok

    vt = tn.get_volumetric_tensor(a_full, sym_storage=False)
    _ok = nm.allclose(vt, _vt_full, rtol=0.0, atol=1e-14)
    tst.report('volumetric tensor full: %s' % _ok)
    ok = ok and _ok

    vt = tn.get_volumetric_tensor(a_sym, sym_storage=True)
    _ok = nm.allclose(vt, _vt_sym, rtol=0.0, atol=1e-14)
    tst.report('volumetric tensor sym: %s' % _ok)
    ok = ok and _ok

    dev = tn.get_deviator(a_full, sym_storage=False)
    _ok = nm.allclose(dev, _dev_full, rtol=0.0, atol=1e-14)
    tst.report('deviator full: %s' % _ok)
    ok = ok and _ok

    aux = (dev * nm.transpose(dev, (0, 2, 1))).sum(axis=1).sum(axis=1)
    vms2 = nm.sqrt((3.0/2.0) * aux)[:,None]

    dev = tn.get_deviator(a_sym, sym_storage=True)
    _ok = nm.allclose(dev, _dev_sym, rtol=0.0, atol=1e-14)
    tst.report('deviator sym: %s' % _ok)
    ok = ok and _ok

    vms = tn.get_von_mises_stress(a_full, sym_storage=False)
    _ok = nm.allclose(vms, _vms, rtol=0.0, atol=1e-14)
    tst.report('von Mises stress full: %s' % _ok)
    ok = ok and _ok

    vms = tn.get_von_mises_stress(a_sym, sym_storage=True)
    _ok = nm.allclose(vms, _vms, rtol=0.0, atol=1e-14)
    tst.report('von Mises stress sym: %s' % _ok)
    ok = ok and _ok

    _ok = nm.allclose(vms2, _vms, rtol=0.0, atol=1e-14)
    tst.report('von Mises stress via deviator: %s' % _ok)
    ok = ok and _ok

    t2s = nm.arange(9).reshape(3, 3)
    t2s = (t2s + t2s.T) / 2
    t4 = tn.get_t4_from_t2s(t2s)
    expected = nm.array([[[[0, 4], [4, 2]],
                          [[4, 8], [8, 6]]],
                         [[[4, 8], [8, 6]],
                          [[2, 6],  [6, 4]]]])
    _ok = nm.allclose(t4, expected, rtol=0.0, atol=1e-14)
    tst.report('full 4D tensor from 2D matrix, 2D space: %s' % _ok)
    ok = ok and _ok

    assert ok

def test_transform_data():
    import numpy as nm
    from sfepy.mechanics.tensors import transform_data

    ok = True

    coors = nm.eye(3)

    data = nm.eye(3)
    expected = nm.zeros((3, 3))
    expected[[0, 1, 2], [0, 0, 2]] = 1.0

    out = transform_data(data, coors)

    _ok = nm.allclose(out, expected, rtol=0.0, atol=1e-14)
    tst.report('vectors in cylindrical coordinates: %s' % _ok)
    ok = ok and _ok

    data = nm.zeros((3, 6))
    data[:, :3] = [[1, 2, 3]]
    expected = data.copy()
    expected[1, [0, 1]] = expected[1, [1, 0]]

    out = transform_data(data, coors)

    _ok = nm.allclose(out, expected, rtol=0.0, atol=1e-14)
    tst.report('sym. tensors in cylindrical coordinates: %s' % _ok)
    ok = ok and _ok

    assert ok

def test_transform_data4():
    import numpy as nm
    import sfepy.mechanics.tensors as tn

    ok = True

    if not hasattr(nm, 'einsum'):
        tst.report('no numpy.einsum(), skipping!')
        return

    expected = nm.zeros((6, 6), dtype=nm.float64)
    expected[0, 0] = expected[1, 1] = 1.0

    phi = nm.deg2rad(30.)
    dr, v1, v2, om1, om2 = get_ortho_d(phi, phi + nm.deg2rad(90.))

    # Rotate coordinate system by phi.
    mtx = tn.make_axis_rotation_matrix([0., 0., 1.], phi)
    do = tn.transform_data(dr[None, ...], mtx=mtx[None, ...])

    _ok = nm.allclose(do, expected, rtol=0.0, atol=1e-14)
    tst.report('sym. 4th-th order tensor rotation: %s' % _ok)
    ok = ok and _ok

    dt, vt1, vt2, omt1, omt2 = get_ortho_d(0, nm.deg2rad(90.))

    expected1 = nm.zeros((3, 3), dtype=nm.float64)
    expected1[0, 0] = 1.0

    expected2 = nm.zeros((3, 3), dtype=nm.float64)
    expected2[1, 1] = 1.0

    omr1 = nm.einsum('pq,ip,jq->ij', om1, mtx, mtx)
    omr2 = nm.einsum('pq,ip,jq->ij', om2, mtx, mtx)

    ii = tn.get_sym_indices(3)
    jj = tn.get_full_indices(3)

    o1 = om1.flat[ii]
    o2 = om2.flat[ii]

    omr12 = tn.transform_data(o1[None,...], mtx=mtx[None, ...])[0, jj]
    omr22 = tn.transform_data(o2[None,...], mtx=mtx[None, ...])[0, jj]

    _ok1 = nm.allclose(omr1, expected1, rtol=0.0, atol=1e-14)
    _ok2 = nm.allclose(omr12, expected1, rtol=0.0, atol=1e-14)
    tst.report('einsum-transform_data compatibility 1: %s %s'
               % (_ok1, _ok2))
    ok = ok and _ok1 and _ok2

    _ok1 = nm.allclose(omr2, expected2, rtol=0.0, atol=1e-14)
    _ok2 = nm.allclose(omr22, expected2, rtol=0.0, atol=1e-14)
    tst.report('einsum-transform_data compatibility 2: %s %s'
               % (_ok1, _ok2))
    ok = ok and _ok1 and _ok2

    assert ok

def test_stress_transform():
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
    tst.report('stress: Cauchy from second Piola-Kirchhoff: %s' % ok)

    assert ok
