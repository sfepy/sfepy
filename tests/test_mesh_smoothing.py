import os.path as op

from sfepy.base.testing import TestCommon

def dets_fast(a):
    import numpy as nm
    from numpy.linalg import lapack_lite
    from numpy.core import intc

    m = a.shape[0]
    n = a.shape[1]
    lapack_routine = lapack_lite.dgetrf
    pivots = nm.zeros((m, n), intc)
    flags = nm.arange(1, n + 1).reshape(1, -1)
    for i in xrange(m):
        tmp = a[i]
        lapack_routine(n, n, tmp, n, pivots[i], 0)

    sign = 1. - 2. * (nm.add.reduce(pivots != flags, axis=1) % 2)
    idx = nm.arange(n)
    d = a[:, idx, idx]
    absd = nm.absolute(d)
    sign *= nm.multiply.reduce(d / absd, axis=1)
    nm.log(absd, absd)
    logdet = nm.add.reduce(absd, axis=-1)

    return sign * nm.exp(logdet)

def get_volume(el, nd):
    from sfepy.mesh.mesh_tools import elems_q2t
    from sfepy.base.compat import factorial
    import numpy as nm

    dim = nd.shape[1]
    nnd = el.shape[1]

    etype = '%d_%d' % (dim, nnd)
    if etype == '2_4' or etype == '3_8':
        el = elems_q2t(el)

    nel = el.shape[0]

    mul = 1.0 / factorial(dim)
    if dim == 3:
        mul *= -1.0

    mtx = nm.ones((nel, dim + 1, dim + 1), dtype=nm.double)
    mtx[:,:,:-1] = nd[el,:]
    vols = mul * dets_fast(mtx.copy())
    vol = vols.sum()
    bc = nm.dot(vols, mtx.sum(1)[:,:-1] / nnd)

    bc /= vol

    return vol

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_mesh_smoothing(self):
        from sfepy.mesh.mesh_tools import smooth_mesh
        from sfepy.fem.mesh import Mesh

        mesh = Mesh.from_file('meshes/3d/cylinder.vtk')
        vol0 = get_volume(mesh.conns[0], mesh.coors)
        mesh.coors = smooth_mesh(mesh, n_iter=10)
        vol1 = get_volume(mesh.conns[0], mesh.coors)
        filename = op.join(self.options.out_dir, 'smoothed_cylinder.vtk')
        mesh.write(filename)
        frac = vol1 / vol0

        if (frac < 0.967) and (frac > 0.966):
            self.report('mesh smoothed')
            return True

        else:
            self.report('mesh smoothed, volume mismatch!')
            return False
