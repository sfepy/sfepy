import os.path as op

import sfepy.base.testing as tst

def get_volume(el, nd):
    from sfepy.mesh.mesh_tools import elems_q2t
    from sfepy.base.compat import factorial
    import numpy as nm
    from sfepy.linalg.utils import dets_fast

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
    vols = mul * dets_fast(mtx)
    vol = vols.sum()

    return vol

def test_mesh_smoothing(output_dir):
    from sfepy.mesh.mesh_tools import smooth_mesh
    from sfepy.discrete.fem.mesh import Mesh
    from sfepy import data_dir

    mesh = Mesh.from_file(data_dir + '/meshes/3d/cylinder.vtk')
    conn = mesh.get_conn('3_4')
    vol0 = get_volume(conn, mesh.coors)
    mesh.coors[:] = smooth_mesh(mesh, n_iter=10)
    vol1 = get_volume(conn, mesh.coors)
    filename = op.join(output_dir, 'smoothed_cylinder.vtk')
    mesh.write(filename)
    frac = vol1 / vol0

    if (frac < 0.967) and (frac > 0.966):
        tst.report('mesh smoothed')
        ok = True

    else:
        tst.report('mesh smoothed, volume mismatch!')
        ok = False

    assert ok
