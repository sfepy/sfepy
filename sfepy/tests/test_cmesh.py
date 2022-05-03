import os

import numpy as nm
import pytest

import sfepy.base.testing as tst
from sfepy import data_dir

# n_vertex, n_edge, n_face, n_cell
# d1 -> d2 : num, n_incident
expected = {
    '1_2_2.mesh' : ([3, 2, 0, 0], {
        (0, 0) : (3, 4),
        (0, 1) : (3, 4),
        (1, 0) : (2, 4),
        (1, 1) : (2, 2),
        }),
    '2_3_2.mesh' : ([4, 5, 2, 0], {
        (0, 0) : (4, 10),
        (0, 1) : (4, 10),
        (0, 2) : (4, 6),
        (1, 0) : (5, 10),
        (1, 1) : (5, 16),
        (1, 2) : (5, 6),
        (2, 0) : (2, 6),
        (2, 1) : (2, 6),
        (2, 2) : (2, 2),
        }),
    '2_4_2.mesh' : ([6, 7, 2, 0], {
        (0, 0) : (6, 22),
        (0, 1) : (6, 14),
        (0, 2) : (6, 8),
        (1, 0) : (7, 14),
        (1, 1) : (7, 20),
        (1, 2) : (7, 8),
        (2, 0) : (2, 8),
        (2, 1) : (2, 8),
        (2, 2) : (2, 2),
        }),
    '3_4_2.mesh' : ([5, 9, 7, 2], {
        (0, 0) : (5, 18),
        (0, 1) : (5, 18),
        (0, 2) : (5, 21),
        (0, 3) : (5, 8),
        (1, 0) : (9, 18),
        (1, 1) : (9, 48),
        (1, 2) : (9, 21),
        (1, 3) : (9, 12),
        (2, 0) : (7, 21),
        (2, 1) : (7, 21),
        (2, 2) : (7, 42),
        (2, 3) : (7, 8),
        (3, 0) : (2, 8),
        (3, 1) : (2, 12),
        (3, 2) : (2, 8),
        (3, 3) : (2, 2),
        }),
    '3_8_2.mesh' : ([12, 20, 11, 2], {
        (0, 0) : (12, 100),
        (0, 1) : (12, 40),
        (0, 2) : (12, 44),
        (0, 3) : (12, 16),
        (1, 0) : (20, 40),
        (1, 1) : (20, 96),
        (1, 2) : (20, 44),
        (1, 3) : (20, 24),
        (2, 0) : (11, 44),
        (2, 1) : (11, 44),
        (2, 2) : (11, 72),
        (2, 3) : (11, 12),
        (3, 0) : (2, 16),
        (3, 1) : (2, 24),
        (3, 2) : (2, 12),
        (3, 3) : (2, 2),
        }),
    'square_triquad.mesh' : ([470, 1127, 658, 0], {
        (0, 0) : (470, 3054),
        (0, 1) : (470, 2254),
        (0, 2) : (470, 2174),
        (1, 0) : (1127, 2254),
        (1, 1) : (1127, 9174),
        (1, 2) : (1127, 2174),
        (2, 0) : (658, 2174),
        (2, 1) : (658, 2174),
        (2, 2) : (658, 6686),
        }),
}

@pytest.fixture(scope='module')
def filename_meshes():
    meshes = [data_dir + '/meshes/elements/%s_2.mesh' % geom
              for geom in ['1_2', '2_3', '2_4', '3_4', '3_8']]
    meshes.append(data_dir + '/meshes/2d/special/square_triquad.mesh')
    return meshes

def test_cmesh_counts(filename_meshes):
    from sfepy.discrete.fem import Mesh
    from sfepy.discrete.fem.geometry_element import create_geometry_elements
    from sfepy.discrete.common.extmods.cmesh import get_cmem_usage

    gels = create_geometry_elements()

    ok = True

    for filename in filename_meshes:
        basename = os.path.basename(filename)
        enum, esizes = expected[basename]

        tst.report('mesh: %s' % basename)

        mesh = Mesh.from_file(filename)
        cmesh = mesh.cmesh
        cmesh.set_local_entities(gels)

        cmesh.setup_entities()

        tst.report('dim:', cmesh.dim)
        tst.report('n_vertex: %d, n_edge: %d, n_face: %d, n_cell: %d' %
                   tuple(cmesh.num))

        _ok = (enum == cmesh.num).all()
        if not _ok:
            tst.report('%s == %s failed!' % (enum, cmesh.num))
        ok = ok and _ok

        dim = cmesh.dim
        for ir in range(dim + 1):
            for ic in range(dim + 1):
                cmesh.setup_connectivity(ir, ic)
                mem_usage1 = get_cmem_usage()[0]

                if (ir == dim) and (ic == 0):
                    continue

                cmesh.free_connectivity(ir, ic)
                mem_usage2 = get_cmem_usage()[0]

                cmesh.setup_connectivity(ir, ic)
                mem_usage3 = get_cmem_usage()[0]

                conn = cmesh.get_conn(ir, ic)

                tst.report('(%d, %d) : (%d, %d)'
                           % (ir, ic, conn.num, conn.n_incident))
                sizes = nm.array([conn.num, conn.n_incident])

                _ok = (esizes[ir, ic] == sizes).all()
                if not _ok:
                    tst.report('%s == %s failed!' % (esizes, sizes))
                ok = ok and _ok

                _ok1 = mem_usage3 == mem_usage1
                _ok2 = mem_usage3 > mem_usage2
                if not (_ok1 and _ok2):
                    tst.report('unexpected memory usage! (%s)'
                               % ([mem_usage1, mem_usage2, mem_usage3],))
                ok = ok and (_ok1 and _ok2)

    assert ok

def test_entity_volumes():
    import sfepy
    from sfepy.discrete.fem import Mesh, FEDomain
    from sfepy.discrete.common import Field
    from sfepy.discrete import Integral

    mesh = Mesh.from_file('meshes/3d/special/cross3d.mesh',
                          prefix_dir=sfepy.data_dir)
    domain = FEDomain('domain', mesh)

    omega = domain.create_region('Omega', 'all')
    gamma = domain.create_region('Gamma', 'vertices of surface', 'facet')
    top = domain.create_region('Top', 'cell 2')

    vfield = Field.from_args('v', nm.float64, 'scalar', omega,
                             approx_order=1)
    sfield = Field.from_args('s', nm.float64, 'scalar', gamma,
                             approx_order=1)

    integral = Integral('i', order=3)
    vgeo, _ = vfield.get_mapping(omega, integral, 'volume')
    domain.create_surface_group(gamma)
    sgeo, _ = sfield.get_mapping(gamma, integral, 'surface')

    evols = mesh.cmesh.get_volumes(1)
    fvols = mesh.cmesh.get_volumes(2) # Approximate for non-planar faces.
    cvols = mesh.cmesh.get_volumes(3)

    ok = True
    _ok = abs(cvols.sum() - vgeo.volume.sum()) < 1e-15
    tst.report('total cell volume: %s (ok: %s)' % (cvols.sum(), _ok))
    ok = _ok and ok

    top_evols = nm.array([ 1.                ,  1.                ,
                           1.                ,  1.                ,
                           0.7211102550927979,  0.7211102550927979,
                           0.7211102550927979,  0.7211102550927979,
                           1.16619037896906  ,  1.16619037896906  ,
                           1.16619037896906  ,  1.16619037896906  ])

    _ok = nm.allclose(top_evols, evols[top.edges], rtol=0.0, atol=1e-15)
    tst.report('total top cell edge length: %s (ok: %s)'
               % (evols[top.edges].sum(), _ok))
    ok = _ok and ok

    i1 = [5, 6, 8, 9]
    i2 = nm.setdiff1d(nm.arange(len(gamma.faces)), i1)
    aux = fvols[gamma.faces] - sgeo.volume.ravel()

    _ok = nm.allclose(aux[i1], 0.10560208437556773, rtol=0.0, atol=1e-15)
    ok = _ok and ok
    tst.report('non-planar faces diff: %s (ok: %s)' % (aux[i1], _ok))

    _ok = (nm.abs(aux[i2]) < 1e-15).all()
    tst.report('max. planar faces diff: %s (ok: %s)'
               % (nm.abs(aux[i2]).max(), _ok))
    ok = _ok and ok

    assert ok
