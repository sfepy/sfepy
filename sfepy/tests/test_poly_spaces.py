"""
Test continuity of polynomial basis and its gradients along an edge on
:math:`y` line (2D) or on a face in :math:`x`-:math:`y` plane (3D) between two
elements aligned with the coordinate system, stack one on top of the other. The
evaluation occurs in several points shifted by a very small amount from the
boundary between the elements into the top and the bottom element.

For H1 space, the basis should be continuous. The components of its gradient
parallel to the edge/face should be continuous as well, while the perpendicular
component should have the same absolute value, but different sign in the top
and the bottom element.

All connectivity permutations of the two elements are tested.

The serendipity basis implementation is a pure python proof-of-concept. Its
order in continuity tests is limited to 2 on 3_8 elements to decrease the tests
run time.
"""
from itertools import product
import numpy as nm
import pytest

from sfepy.base.base import assert_
import sfepy.base.testing as tst

rsels = {
    '2_3' : 'vertices in (y > -0.1) & (y < 0.1)',
    '2_4' : 'vertices in (y > 0.9) & (y < 1.1)',
    '3_4' : 'vertices in (z > -0.1) & (z < 0.1)',
    '3_8' : 'vertices in (z > 0.9) & (z < 1.1)',
}

eps = 1e-5

shifts = {
    '2_3' : nm.array([[0.0, 0.0], [0.0, eps]], dtype=nm.float64),
    '2_4' : nm.array([[0.0, 1.0], [0.0, eps]], dtype=nm.float64),
    '3_4' : nm.array([[0.0, 0.0, 0.0], [0.0, 0.0, eps]], dtype=nm.float64),
    '3_8' : nm.array([[0.0, 0.0, 1.0], [0.0, 0.0, eps]], dtype=nm.float64),
}

def _permute_quad_face(mesh0):
    from sfepy.discrete.fem import Mesh

    coors0, ngroups0, conns0, mat_ids0, desc0 = mesh0._get_io_data()
    im = nm.arange(4, 8)
    cperms = [[0, 1, 2, 3],
              [0, 2, 1, 3],
              [0, 1, 3, 2]]
    meshes = []
    for ip, cperm in enumerate(cperms):
        tr = nm.arange(12)
        tr[im] = tr[im[cperm]]

        coors = coors0.copy()
        coors[im] = coors0[im[cperm]]
        conn = conns0[0].copy()
        conn = tr[conn]

        mesh = Mesh.from_data(mesh0.name + str(ip), coors, ngroups0,
                              [conn], mat_ids0, desc0)
        meshes.append(mesh)

    return meshes

def _get_possible_oris(geom):
    import sfepy.discrete.fem.facets as facets

    if geom in ('2_3', '2_4'):
        oris = facets.ori_line_to_iter.keys()

    elif geom == '3_4':
        oris = facets.ori_triangle_to_iter.keys()

    elif geom == '3_8':
        oris = facets._quad_ori_groups.keys() # Not ori_square_to_iter!

    else:
        oris = []

    return set(oris)

def _gen_common_data(orders, gels):
    import sfepy
    from sfepy.base.base import Struct
    from sfepy.discrete import FieldVariable, Integral
    from sfepy.discrete.fem import Mesh, FEDomain, Field
    from sfepy.discrete.common.global_interp import get_ref_coors

    bases = ([ii for ii in product(['2_4', '3_8'],
                                   ['lagrange', 'serendipity', 'bernstein',
                                    'lobatto', 'sem'])]
             + [ii for ii in product(['2_3', '3_4'],
                                     ['lagrange', 'bernstein'])])
    for geom, poly_space_basis in bases:
        order = orders[geom]
        if (geom == '3_8') and (poly_space_basis == 'serendipity'):
            order = 2

        tst.report('geometry: %s, base: %s, order: %d'
                   % (geom, poly_space_basis, order))

        integral = Integral('i', order=order)

        aux = '' if geom in ['2_4', '3_8'] else 'z'
        mesh0 = Mesh.from_file('meshes/elements/%s_2%s.mesh' % (geom, aux),
                               prefix_dir=sfepy.data_dir)
        if (geom == '3_8'):
            meshes = _permute_quad_face(mesh0)

        else:
            meshes = [mesh0]

        gel = gels[geom]

        perms = gel.get_conn_permutations()

        qps, qp_weights = integral.get_qp(gel.surface_facet.name)
        zz = nm.zeros_like(qps[:, :1])
        qps = nm.hstack(([qps] + [zz]))

        shift = shifts[geom]
        rcoors = nm.ascontiguousarray(qps
                                      + shift[:1, :] - shift[1:, :])
        ccoors = nm.ascontiguousarray(qps
                                      + shift[:1, :] + shift[1:, :])

        all_oris = _get_possible_oris(geom)
        oris = set()

        for (ir, pr), (ic, pc), (im, mesh0) in product(
                enumerate(perms), enumerate(perms), enumerate(meshes),
        ):
            tst.report('im: %d, ir: %d, ic: %d' % (im, ir, ic))
            tst.report('pr: %s, pc: %s' % (pr, pc))

            mesh = mesh0.copy()
            conn = mesh.cmesh.get_conn(mesh0.cmesh.tdim, 0).indices
            conn = conn.reshape((mesh0.n_el, -1))
            conn[0, :] = conn[0, pr]
            conn[1, :] = conn[1, pc]

            conn2 = mesh.get_conn(gel.name)
            assert_((conn == conn2).all())

            cache = Struct(mesh=mesh)

            domain = FEDomain('domain', mesh)
            omega = domain.create_region('Omega', 'all')
            region = domain.create_region('Facet', rsels[geom], 'facet')
            field = Field.from_args('f', nm.float64, shape=1,
                                    region=omega, approx_order=order,
                                    poly_space_basis=poly_space_basis)
            fis = region.get_facet_indices()
            conn = mesh.cmesh.get_conn_as_graph(region.dim,
                                                region.dim - 1)
            _oris = mesh.cmesh.facet_oris[conn.indptr[fis[:, 0]]
                                          + fis[:, 1]]
            oris |= set(_oris)
            if oris == all_oris:
                break

            var = FieldVariable('u', 'unknown', field)
            tst.report('# dofs: %d' % var.n_dof)

            vec = nm.empty(var.n_dof, dtype=var.dtype)

            ps = field.poly_space

            dofs = field.get_dofs_in_region(region, merge=False)
            edofs, fdofs = nm.unique(dofs[1]), nm.unique(dofs[2])

            rrc, rcells, rstatus = get_ref_coors(field, rcoors,
                                                 cache=cache)
            crc, ccells, cstatus = get_ref_coors(field, ccoors,
                                                 cache=cache)
            assert_((rstatus == 0).all() and (cstatus == 0).all())

            yield (geom, poly_space_basis, qp_weights, mesh, im, ir, ic,
                   field, ps, rrc, rcells[0], crc, ccells[0],
                   vec, edofs, fdofs)

@pytest.fixture(scope='module')
def gels():
    from sfepy.discrete.fem.geometry_element import GeometryElement

    gels = {}
    for key in ['2_3', '2_4', '3_4', '3_8']:
        gel = GeometryElement(key)
        gel.create_surface_facet()
        gels[key] = gel

    return gels

def test_partition_of_unity(gels):
    from sfepy.discrete import Integral, PolySpace

    ok = True
    orders = {'2_3' : 5, '2_4' : 5, '3_4' : 5, '3_8' : 5}
    bases = (
        [ii for ii in product(['2_4', '3_8'],
                              ['lagrange', 'serendipity', 'bernstein', 'sem']
        )]
        + [ii for ii in product(['2_3', '3_4'],
                                ['lagrange', 'bernstein'])]
    )

    for geom, poly_space_basis in bases:
        max_order = orders[geom]
        for order in range(max_order + 1):
            if (poly_space_basis == 'serendipity') and not (0 < order < 4):
                continue
            if (poly_space_basis == 'sem') and not (0 < order):
                continue
            tst.report('geometry: %s, base: %s, order: %d'
                       % (geom, poly_space_basis, order))

            integral = Integral('i', order=2 * order)
            coors, _ = integral.get_qp(geom)

            ps = PolySpace.any_from_args('ps', gels[geom], order,
                                         basis=poly_space_basis)
            vals = ps.eval_basis(coors)
            _ok = nm.allclose(vals.sum(axis=-1), 1, atol=1e-14, rtol=0.0)
            tst.report('partition of unity:', _ok)
            ok = ok and _ok

    assert_(ok)

@pytest.mark.slow
def test_continuity(gels):
    ok = True
    orders = {'2_3' : 3, '2_4' : 3, '3_4' : 4, '3_8' : 3}

    bads = []
    bad_families = set()
    for (geom, poly_space_basis, qp_weights, mesh, im, ir, ic,
         field, ps, rrc, rcell, crc, ccell, vec,
         edofs, fdofs) in _gen_common_data(orders, gels):

        if poly_space_basis in ('lagrange', 'serendipity', 'bernstein', 'sem'):
            rbf = ps.eval_basis(rrc)
            cbf = ps.eval_basis(crc)

        else:
            rbf = ps.eval_basis(rrc, ori=field.ori[:1])
            cbf = ps.eval_basis(crc, ori=field.ori[1:])

        dofs = nm.r_[edofs, fdofs]

        res = nm.zeros((2, dofs.shape[0]), dtype=nm.int32)
        res[0, :] = dofs
        for ii, ip in enumerate(dofs):
            vec.fill(0.0)
            vec[ip] = 1.0

            evec = vec[field.econn]

            rvals = nm.dot(rbf, evec[rcell])
            cvals = nm.dot(cbf, evec[ccell])

            _ok = nm.allclose(rvals, cvals, atol=1e-14, rtol=0.0)
            res[1, ii] = _ok
            if not _ok:
                bads.append([geom, poly_space_basis, im, ir, ic, ip])
                bad_families.add((geom, poly_space_basis))

            ok = ok and _ok

        tst.report('results (dofs, status: 1 ok, 0 failure):\n%s' % res)

    if not ok:
        tst.report('continuity errors:\n', bads)
        tst.report('%d in total!' % len(bads))
        tst.report('continuity errors occurred in these spaces:\n',
                   bad_families)

    assert_(ok)

@pytest.mark.slow
def test_gradients(gels):
    from sfepy.discrete.fem.mappings import FEMapping

    ok = True
    orders = {'2_3' : 3, '2_4' : 3, '3_4' : 4, '3_8' : 3}

    bads = []
    bad_families = set()
    for (geom, poly_space_basis, qp_weights, mesh, im, ir, ic,
         field, ps, rrc, rcell, crc, ccell, vec,
         edofs, fdofs) in _gen_common_data(orders, gels):
        gel = gels[geom]
        conn = mesh.get_conn(gel.name)

        geo_ps = field.gel.poly_space
        rmapping = FEMapping(mesh.coors, conn[rcell:rcell+1],
                             poly_space=geo_ps)
        rori = field.ori[:1] if field.ori is not None else None
        rvg = rmapping.get_mapping(rrc, qp_weights, poly_space=ps, ori=rori)
        rbfg = rvg.bfg

        cmapping = FEMapping(mesh.coors, conn[ccell:ccell+1],
                             poly_space=geo_ps)
        cori = field.ori[1:] if field.ori is not None else None
        cvg = cmapping.get_mapping(crc, qp_weights, poly_space=ps, ori=cori)
        cbfg = cvg.bfg

        dofs = nm.r_[edofs, fdofs]

        res = nm.zeros((2, dofs.shape[0]), dtype=nm.int32)
        res[0, :] = dofs
        for ii, ip in enumerate(dofs):
            vec.fill(0.0)
            vec[ip] = 1.0

            evec = vec[field.econn]

            rvals = nm.dot(rbfg, evec[rcell])[0]
            cvals = nm.dot(cbfg, evec[ccell])[0]

            okx = nm.allclose(rvals[:, 0], cvals[:, 0],
                              atol=1e-12, rtol=0.0)
            if gel.dim == 2:
                oky = nm.allclose(rvals[:, 1], -cvals[:, 1],
                                  atol=1e-12, rtol=0.0)
                _ok = okx and oky

            else:
                oky = nm.allclose(rvals[:, 1], cvals[:, 1],
                                  atol=1e-12, rtol=0.0)
                okz = nm.allclose(rvals[:, 2], -cvals[:, 2],
                                  atol=1e-12, rtol=0.0)
                _ok = okx and oky and okz

            res[1, ii] = _ok
            if not _ok:
                bads.append([geom, poly_space_basis, im, ir, ic, ip])
                bad_families.add((geom, poly_space_basis))

            ok = ok and _ok

        tst.report('results (dofs, status: 1 ok, 0 failure):\n%s' % res)

    if not ok:
        tst.report('gradient continuity errors:\n', bads)
        tst.report('%d in total!' % len(bads))
        tst.report('gradient continuity errors occurred in these'
                   ' spaces:\n', bad_families)

    assert_(ok)

def test_hessians(gels):
    """
    Test the second partial derivatives of basis functions using finite
    differences.
    """
    from sfepy.discrete import Integral, PolySpace

    ok = True
    orders = {'2_3' : 3, '2_4' : 3, '3_4' : 4, '3_8' : 3}
    bases = ([ii for ii in product(['2_3', '2_4', '3_4', '3_8'],
                                   ['lagrange'])])

    for geom, poly_space_basis in bases:
        tst.report('geometry: %s, base: %s' % (geom, poly_space_basis))
        order = orders[geom]

        integral = Integral('i', order=order)
        coors, _ = integral.get_qp(geom)

        ps = PolySpace.any_from_args('ps', gels[geom], order,
                                     basis=poly_space_basis)

        dim = coors.shape[1]
        h1 = nm.zeros((coors.shape[0], dim, dim, ps.n_nod), nm.float64)
        eps = 1e-8
        for ir in range(dim):
            cc = coors.copy()

            cc[:, ir] -= eps
            aux0 = ps.eval_basis(cc, diff=1)

            cc[:, ir] += 2 * eps
            aux1 = ps.eval_basis(cc, diff=1)

            h1[:, :, ir, :] = 0.5 * (aux1 - aux0) / eps

        h2 = ps.eval_basis(coors, diff=2)

        _ok = nm.allclose(h1, h2, rtol=0, atol=50*eps)
        tst.report('hessians: error: %.2e ok: %s'
                   % (nm.abs(h1 - h2).max(), _ok))
        ok = ok and _ok

    assert_(ok)
