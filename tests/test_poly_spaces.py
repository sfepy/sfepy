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

WARNING: Lagrange basis on 3_8 elements fails the test for order >= 3 for many
connectivity permutations!
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy.base.base import assert_

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

def _gen_common_data(orders, gels, report):
    import sfepy
    from sfepy.base.base import Struct
    from sfepy.linalg import combine
    from sfepy.discrete import FieldVariable, Integral
    from sfepy.discrete.fem import Mesh, FEDomain, Field
    from sfepy.discrete.common.global_interp import get_ref_coors

    bases = ([ii for ii in combine([['2_4', '3_8'],
                                    ['lagrange', 'lobatto']])]
             + [ii for ii in combine([['2_3', '3_4'],
                                      ['lagrange']])])
    for geom, poly_space_base in bases:
        report('geometry: %s, base: %s' % (geom, poly_space_base))

        order = orders[geom]
        integral = Integral('i', order=order)

        aux = '' if geom in ['2_4', '3_8'] else 'z'
        mesh0 = Mesh.from_file('meshes/elements/%s_2%s.mesh' % (geom, aux),
                               prefix_dir=sfepy.data_dir)
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

        for ir, pr in enumerate(perms):
            for ic, pc in enumerate(perms):
                report('ir: %d, ic: %d' % (ir, ic))
                report('pr: %s, pc: %s' % (pr, pc))

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
                                        poly_space_base=poly_space_base)
                var = FieldVariable('u', 'unknown', field)
                report('# dofs: %d' % var.n_dof)

                vec = nm.empty(var.n_dof, dtype=var.dtype)

                ps = field.poly_space

                dofs = field.get_dofs_in_region(region, merge=False)
                edofs, fdofs = nm.unique(dofs[1]), nm.unique(dofs[2])

                rrc, rcells, rstatus = get_ref_coors(field, rcoors,
                                                     cache=cache)
                crc, ccells, cstatus = get_ref_coors(field, ccoors,
                                                     cache=cache)
                assert_((rstatus == 0).all() and (cstatus == 0).all())

                yield (geom, poly_space_base, qp_weights, mesh, ir, ic,
                       field, ps, rrc, rcells[0], crc, ccells[0],
                       vec, edofs, fdofs)

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete.fem.geometry_element import GeometryElement

        gels = {}
        for key in ['2_3', '2_4', '3_4', '3_8']:
            gel = GeometryElement(key)
            gel.create_surface_facet()
            gels[key] = gel

        return Test(conf=conf, options=options, gels=gels)

    def test_continuity(self):
        ok = True
        orders = {'2_3' : 3, '2_4' : 3, '3_4' : 4, '3_8' : 3}

        bads = []
        bad_families = set()
        for (geom, poly_space_base, qp_weights, mesh, ir, ic,
             field, ps, rrc, rcell, crc, ccell, vec,
             edofs, fdofs) in _gen_common_data(orders, self.gels, self.report):

            if poly_space_base == 'lagrange':
                rbf = ps.eval_base(rrc)
                cbf = ps.eval_base(crc)

            else:
                rbf = ps.eval_base(rrc, ori=field.ori[:1])
                cbf = ps.eval_base(crc, ori=field.ori[1:])

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
                    bads.append([geom, poly_space_base, ir, ic, ip])
                    bad_families.add((geom, poly_space_base))

                ok = ok and _ok

            self.report('results (dofs, status: 1 ok, 0 failure):\n%s' % res)

        if not ok:
            self.report('continuity errors:\n', bads)
            self.report('%d in total!' % len(bads))
            self.report('continuity errors occurred in these spaces:\n',
                        bad_families)

        return ok

    def test_gradients(self):
        from sfepy.discrete.fem.mappings import VolumeMapping

        ok = True
        orders = {'2_3' : 3, '2_4' : 3, '3_4' : 4, '3_8' : 3}

        bads = []
        bad_families = set()
        for (geom, poly_space_base, qp_weights, mesh, ir, ic,
             field, ps, rrc, rcell, crc, ccell, vec,
             edofs, fdofs) in _gen_common_data(orders, self.gels, self.report):
            gel = self.gels[geom]
            conn = mesh.get_conn(gel.name)

            geo_ps = field.gel.poly_space
            rmapping = VolumeMapping(mesh.coors, conn[rcell:rcell+1],
                                     poly_space=geo_ps)
            rori = field.ori[:1] if field.ori is not None else None
            rvg = rmapping.get_mapping(rrc, qp_weights,
                                       poly_space=ps, ori=rori)
            rbfg = rvg.bfg

            cmapping = VolumeMapping(mesh.coors, conn[ccell:ccell+1],
                                     poly_space=geo_ps)
            cori = field.ori[1:] if field.ori is not None else None
            cvg = cmapping.get_mapping(crc, qp_weights,
                                       poly_space=ps, ori=cori)
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
                    bads.append([geom, poly_space_base, ir, ic, ip])
                    bad_families.add((geom, poly_space_base))

                ok = ok and _ok

            self.report('results (dofs, status: 1 ok, 0 failure):\n%s' % res)

        if not ok:
            self.report('gradient continuity errors:\n', bads)
            self.report('%d in total!' % len(bads))
            self.report('gradient continuity errors occurred in these'
                        ' spaces:\n', bad_families)

        return ok

    def test_hessians(self):
        """
        Test the second partial derivatives of basis functions using finite
        differences.
        """
        from sfepy.linalg import combine
        from sfepy.discrete import Integral, PolySpace

        ok = True
        orders = {'2_3' : 3, '2_4' : 3, '3_4' : 4, '3_8' : 3}
        bases = ([ii for ii in combine([['2_3', '2_4', '3_4', '3_8'],
                                        ['lagrange']])])

        for geom, poly_space_base in bases:
            self.report('geometry: %s, base: %s' % (geom, poly_space_base))
            order = orders[geom]

            integral = Integral('i', order=order)
            coors, _ = integral.get_qp(geom)

            ps = PolySpace.any_from_args('ps', self.gels[geom], order,
                                         base=poly_space_base)

            dim = coors.shape[1]
            h1 = nm.zeros((coors.shape[0], dim, dim, ps.n_nod), nm.float64)
            eps = 1e-8
            for ir in range(dim):
                cc = coors.copy()

                cc[:, ir] -= eps
                aux0 = ps.eval_base(cc, diff=1)

                cc[:, ir] += 2 * eps
                aux1 = ps.eval_base(cc, diff=1)

                h1[:, :, ir, :] = 0.5 * (aux1 - aux0) / eps

            h2 = ps.eval_base(coors, diff=2)

            _ok = nm.allclose(h1, h2, rtol=0, atol=50*eps)
            self.report('hessians: error: %.2e ok: %s'
                        % (nm.abs(h1 - h2).max(), _ok))
            ok = ok and _ok

        return ok
