from __future__ import absolute_import
import numpy as nm

from sfepy.base.testing import TestCommon

expected_normals = { # Need to be normalized!
    '2_3' : nm.array([[ 0, -1],
                      [ 1,  1],
                      [-1,  0]], dtype=nm.float64),
    '2_4' : nm.array([[ 0, -1],
                      [ 1,  0],
                      [ 0,  1],
                      [-1,  0]], dtype=nm.float64),
    '3_4' : nm.array([[ 0,  0, -1],
                      [-1,  0,  0],
                      [ 0, -1,  0],
                      [ 1,  1,  1]], dtype=nm.float64),
    '3_8' : nm.array([[ 0,  0, -1],
                      [-1,  0,  0],
                      [ 0, -1,  0],
                      [ 0,  0,  1],
                      [ 1,  0,  0],
                      [ 0,  1,  0]], dtype=nm.float64),
}

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_normals(self):
        """
        Check orientations of surface normals on the reference elements.
        """
        import sfepy
        from sfepy.discrete import Integral, PolySpace
        from sfepy.discrete.fem import Mesh, FEDomain
        from sfepy.discrete.fem.mappings import SurfaceMapping
        from sfepy.linalg import normalize_vectors

        ok = True

        for geom in ['2_3', '2_4', '3_4', '3_8']:
            mesh = Mesh.from_file('meshes/elements/%s_1.mesh' % geom,
                                  prefix_dir=sfepy.data_dir)
            domain = FEDomain('domain', mesh)
            surface = domain.create_region('Surface', 'vertices of surface',
                                           'facet')
            domain.create_surface_group(surface)

            sd = domain.surface_groups[surface.name]

            coors = domain.get_mesh_coors()
            gel = domain.geom_els[geom].surface_facet
            ps = PolySpace.any_from_args('aux', gel, 1)

            mapping = SurfaceMapping(coors, sd.get_connectivity(), ps)

            integral = Integral('i', order=1)
            vals, weights = integral.get_qp(gel.name)

            # Evaluate just in the first quadrature point...
            geo = mapping.get_mapping(vals[:1], weights[:1])

            expected = expected_normals[geom].copy()
            normalize_vectors(expected)

            _ok = nm.allclose(expected, geo.normal[:, 0, :, 0],
                              rtol=0.0, atol=1e-14)
            self.report('%s: %s' % (geom, _ok))

            if not _ok:
                self.report('expected:')
                self.report(expected)
                self.report('actual:')
                self.report(geo.normal[:, 0, :, 0])

            ok = ok and _ok

        return ok
