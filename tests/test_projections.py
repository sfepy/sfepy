import os.path as op
import numpy as nm

import sfepy
from sfepy.fem import Mesh, Domain, Field, FieldVariable

from sfepy.base.testing import TestCommon, assert_

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        mesh = Mesh.from_file('meshes/2d/square_unit_tri.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = Domain('domain', mesh)

        omega = domain.create_region('Omega', 'all')

        field = Field('linear', nm.float64, 'scalar', omega,
                      space='H1', poly_space_base='lagrange', approx_order=1)

        test = Test(conf=conf, options=options, omega=omega, field=field)
        return test

    def test_mass_matrix(self):
        from sfepy.fem.projections import create_mass_matrix

        field = self.field

        mtx = create_mass_matrix(field)

        assert_(mtx.shape == (field.n_nod, field.n_nod))
        assert_(abs(mtx.sum() - 1.0) < 1e-14)

        return True

    def test_projection_tri_quad(self):
        from sfepy.fem.projections import make_l2_projection

        source = FieldVariable('us', 'unknown', self.field, 1)

        coors = self.field.get_coor()
        vals = nm.sin(2.0 * nm.pi * coors[:,0] * coors[:,1])
        source.data_from_any(vals)

        name = op.join(self.options.out_dir,
                       'test_projection_tri_quad_source.vtk')
        source.save_as_mesh(name)

        mesh = Mesh.from_file('meshes/2d/square_quad.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = Domain('domain', mesh)

        omega = domain.create_region('Omega', 'all')


        field = Field('bilinear', nm.float64, 'scalar', omega,
                      space='H1', poly_space_base='lagrange', approx_order=1)

        target = FieldVariable('ut', 'unknown', field, 1)

        make_l2_projection(target, source)

        name = op.join(self.options.out_dir,
                       'test_projection_tri_quad_target.vtk')
        target.save_as_mesh(name)

        bbox = self.field.domain.get_mesh_bounding_box()
        x = nm.linspace(bbox[0, 0] + 0.001, bbox[1, 0] - 0.001, 20)
        y = nm.linspace(bbox[0, 1] + 0.001, bbox[1, 1] - 0.001, 20)

        xx, yy = nm.meshgrid(x, y)
        test_coors = nm.c_[xx.ravel(), yy.ravel()].copy()

        vec1 = source.evaluate_at(test_coors)
        vec2 = target.evaluate_at(test_coors)

        ok = (nm.abs(vec1 - vec2) < 0.14).all()

        return ok
