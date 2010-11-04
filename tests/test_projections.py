import numpy as nm

from sfepy.base.testing import TestCommon, assert_

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        import sfepy
        from sfepy.fem import Mesh, Domain, Field
        mesh = Mesh.from_file('meshes/2d/rectangle_tri.mesh',
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
        assert_(abs(mtx.sum() - 200.0) < 1e-12)

        return True
