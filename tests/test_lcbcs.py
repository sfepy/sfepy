from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_laplace_shifted_periodic(self):
        import numpy as nm
        from sfepy.mesh.mesh_generators import gen_block_mesh
        from sfepy.discrete.fem import FEDomain
        from examples.standalone.interactive.laplace_shifted_periodic import run

        dims = [2.0, 1.0]
        shape = [21, 11]
        centre = [0.0, 0.0]
        mesh = gen_block_mesh(dims, shape, centre, name='block-fem')
        fe_domain = FEDomain('domain', mesh)

        pb, state = run(fe_domain, 3)

        gamma3 = pb.domain.regions['Gamma3']
        gamma4 = pb.domain.regions['Gamma4']

        field = pb.fields['fu']

        # Check that the shift equals to one.
        i3 = field.get_dofs_in_region(gamma3, merge=True)
        i4 = field.get_dofs_in_region(gamma4, merge=True)

        i_corners = nm.array([0, shape[0] - 1])
        ii = nm.setdiff1d(nm.arange(len(i3)), i_corners)

        vals = state()

        shift = vals[i3] - vals[i4]

        ok = (shift[i_corners] == 0.0).all()

        ok = ok and nm.allclose(shift[ii], 1.0, rtol=0.0, atol=1e-14)

        if not ok:
            self.report('wrong shift:', shift)

        return ok
