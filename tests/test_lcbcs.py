from __future__ import absolute_import
import numpy as nm

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_stokes_slip_bc(self):
        import scipy.sparse as sp

        from sfepy.base.conf import ProblemConf
        from sfepy.discrete import Problem
        import examples.navier_stokes.stokes_slip_bc as ssb

        conf = ProblemConf.from_module(ssb)
        pb = Problem.from_conf(conf, init_solvers=False)
        pb.time_update()
        variables = pb.get_variables()

        adi = variables.adi
        lcdi = variables.lcdi
        mtx = variables.mtx_lcbc

        ok = adi.var_names == lcdi.var_names
        self.report('same adi-lcdi ordering:', ok)

        ublock = mtx[adi.indx['u']]
        ir, ic = ublock.nonzero()
        ir += adi.indx['u'].start

        i0, i1 = adi.indx['u'].start, adi.indx['u'].stop
        _ok0 = (i0 <= ir).all() and  (ir < i1).all()
        self.report('u block rows in [%d %d[: %s' % (i0, i1, _ok0))

        i0, i1 = lcdi.indx['u'].start, lcdi.indx['u'].stop
        _ok1 = (i0 <= ic).all() and  (ic < i1).all()
        self.report('u block cols in [%d %d[: %s' % (i0, i1, _ok1))

        ok = ok and _ok0 and _ok1

        pblock = mtx[adi.indx['p']]
        ir, ic, iv = sp.find(pblock)
        ir += adi.indx['p'].start

        i0, i1 = adi.indx['p'].start, adi.indx['p'].stop
        _ok0 = (i0 <= ir).all() and  (ir < i1).all()
        self.report('p block rows in [%d %d[: %s' % (i0, i1, _ok0))

        i0, i1 = lcdi.indx['p'].start, lcdi.indx['p'].stop
        _ok1 = (i0 <= ic).all() and  (ic < i1).all()
        self.report('p block cols in [%d %d[: %s' % (i0, i1, _ok1))

        ok = ok and _ok0 and _ok1

        _ok0 = (len(ir) == adi.n_dof['p'])
        self.report('p block size correct:', _ok0)
        _ok1 = ((ir - adi.indx['p'].start) == (ic - lcdi.indx['p'].start)).all()
        self.report('p block diagonal:', _ok1)
        _ok2 = (iv == 1.0).all()
        self.report('p block identity:', _ok2)

        ok = ok and _ok0 and _ok1 and _ok2

        return ok

    def test_laplace_shifted_periodic(self):
        from sfepy.mesh.mesh_generators import gen_block_mesh
        from sfepy.discrete.fem import FEDomain
        from examples.diffusion.laplace_shifted_periodic import run

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
