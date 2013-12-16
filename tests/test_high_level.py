import time
import os.path as op
import numpy as nm

from sfepy.base.testing import TestCommon

def fix_u_fun(ts, coors, bc=None, problem=None, extra_arg=None):
    return nm.zeros_like(coors)

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        import sfepy
        from sfepy.fem import Mesh, Domain, H1NodalVolumeField
        mesh = Mesh.from_file('meshes/2d/rectangle_tri.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = Domain('domain', mesh)
        dim = domain.shape.dim

        min_x, max_x = domain.get_mesh_bounding_box()[:,0]
        eps = 1e-8 * (max_x - min_x)

        omega = domain.create_region('Omega', 'all')
        gamma1 = domain.create_region('Gamma1',
                                      'vertices in x < %.10f' % (min_x + eps),
                                      'facet')
        gamma2 = domain.create_region('Gamma2',
                                      'vertices in x > %.10f' % (max_x - eps),
                                      'facet')

        field = H1NodalVolumeField('fu', nm.float64, 'vector', omega,
                                   approx_order=2)

        test = Test(conf=conf, options=options, dim=dim,
                    omega=omega, gamma1=gamma1, gamma2=gamma2,
                    field=field)
        return test

    def test_term_evaluation(self):
        from sfepy.fem import Integral, FieldVariable
        from sfepy.terms.terms import Term

        integral = Integral('i', order=3)

        u = FieldVariable('u', 'parameter', self.field,
                          primary_var_name='(set-to-None)')

        term = Term.new('d_volume(u)', integral, self.omega, u=u)
        term *= 10.0

        term.setup()

        vol = term.evaluate()

        self.report('volume: %.8f == 2000.0' % vol)
        ok = nm.allclose(vol, 2000.0, rtol=1e-15, atol=0)

        ## vec = t1.evaluate() # Returns vector.
        ## vec = t1.evaluate(u=u_vec) # Returns the same vector.
        ## mtx = t1.evaluate(diff_var='u') # Returns matrix.
        ## val = t1.evaluate(v=u_vec, u=u_vec) # Forbidden - virtual variable
        ##                                     # cannot have value.

        return ok

    def test_term_arithmetics(self):
        from sfepy.fem import FieldVariable, Integral
        from sfepy.terms.terms import Term

        integral = Integral('i', order=3)

        u = FieldVariable('u', 'parameter', self.field,
                          primary_var_name='(set-to-None)')

        t1 = Term.new('d_volume(u)', integral, self.omega, u=u)
        t2 = Term.new('d_surface(u)', integral, self.gamma1, u=u)

        expr = 2.2j * (t1 * 5.5 - 3j * t2) * 0.25

        ok = len(expr) == 2
        if not ok:
            self.report('wrong expression length!')

        _ok = nm.allclose(expr[0].sign, 3.025j, rtol=1e-15, atol=0)
        if not _ok:
            self.report('wrong sign of the first term!')
        ok = ok and _ok

        _ok = nm.allclose(expr[1].sign, 1.65, rtol=1e-15, atol=0)
        if not _ok:
            self.report('wrong sign of the second term!')
        ok = ok and _ok

        return ok

    def test_variables(self):
        from sfepy.fem import FieldVariable

        u = FieldVariable('u', 'parameter', self.field,
                          primary_var_name='(set-to-None)')

        u.set_constant(1.0)

        vec = u() # Nodal values.

        ok = nm.allclose(vec, 1.0)

        ## print u() 
        ## print u.get_vector() # Coefficient vector w.r.t. the field space basis.
        ## print u(gamma1)
        ## print u.get_vector(gamma2)

        return ok

    def test_solving(self):
        from sfepy.base.base import IndexedStruct
        from sfepy.fem \
             import FieldVariable, Material, ProblemDefinition, \
                    Function, Equation, Equations, Integral
        from sfepy.fem.conditions import Conditions, EssentialBC
        from sfepy.terms import Term
        from sfepy.solvers.ls import ScipyDirect
        from sfepy.solvers.nls import Newton

        u = FieldVariable('u', 'unknown', self.field)
        v = FieldVariable('v', 'test', self.field, primary_var_name='u')

        m = Material('m', lam=1.0, mu=1.0)
        f = Material('f', val=[[0.02], [0.01]])

        bc_fun = Function('fix_u_fun', fix_u_fun,
                          extra_args={'extra_arg' : 'hello'})

        fix_u = EssentialBC('fix_u', self.gamma1, {'u.all' : bc_fun})
        shift_u = EssentialBC('shift_u', self.gamma2, {'u.0' : 0.1})

        integral = Integral('i', order=3)

        t1 = Term.new('dw_lin_elastic_iso(m.lam, m.mu, v, u)',
                      integral, self.omega, m=m, v=v, u=u)

        t2 = Term.new('dw_volume_lvf(f.val, v)', integral, self.omega, f=f, v=v)

        eq = Equation('balance', t1 + t2)
        eqs = Equations([eq])

        ls = ScipyDirect({})

        nls_status = IndexedStruct()
        nls = Newton({}, lin_solver=ls, status=nls_status)

        pb = ProblemDefinition('elasticity', equations=eqs, nls=nls, ls=ls)
        ## pb.save_regions_as_groups('regions')

        pb.time_update(ebcs=Conditions([fix_u, shift_u]))

        state = pb.solve()

        name = op.join(self.options.out_dir, 'test_high_level_solving.vtk')
        pb.save_state(name, state)

        ok = nls_status.condition == 0
        if not ok:
            self.report('solver did not converge!')

        _ok = state.has_ebc()
        if not _ok:
            self.report('EBCs violated!')

        ok = ok and _ok

        return ok
