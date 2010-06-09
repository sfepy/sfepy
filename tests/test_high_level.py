import time
import os.path as op
import numpy as nm

from sfepy.base.testing import TestCommon, assert_
from sfepy.base.base import pause, debug

def fix_u_fun(ts, coors, bc=None, extra_arg=None):
    return nm.zeros_like(coors)

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    @TestCommon.xfail
    def test_solving(self):
        import sfepy
        from sfepy.fem \
             import Mesh, Domain, Field, FieldVariable, Material, Function
        from sfepy.fem.conditions import EssentialBC
        from sfepy.terms.terms import Term

        ok = True

        mesh = Mesh.from_file('meshes/2d/special/square_triquad.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = Domain('domain', mesh)
        dim = domain.shape.dim

        min_x, max_x = domain.get_mesh_bounding_box()[:,0]
        eps = 1e-8 * (max_x - min_x)

        omega = domain.create_region('Omega', 'all')
        gamma1 = domain.create_region('Gamma1',
                                      'nodes in x < %.10f' % (min_x + eps))
        gamma2 = domain.create_region('Gamma2',
                                      'nodes in x > %.10f' % (max_x - eps))

        print domain

        field = Field('fu', nm.float64, 'scalar', omega,
                      space='H1', poly_space_base='lagrange', approx_order=2)

        print field

        u = FieldVariable('u', 'unknown', field, dim)
        v = FieldVariable('v', 'test', field, dim, primary_var_name='u')

        print u
        print v

        m = Material('m', lam=1.0, mu=1.0)
        f = Material('f', val=1.0)

        print m
        print f

        bc_fun = Function('fix_u_fun', fix_u_fun,
                          extra_args={'etra_arg' : 'hello'})

        print bc_fun

        fix_u = EssentialBC('fix_u', gamma1, {'u.all' : bc_fun})
        shift_u = EssentialBC('shift_u', gamma2, {'u.0' : 0.1})

        print fix_u
        print shift_u

        # Does not work below...

        t1 = Term('dw_lin_elastic', 'auto', omega, m, v, u)

        print t1
        print t1 - (3 * t1 - 2 * t1) * 4 + t1 * 2.4j

        u.set_nodal_values(1.0)
        u_vec = u.get_vector() # coefficient vector w.r.t. the field space basis.
        vec = t1.evaluate() # Returns vector.
        vec = t1.evaluate(u=u_vec) # Returns the same vector.
        mtx = t1.evaluate(diff_var='u') # Returns matrix.
        val = t1.evaluate(v=u_vec, u=u_vec) # Forbidden - virtual variable
                                            # cannot have value.

        pause()

        eq = Equation(t1 + Term('lin_volume_force', 1, gamma2, f, v))

        eqs = Equations([eq], [fix_u])

        print u() # Nodal values.
        print u.get_vector() # Coefficient vector w.r.t. the field space basis.
        print u(gamma1)
        print u.get_vector(gamma2)

        pb = ProblemDefinition('problem', eqs, nls, ls, ts=None)
        pb.time_update(ebcs=[fix_u, shift_u])
        pb.solve()
        pb.save_solution('results.vtk')

        return ok
