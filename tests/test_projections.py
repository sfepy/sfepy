from __future__ import absolute_import
import os.path as op
import numpy as nm

import sfepy
from sfepy.discrete import FieldVariable
from sfepy.discrete.fem import Mesh, FEDomain, Field

from sfepy.base.base import assert_
from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        mesh = Mesh.from_file('meshes/2d/square_unit_tri.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = FEDomain('domain', mesh)

        omega = domain.create_region('Omega', 'all')

        field = Field.from_args('linear', nm.float64, 'scalar', omega,
                                approx_order=1)

        test = Test(conf=conf, options=options, omega=omega, field=field)
        return test

    def test_mass_matrix(self):
        from sfepy.discrete.projections import create_mass_matrix

        field = self.field

        mtx = create_mass_matrix(field)

        assert_(mtx.shape == (field.n_nod, field.n_nod))
        assert_(abs(mtx.sum() - 1.0) < 1e-14)

        return True

    def test_projection_tri_quad(self):
        from sfepy.discrete.projections import make_l2_projection

        source = FieldVariable('us', 'unknown', self.field)

        coors = self.field.get_coor()
        vals = nm.sin(2.0 * nm.pi * coors[:,0] * coors[:,1])
        source.set_data(vals)

        name = op.join(self.options.out_dir,
                       'test_projection_tri_quad_source.vtk')
        source.save_as_mesh(name)

        mesh = Mesh.from_file('meshes/2d/square_quad.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = FEDomain('domain', mesh)

        omega = domain.create_region('Omega', 'all')


        field = Field.from_args('bilinear', nm.float64, 'scalar', omega,
                                approx_order=1)

        target = FieldVariable('ut', 'unknown', field)

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

        ok = (nm.abs(vec1 - vec2) < 0.01).all()

        return ok

    def test_projection_iga_fem(self):
        try:
            from igakit import igalib
        except ImportError:
            self.report('iga-fem projection not-tested (missing igalib module)!')
            return True

        from sfepy.discrete import FieldVariable
        from sfepy.discrete.fem import FEDomain, Field
        from sfepy.discrete.iga.domain import IGDomain
        from sfepy.mesh.mesh_generators import gen_block_mesh
        from sfepy.discrete.iga.domain_generators import gen_patch_block_domain
        from sfepy.discrete.projections import (make_l2_projection,
                                                make_l2_projection_data)

        shape = [10, 12, 12]
        dims = [5, 6, 6]
        centre = [0, 0, 0]
        degrees = [2, 2, 2]

        nurbs, bmesh, regions = gen_patch_block_domain(dims, shape, centre,
                                                       degrees,
                                                       cp_mode='greville',
                                                       name='iga')
        ig_domain = IGDomain('iga', nurbs, bmesh, regions=regions)

        ig_omega = ig_domain.create_region('Omega', 'all')
        ig_field = Field.from_args('iga', nm.float64, 1, ig_omega,
                                   approx_order='iga', poly_space_base='iga')
        ig_u = FieldVariable('ig_u', 'parameter', ig_field,
                             primary_var_name='(set-to-None)')

        mesh = gen_block_mesh(dims, shape, centre, name='fem')
        fe_domain = FEDomain('fem', mesh)

        fe_omega = fe_domain.create_region('Omega', 'all')
        fe_field = Field.from_args('fem', nm.float64, 1, fe_omega,
                                   approx_order=2)
        fe_u = FieldVariable('fe_u', 'parameter', fe_field,
                             primary_var_name='(set-to-None)')

        def _eval_data(ts, coors, mode, **kwargs):
            return nm.prod(coors**2, axis=1)[:, None, None]

        make_l2_projection_data(ig_u, _eval_data)

        make_l2_projection(fe_u, ig_u) # This calls ig_u.evaluate_at().

        coors = 0.5 * nm.random.rand(20, 3) * dims

        ig_vals = ig_u.evaluate_at(coors)
        fe_vals = fe_u.evaluate_at(coors)

        ok = nm.allclose(ig_vals, fe_vals, rtol=0.0, atol=1e-12)
        if not ok:
            self.report('iga-fem projection failed!')
            self.report('coors:')
            self.report(coors)
            self.report('iga fem diff:')
            self.report(nm.c_[ig_vals, fe_vals, nm.abs(ig_vals - fe_vals)])

        return ok

    def test_project_tensors(self):
        from sfepy.discrete import FieldVariable
        from sfepy.discrete.projections import project_by_component

        ok = True

        u = FieldVariable('u', 'parameter', self.field,
                          primary_var_name='(set-to-None)')
        u.set_constant(1.0)

        component = FieldVariable('component', 'parameter', self.field,
                                  primary_var_name='(set-to-None)')

        nls_options = {'eps_a' : 1e-16, 'i_max' : 1}

        u_qp = u.evaluate()
        u2 = FieldVariable('u2', 'parameter', self.field,
                           primary_var_name='(set-to-None)')
        project_by_component(u2, u_qp, component, self.field.approx_order,
                             nls_options=nls_options)

        _ok = self.compare_vectors(u(), u2())
        ok = ok and _ok

        gu_qp = u.evaluate(mode='grad')

        gfield = Field.from_args('gu', nm.float64, 2, self.field.region,
                                 approx_order=self.field.approx_order)
        gu = FieldVariable('gu', 'parameter', gfield,
                           primary_var_name='(set-to-None)')

        project_by_component(gu, gu_qp, component, gfield.approx_order,
                             nls_options=nls_options)

        _ok = self.compare_vectors(gu(), nm.zeros_like(gu()))
        ok = ok and _ok

        return ok
