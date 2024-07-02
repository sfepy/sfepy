import os.path as op
import numpy as nm
import pytest

import sfepy
from sfepy.discrete import FieldVariable
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.base.base import assert_
import sfepy.base.testing as tst

@pytest.fixture(scope='module')
def data():
    from sfepy.base.base import Struct

    mesh = Mesh.from_file('meshes/2d/square_unit_tri.mesh',
                          prefix_dir=sfepy.data_dir)
    domain = FEDomain('domain', mesh)

    omega = domain.create_region('Omega', 'all')

    field = Field.from_args('linear', nm.float64, 'scalar', omega,
                            approx_order=1)
    return Struct(omega=omega, field=field)

def test_mass_matrix(data):
    from sfepy.discrete.projections import create_mass_matrix

    field = data.field

    mtx = create_mass_matrix(field)

    assert_(mtx.shape == (field.n_nod, field.n_nod))
    assert_(abs(mtx.sum() - 1.0) < 1e-14)

def test_projection_tri_quad(data, output_dir):
    from sfepy.discrete.projections import make_l2_projection

    source = FieldVariable('us', 'unknown', data.field)

    coors = data.field.get_coor()
    vals = nm.sin(2.0 * nm.pi * coors[:,0] * coors[:,1])
    source.set_data(vals)

    name = op.join(output_dir,
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

    name = op.join(output_dir,
                   'test_projection_tri_quad_target.vtk')
    target.save_as_mesh(name)

    bbox = data.field.domain.get_mesh_bounding_box()
    x = nm.linspace(bbox[0, 0] + 0.001, bbox[1, 0] - 0.001, 20)
    y = nm.linspace(bbox[0, 1] + 0.001, bbox[1, 1] - 0.001, 20)

    xx, yy = nm.meshgrid(x, y)
    test_coors = nm.c_[xx.ravel(), yy.ravel()].copy()

    vec1 = source.evaluate_at(test_coors)
    vec2 = target.evaluate_at(test_coors)

    ok = (nm.abs(vec1 - vec2) < 0.01).all()
    assert_(ok)

def test_projection_iga_fem():
    try:
        from igakit import igalib; igalib
    except ImportError:
        tst.report('iga-fem projection not-tested (missing igalib module)!')
        return

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
                               approx_order='iga', poly_space_basis='iga')
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
        tst.report('iga-fem projection failed!')
        tst.report('coors:')
        tst.report(coors)
        tst.report('iga fem diff:')
        tst.report(nm.c_[ig_vals, fe_vals, nm.abs(ig_vals - fe_vals)])

    assert_(ok)

def test_project_tensors(data):
    from sfepy.discrete import FieldVariable
    from sfepy.discrete.projections import project_by_component

    ok = True

    u = FieldVariable('u', 'parameter', data.field,
                      primary_var_name='(set-to-None)')
    u.set_constant(1.0)

    component = FieldVariable('component', 'parameter', data.field,
                              primary_var_name='(set-to-None)')

    nls_options = {'eps_a' : 1e-16, 'i_max' : 1}

    u_qp = u.evaluate()
    u2 = FieldVariable('u2', 'parameter', data.field,
                       primary_var_name='(set-to-None)')
    project_by_component(u2, u_qp, component, data.field.approx_order,
                         nls_options=nls_options)

    _ok = tst.compare_vectors(u(), u2())
    ok = ok and _ok

    gu_qp = u.evaluate(mode='grad')

    gfield = Field.from_args('gu', nm.float64, 2, data.field.region,
                             approx_order=data.field.approx_order)
    gu = FieldVariable('gu', 'parameter', gfield,
                       primary_var_name='(set-to-None)')

    project_by_component(gu, gu_qp, component, gfield.approx_order,
                         nls_options=nls_options)

    _ok = tst.compare_vectors(gu(), nm.zeros_like(gu()))
    ok = ok and _ok

    assert_(ok)
