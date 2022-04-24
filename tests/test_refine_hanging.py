"""
Test continuity along a boundary with hanging nodes due to a mesh refinement.
"""
import os.path as op

import numpy as nm
import pytest

from sfepy.base.base import assert_, Struct
import sfepy.base.testing as tst

def eval_fun(ts, coors, mode, **kwargs):
    val = nm.sin(nm.sum(coors**2, axis=1) * nm.pi)

    val = val.reshape((coors.shape[0], 1, 1))
    return val

def _gen_lines_2_4(bbox, eps):
    from sfepy.discrete.probes import LineProbe

    LineProbe.__base__.cache = Struct(name='probe_shared_evaluate_cache')

    n_point = 100

    line0 = LineProbe([bbox[0, 0], -eps], [bbox[1, 0], -eps], n_point)
    line1 = LineProbe([bbox[0, 0], eps], [bbox[1, 0], eps], n_point)

    yield line0, line1

    line0 = LineProbe([-eps, bbox[0, 1]], [-eps, bbox[1, 1]], n_point)
    line1 = LineProbe([eps, bbox[0, 1]], [eps, bbox[1, 1]], n_point)

    yield line0, line1

def _gen_grid_3_8(bbox, eps):
    from sfepy.discrete.probes import PointsProbe

    PointsProbe.__base__.cache = Struct(name='probe_shared_evaluate_cache')

    n_point = 100

    px = nm.linspace(bbox[0, 0], bbox[1, 0], n_point)
    py = nm.linspace(bbox[0, 1], bbox[1, 1], n_point)
    pz = nm.linspace(bbox[0, 2], bbox[1, 2], n_point)

    # Test face substitutions.
    pps = nm.meshgrid(px, py, -eps)
    points0 = nm.array([ii.ravel() for ii in pps]).T
    pps = nm.meshgrid(px, py, eps)
    points1 = nm.array([ii.ravel() for ii in pps]).T

    grid0 = PointsProbe(points0)
    grid1 = PointsProbe(points1)

    yield grid0, grid1

    pps = nm.meshgrid(px, -eps, pz)
    points0 = nm.array([ii.ravel() for ii in pps]).T
    pps = nm.meshgrid(px, eps, pz)
    points1 = nm.array([ii.ravel() for ii in pps]).T

    grid0 = PointsProbe(points0)
    grid1 = PointsProbe(points1)

    yield grid0, grid1

    pps = nm.meshgrid(-eps, py, pz)
    points0 = nm.array([ii.ravel() for ii in pps]).T
    pps = nm.meshgrid(eps, py, pz)
    points1 = nm.array([ii.ravel() for ii in pps]).T

    grid0 = PointsProbe(points0)
    grid1 = PointsProbe(points1)

    yield grid0, grid1

    # Test edge substitutions.
    pps = nm.meshgrid(px, -eps, -eps)
    points0 = nm.array([ii.ravel() for ii in pps]).T
    pps = nm.meshgrid(px, eps, eps)
    points1 = nm.array([ii.ravel() for ii in pps]).T

    grid0 = PointsProbe(points0)
    grid1 = PointsProbe(points1)

    yield grid0, grid1

    pps = nm.meshgrid(-eps, py, -eps)
    points0 = nm.array([ii.ravel() for ii in pps]).T
    pps = nm.meshgrid(eps, py, eps)
    points1 = nm.array([ii.ravel() for ii in pps]).T

    grid0 = PointsProbe(points0)
    grid1 = PointsProbe(points1)

    yield grid0, grid1

    pps = nm.meshgrid(-eps, -eps, pz)
    points0 = nm.array([ii.ravel() for ii in pps]).T
    pps = nm.meshgrid(eps, eps, pz)
    points1 = nm.array([ii.ravel() for ii in pps]).T

    grid0 = PointsProbe(points0)
    grid1 = PointsProbe(points1)

    yield grid0, grid1

def _build_filenames(output_dir, key, order, ip):
    return [op.join(output_dir,
                   'test_refine_hanging_lin_%s_%d_%d.vtk' % (key, order, ip)),
            op.join(output_dir,
                    'test_refine_hanging_%s_%d_%d.vtk' % (key, order, ip))]

@pytest.fixture(scope='module')
def gels():
    from sfepy.discrete.fem.geometry_element import GeometryElement

    gels = {}
    for key in ['2_4', '3_8']:
        gel = GeometryElement(key)
        gel.create_surface_facet()
        gels[key] = gel

    return gels

@pytest.mark.slow
def test_continuity(gels, output_dir):
    from sfepy.base.base import Struct
    from sfepy.mesh.mesh_generators import gen_block_mesh
    from sfepy.discrete import FieldVariable
    from sfepy.discrete.fem import FEDomain, Field
    from sfepy.discrete.projections import make_l2_projection_data
    import sfepy.discrete.fem.refine_hanging as rh

    dims = [1.5, 2.0, 1.3]
    shape = [3, 3, 3]
    centre = [0.0, 0.0, 0.0]
    probe_gens = {'2_4' : _gen_lines_2_4, '3_8' : _gen_grid_3_8}

    ok = True
    for key, gel in gels.items():
        probe_gen = probe_gens[key]

        perms = gel.get_conn_permutations()

        dim = gel.dim
        for io, order in enumerate(range(1, 4)):
            mesh00 = gen_block_mesh(dims[:dim], shape[:dim], centre[:dim],
                                    name='block')

            for ip, perm in enumerate(perms):
                tst.report('geometry: %s, order: %d, permutation: %d: %s'
                           % (key, order, ip, perm))
                mesh0 = mesh00.copy()
                conn = mesh0.cmesh.get_conn(dim, 0).indices
                conn = conn.reshape((mesh0.n_el, -1))
                conn[-1, :] = conn[-1, perm]

                domain0 = FEDomain('d', mesh0)

                refine = nm.zeros(mesh0.n_el, dtype=nm.uint8)
                refine[:-1] = 1

                subs = None
                domain, subs = rh.refine(domain0, refine, subs=subs)

                omega = domain.create_region('Omega', 'all')
                field = Field.from_args('fu', nm.float64, 1, omega,
                                        approx_order=order)

                field.substitute_dofs(subs)

                uvar = FieldVariable('u', 'parameter', field,
                                     primary_var_name='(set-to-None)')

                field.restore_dofs(store=True)
                field.substitute_dofs(subs=None, restore=True)

                make_l2_projection_data(uvar, eval_fun)

                field.restore_dofs()

                bbox = domain.get_mesh_bounding_box()
                eps = 1e-7

                save = False
                for ii, (probe0, probe1) in enumerate(probe_gen(bbox, eps)):
                    probe0.set_options(close_limit=0.0)
                    probe1.set_options(close_limit=0.0)

                    pars0, vals0 = probe0(uvar)
                    pars1, vals1 = probe1(uvar)

                    assert_(nm.allclose(pars0, pars1, atol=1e-14, rtol=0.0))

                    _ok = nm.allclose(vals0, vals1, atol=20.0 * eps,
                                      rtol=0.0)
                    if not _ok:
                        save = True
                        tst.report('probe %d failed! (max. error: %e)'
                                   % (ii, nm.abs(vals0 - vals1).max()))

                    ok = ok and _ok

                if (ip == 0) or save:
                    out = uvar.create_output()
                    filenames = _build_filenames(output_dir,
                                                 key, order, ip)

                    domain.mesh.write(filenames[0], out=out)

                    linearization = Struct(kind='adaptive',
                                           min_level=0,
                                           max_level=4,
                                           eps=1e-2)

                    out = uvar.create_output(linearization=linearization)
                    val = out['u']
                    mesh = val.get('mesh', domain.mesh)
                    mesh.write(filenames[1], out=out)

    assert_(ok)

def test_preserve_coarse_entities(output_dir):
    from sfepy.mesh.mesh_generators import gen_block_mesh
    from sfepy.discrete.fem import FEDomain
    import sfepy.discrete.fem.refine_hanging as rh

    dims = [1.5, 2.0]
    shape = [11, 11]
    centre = [0.0, 0.0]

    mesh0 = gen_block_mesh(dims, shape, centre, name='block')
    domain0 = FEDomain('d', mesh0)
    reg = domain0.create_region('surface', 'vertices of surface', 'facet',
                                add_to_regions=False)
    cmesh0 = mesh0.cmesh

    cmesh0.vertex_groups[5::11] = 2
    cmesh0.vertex_groups[reg.vertices] = 1
    cmesh0.cell_groups[0] = 1
    cmesh0.cell_groups[50:60] = 2
    mesh0.write(op.join(output_dir,
                        'test_refine_hanging_ids0.vtk'), io='auto')

    refine = nm.zeros(mesh0.n_el, dtype=nm.uint8)
    refine[0:10] = 1
    refine[5::10] = 1

    domain, _, sub_cells = rh.refine(domain0, refine, subs=None,
                                     ret_sub_cells=True)
    domain.mesh.write(op.join(output_dir,
                              'test_refine_hanging_ids1.vtk'), io='auto')
    cmesh1 = domain.mesh.cmesh

    ii = nm.where(refine == 0)[0]
    conn0 = mesh0.get_conn('2_4')
    v0 = conn0[ii]

    conn1 = domain.mesh.get_conn('2_4')
    v1 = conn1[ii]

    ok = (v0 == v1).all()
    tst.report('coarse cells positions preserved:', ok)

    cgs0 = cmesh0.cell_groups[ii]
    cgs1 = cmesh1.cell_groups[ii]

    _ok = (cgs0 == cgs1).all()
    tst.report('coarse cells cell groups preserved:', _ok)
    ok = ok and _ok

    vgs0 = cmesh0.vertex_groups[v0]
    vgs1 = cmesh1.vertex_groups[v1]

    _ok = (vgs0 == vgs1).all()
    tst.report('coarse cells vertex groups preserved:', _ok)
    ok = ok and _ok

    ii = nm.where(refine == 1)[0]
    cgs0 = cmesh0.cell_groups[ii]
    cgs1 = cmesh1.cell_groups[sub_cells[:, 1:]]

    _ok = (cgs0[:, None] == cgs1).all()
    tst.report('refined cells cell groups preserved:', _ok)
    ok = ok and _ok

    assert_(ok)
