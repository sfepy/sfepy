"""
Test continuity along a boundary with hanging nodes due to a mesh refinement.
"""
from __future__ import absolute_import
import os.path as op

import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy.base.base import assert_, Struct

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

    n_point = 10

    px = nm.linspace(bbox[0, 0], bbox[1, 0], n_point)
    py = nm.linspace(bbox[0, 1], bbox[1, 1], n_point)
    pz = nm.linspace(bbox[0, 2], bbox[1, 2], n_point)

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

def _build_filenames(output_dir, key, order, ip):
    return [op.join(output_dir,
                   'test_refine_hanging_lin_%s_%d_%d.vtk' % (key, order, ip)),
            op.join(output_dir,
                    'test_refine_hanging_%s_%d_%d.vtk' % (key, order, ip))]

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete.fem.geometry_element import GeometryElement

        gels = {}
        gel_names = ['2_4', '3_8']
        for key in gel_names:
            gel = GeometryElement(key)
            gel.create_surface_facet()
            gels[key] = gel

        return Test(conf=conf, options=options, gel_names=gel_names, gels=gels)

    def test_continuity(self):
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
        for key in self.gel_names:
            gel = self.gels[key]
            probe_gen = probe_gens[key]

            perms = gel.get_conn_permutations()

            dim = gel.dim
            for io, order in enumerate(range(1, 4)):
                mesh00 = gen_block_mesh(dims[:dim], shape[:dim], centre[:dim],
                                       name='block')

                for ip, perm in enumerate(perms):
                    self.report('geometry: %s, order: %d, permutation: %d: %s'
                                % (key, order, ip, perm))
                    mesh0 = mesh00.copy()
                    conn = mesh0.get_conn(gel.name)
                    conn[-1, :] = conn[0, perm]

                    domain0 = FEDomain('d', mesh0)

                    refine = nm.zeros(mesh0.n_el, dtype=nm.uint8)
                    refine[:-1] = 1

                    gsubs = None
                    domain, gsubs = rh.refine(domain0, refine, gsubs=gsubs)

                    omega = domain.create_region('Omega', 'all')
                    field = Field.from_args('fu', nm.float64, 1, omega,
                                            approx_order=order)

                    basis_transform = rh.eval_basis_transform(field, gsubs)
                    field.substitute_dofs(gsubs)
                    field.set_basis_transform(basis_transform)

                    uvar = FieldVariable('u', 'parameter', field,
                                         primary_var_name='(set-to-None)')

                    make_l2_projection_data(uvar, eval_fun)

                    field.restore_dofs()

                    out = uvar.create_output()
                    filenames = _build_filenames(self.options.out_dir,
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

                    bbox = domain.get_mesh_bounding_box()
                    eps = 1e-7

                    for ii, (probe0, probe1) in enumerate(probe_gen(bbox, eps)):
                        probe0.set_options(close_limit=0.0)
                        probe1.set_options(close_limit=0.0)

                        pars0, vals0 = probe0(uvar)
                        pars1, vals1 = probe1(uvar)

                        assert_(nm.allclose(pars0, pars1, atol=1e-14, rtol=0.0))

                        _ok = nm.allclose(vals0, vals1, atol=10.0 * eps,
                                          rtol=0.0)
                        if not _ok:
                            self.report('probe %d failed! (max. error: %e)'
                                        % (ii, nm.abs(vals0 - vals1).max()))

                        ok = ok and _ok

        return ok
