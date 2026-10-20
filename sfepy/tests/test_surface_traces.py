import numpy as nm
import pytest

from sfepy.discrete import FieldVariable, Integral
from sfepy.discrete.common.mappings import get_physical_qps
from sfepy.discrete.fem import FEDomain, Field, Mesh
from sfepy.discrete.fem.geometry_element import GeometryElement
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.mesh.mesh_tools import triangulate


def _create_fields(geometry, order, seed):
    dim = int(geometry[0])
    mesh = gen_block_mesh([2] * dim, [3] * dim, [0] * dim,
                          verbose=False)
    if geometry in ('2_3', '3_4'):
        mesh = triangulate(mesh)

    rng = nm.random.default_rng(seed)
    labels = rng.permutation(mesh.n_nod)
    coors = mesh.coors[labels]
    conn = nm.argsort(labels)[mesh.get_conn(geometry)]
    permutations = GeometryElement(geometry).get_conn_permutations()
    for cell in conn:
        cell[:] = cell[permutations[rng.integers(len(permutations))]]
    groups = nm.where(coors[conn].mean(axis=1)[:, 0] < 0, 1, 2)
    mesh = Mesh.from_data('split', coors, None, [conn],
                          [groups.astype(nm.int32)], [geometry])
    domain = FEDomain('domain', mesh)
    omega1 = domain.create_region('Omega1', 'cells of group 1')
    omega2 = domain.create_region('Omega2', 'cells of group 2')
    domain.create_region('Interface', 'r.Omega1 *v r.Omega2', 'facet')
    regions = [domain.create_region('Side%d' % side, 'copy r.Interface',
                                    'facet', 'Omega%d' % side)
               for side in (1, 2)]
    regions[0].setup_mirror_region('Side2')
    fields = [Field.from_args('f%d' % side, nm.float64, 'scalar', omega,
                              approx_order=order)
              for side, omega in enumerate((omega1, omega2), 1)]
    return fields, regions


@pytest.mark.parametrize('geometry', ['2_3', '2_4', '3_4', '3_8'])
@pytest.mark.parametrize('order', [1, 2, 3, 4])
@pytest.mark.parametrize('seed', [0, 17])
def test_surface_trace_nodes(geometry, order, seed):
    """A trace must match physical nodes, including edge and face nodes."""
    fields, regions = _create_fields(geometry, order, seed)
    for field, region, mirror in zip(fields, regions, regions[::-1]):
        field.setup_surface_data(region, mirror.name)

    surfaces = [field.extra_data['sd_' + region.name]
                for field, region in zip(fields, regions)]
    for ii, (field, surface) in enumerate(zip(fields, surfaces)):
        other = 1 - ii
        mirror = surfaces[other]
        traced = surface.get_connectivity(trace_region=regions[other].name)
        expected = fields[other].coors[mirror.econn]
        nm.testing.assert_allclose(field.coors[traced], expected, atol=1e-14,
                                   rtol=0)
        nm.testing.assert_array_equal(nm.sort(traced, axis=1),
                                      nm.sort(surface.econn, axis=1))


@pytest.mark.parametrize('geometry', ['2_3', '2_4', '3_4', '3_8'])
@pytest.mark.parametrize('order', [1, 2, 3, 4])
@pytest.mark.parametrize('seed', [0, 17])
def test_surface_values(geometry, order, seed):
    """Facet interpolation must reproduce an affine physical function."""
    fields, regions = _create_fields(geometry, order, seed)
    integral = Integral('i', order=2 * order)
    for field, region in zip(fields, regions):
        field.domain.create_surface_group(region)
        field.setup_surface_data(region)
        variable = FieldVariable('u', 'parameter', field,
                                 primary_var_name='(set-to-None)')
        coefficients = nm.arange(1, field.coors.shape[1] + 1)
        variable.set_data(1 + field.coors @ coefficients)
        values = variable.evaluate('val', region=region, integral=integral,
                                   integration='facet')
        coors = get_physical_qps(region, integral).values
        expected = 1 + coors @ coefficients
        nm.testing.assert_allclose(values.ravel(), expected, atol=1e-13,
                                   rtol=0)
