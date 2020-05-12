import numpy as nm
import numpy.testing as nmts

from base.base import Struct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.mesh.mesh_generators import gen_block_mesh

from sfepy.discrete.dg.dg_poly_spaces import get_n_el_nod


from sfepy.discrete.dg.fields import DGField


def prepare_field(approx_order, mesh):

    domain = FEDomain("test_domain", mesh)
    omega = domain.create_region('Omega', 'all')
    regions = {}
    if mesh.dim > 1:
        left = domain.create_region('left',
                                    'vertices in x == 0',
                                    'edge')

        right = domain.create_region('right',
                                     'vertices in x == 1',
                                     'edge')
        bottom = domain.create_region('bottom',
                                      'vertices in y == 0',
                                      'edge')
        top = domain.create_region('top',
                                   'vertices in y == 1',
                                   'edge')
        regions.update({"top": top, "bottom": bottom})
    else:
        left = domain.create_region('left',
                                    'vertices in x == 0',
                                    'vertex')
        right = domain.create_region('right',
                                     'vertices in x == 1',
                                     'vertex')

    regions.update({"left": left, "right": right, "omega" : omega})

    field = DGField('dgfu', nm.float64, 'scalar', omega,
                    approx_order=approx_order)

    return field, regions


def prepare_field_1D(approx_order):
    mesh = gen_block_mesh((1,), (4,), (.5,))
    return prepare_field(approx_order, mesh), mesh


def prepare_field_2D(approx_order):
    mesh = gen_block_mesh((1, 1), (4, 4), (.5, .5))
    return prepare_field(approx_order, mesh), mesh


def test_get_facet_idx1D():
    mesh = gen_block_mesh((1,), (4,), (.5,))
    field, regions = prepare_field(1, mesh)
    assert nm.all(field.get_bc_facet_idx(regions["left"]) == nm.array([[0, 0]]))
    assert nm.all(field.get_bc_facet_idx(regions["right"]) == nm.array([[2, 1]]))


def test_get_facet_idx2D():
    mesh = gen_block_mesh((1, 1), (4, 4), (.5, .5))
    field, regions = prepare_field(1, mesh)
    assert nm.all(field.get_bc_facet_idx(regions["left"]) == nm.array([[0, 3], [1, 3], [2, 3]]))
    assert nm.all(field.get_bc_facet_idx(regions["top"]) == nm.array([[2, 2], [5, 2], [8, 2]]))


def test_create_output2D():
    mesh = gen_block_mesh((1, 1), (4, 4), (.5, .5))
    approx_order = 2
    n_cell_nod = 6
    field, regions = prepare_field(approx_order, mesh)
    dofs = nm.zeros((n_cell_nod * 9, 1))
    output = field.create_output(dofs, "")
    assert output["u_modal_cell_nodes"].mode == "cell_nodes"
    assert nm.allclose(output["u_modal_cell_nodes"].data, nm.zeros((9, n_cell_nod)))
    assert output["u_modal_cell_nodes"].interpolation_scheme is not None


def test_create_output1D():
    mesh = gen_block_mesh((1,), (4,), (.5,))
    approx_order = 2
    n_cell_nod = approx_order + 1
    field, regions = prepare_field(approx_order, mesh)
    dofs = nm.zeros((n_cell_nod * 3, 1))
    output = field.create_output(dofs, "")
    for i in range(n_cell_nod):
        assert output["u_modal{}".format(i)].mode == "cell"
        assert nm.allclose(output["u_modal{}".format(i)].data, nm.zeros((3, 1)))


def test_get_bc_facet_values_1D():
    # TODO write tests for this for sure
    (field, regions), mesh = prepare_field_1D(3)
    fun = 42
    coor, val = field.get_bc_facet_values(fun, regions["left"], ret_coors=True)
    nmts.assert_equal(nm.array(0).reshape((1, 1, 1)), coor)
    nmts.assert_array_equal(val, nm.array(fun).reshape((1, 1)))

def test_get_bc_facet_values_2D_const():
    (field, regions), mesh = prepare_field_2D(2)
    fun = 42
    coor, val = field.get_bc_facet_values(fun, regions["left"], ret_coors=True)
    qps = nm.array([
        [0, 1/2 * (-nm.sqrt(3/5) + 1)],
        [0, 1/2],
        [0, 1/2 * (nm.sqrt(3/5) + 1)]])
    physical_qps = nm.stack([1/3 * qps,
                             1/3 * qps + 1/3,
                             1/3 * qps + 2/3])
    physical_qps[:, :, 0] = 0
    nmts.assert_allclose(coor, physical_qps)
    nmts.assert_allclose(val,  nm.ones((3, 3)) * fun)

def test_get_bc_facet_values_2D():
    (field, regions), mesh = prepare_field_2D(2)
    fun = lambda x: nm.sum(x ** 2, axis=-1)
    coor, val = field.get_bc_facet_values(fun, regions["left"], ret_coors=True)
    qps = nm.array([
        [0, 1 / 2 * (-nm.sqrt(3 / 5) + 1)],
        [0, 1 / 2],
        [0, 1 / 2 * (nm.sqrt(3 / 5) + 1)]])
    physical_qps = nm.stack([1 / 3 * qps,
                             1 / 3 * qps + 1 / 3,
                             1 / 3 * qps + 2 / 3])
    physical_qps[:, :, 0] = 0
    nmts.assert_allclose(coor, physical_qps)
    nmts.assert_allclose(val, fun(physical_qps))


def test_set_dofs_1D():
    (field, regions), mesh = prepare_field_1D(2)
    fun = lambda x: nm.sum(x ** 2, axis=-1)[..., None]
    nods, vals = field.set_dofs(fun, regions["omega"])
    rnods = nm.arange(3*3)
    qps = nm.array([
        [0, 1 / 2 * (-nm.sqrt(3 / 5) + 1)],
        [0, 1 / 2],
        [0, 1 / 2 * (nm.sqrt(3 / 5) + 1)]])
    physical_qps = nm.stack([1 / 3 * qps,
                             1 / 3 * qps + 1 / 3,
                             1 / 3 * qps + 2 / 3])

    physical_qps = physical_qps[:, :, 1]
    nmts.assert_equal(nods, rnods)
    assert vals.shape == physical_qps.shape



def test_set_dofs_2D():
    (field, regions), mesh = prepare_field_2D(3)

    fun = lambda x: nm.sum(x ** 2, axis=-1)[..., None]

    nods, vals = field.set_dofs(fun, regions["omega"])
    n_el_nod = get_n_el_nod(3, 2)
    rnods = nm.arange(n_el_nod * 9)
    rvals_shape = (9, n_el_nod)
    nmts.assert_equal(nods, rnods)
    assert vals.shape == rvals_shape


# def test_get_cell_normals_per_facet():
#     # TODO do for whole region and at least one boundary
#     assert False


# def test_get_both_facet_base_vals():
#     #  TODO test this, used in
#     assert True


# def test_get_both_facet_state_vals():
#     assert False


def test_get_facet_neighbor_idx_1d():
    (field, regions), mesh = prepare_field_1D(3)
    eq_map = Struct()
    eq_map.n_epbc = 0
    eq_map.dg_epbc = []
    nbr_idx = field.get_facet_neighbor_idx(regions["omega"], eq_map)

    rnbr_idx = [[[-1, -1], [1, 0]],  # 0
                [[ 0,  1], [2, 0]],  # 1
                [[ 1, 1],  [-1, -1]],  # 2
                ]
    rnbr_idx = nm.array(rnbr_idx, dtype=nm.int32)
    nmts.assert_equal(rnbr_idx, nbr_idx)


    # periodic BCs
    eq_map.dg_epbc = [(nm.array([[0, 0]], dtype=nm.int32),
                       nm.array([[2, 1]], dtype=nm.int32))]
    field.clear_facet_neighbour_idx_cache()
    nbr_idx = field.get_facet_neighbor_idx(regions["omega"], eq_map)

    rnbr_idx = [[[2, 1], [1, 0]],  # 0
                [[0, 1], [2, 0]],  # 1
                [[1, 1], [0, 0]],  # 2
                ]
    rnbr_idx = nm.array(rnbr_idx, dtype=nm.int32)
    nmts.assert_equal(rnbr_idx, nbr_idx)


def test_get_facet_neighbor_idx_2d():
    (field, regions), mesh = prepare_field_2D(3)
    eq_map = Struct()
    eq_map.n_epbc = 0
    eq_map.dg_epbc = []
    nbr_idx = field.get_facet_neighbor_idx(regions["omega"], eq_map)

    rnbr_idx = [[[-1, -1], [3, 3], [1, 0], [-1, -1]],  # 0
                [[ 0,  2], [4, 3], [2, 0], [-1, -1]],  # 1
                [[ 1, 2],  [5, 3],[-1, -1],[-1, -1]],  # 2
                [[-1, -1], [6, 3],[ 4,  0],[ 0,  1]],  # 3
                [[ 3,  2], [7, 3],[ 5,  0], [1,  1]],  # 4
                [[ 4, 2],  [8, 3], [-1, -1], [2, 1]],  # 5
                [[-1, -1], [-1, -1], [7, 0], [3, 1]],  # 6
                [ [6, 2],  [-1, -1], [8, 0], [4, 1]],  # 7
                [[7, 2], [-1, -1], [-1, -1], [5, 1]],  # 8
                ]
    rnbr_idx = nm.array(rnbr_idx, dtype=nm.int32)
    nmts.assert_equal(rnbr_idx, nbr_idx)

    # periodic BCs
    eq_map.dg_epbc = [(nm.array([[0, 3], [1, 3], [2, 3]], dtype=nm.int32),
                       nm.array([[6, 1], [7, 1], [8, 1]], dtype=nm.int32))]
    field.clear_facet_neighbour_idx_cache()
    nbr_idx = field.get_facet_neighbor_idx(regions["omega"], eq_map)

    rnbr_idx = [[[-1, -1], [3, 3], [1, 0], [6, 1]],  # 0
                [[0, 2], [4, 3], [2, 0], [7, 1]],  # 1
                [[1, 2], [5, 3], [-1, -1], [8, 1]],  # 2
                [[-1, -1], [6, 3], [4, 0], [0, 1]],  # 3
                [[3, 2], [7, 3], [5, 0], [1, 1]],  # 4
                [[4, 2], [8, 3], [-1, -1], [2, 1]],  # 5
                [[-1, -1], [0, 3], [7, 0], [3, 1]],  # 6
                [[6, 2], [1, 3], [8, 0], [4, 1]],  # 7
                [[7, 2], [2, 3], [-1, -1], [5, 1]],  # 8
                ]
    rnbr_idx = nm.array(rnbr_idx, dtype=nm.int32)
    nmts.assert_equal(rnbr_idx, nbr_idx)


# def test_get_facet_base():
#     # TODO test for simplex and tensor prod mesh
#     assert False