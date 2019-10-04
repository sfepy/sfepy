import numpy as nm


from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.mesh.mesh_generators import gen_block_mesh


from sfepy.discrete.dg.dg_field import DGField


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


def test_get_bc_facet_values():
    # TODO write tests for this for sure
    assert False


def test_set_cell_dofs():
    # TODO use quadratic function for exact projection
    assert False


def test_get_nbrhd_dofs():
    # TODO use with primitive DOFs
    assert False


def test_get_dofs_in_region():
    # TODO test on primitive DOFs
    assert False


def test_get_facet_vols():
    # TODO should be easy
    assert False


def test_get_cell_normals_per_facet():
    # TODO do for whole region and at least one boundary
    assert False

def test_get_both_facet_base_vals():
    #  probably not needed
    assert True

def get_both_facet_state_vals():
    assert False

def get_facet_neighbor_idx():
    # TODO definitely test this
    assert False

def test_get_facet_base():
    # TODO test for simplex and tensor prod mesh
    assert False