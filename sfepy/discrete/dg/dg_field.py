"""
Fields for Discontinous Galerkin
"""
import numpy as nm

# sfepy imports
from sfepy.base.base import assert_, basestr, Struct
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete import Integral, FieldVariable
from sfepy.discrete.iga.mappings import IGMapping
from sfepy.discrete.iga.iga import get_bezier_element_entities
from six.moves import range
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete.fem.poly_spaces import PolySpace


# local imports
from dg_basis import LegendrePolySpace

class DGField(Field):
    # TODO inherit from FEField?
    # - that would require adding element dofs to it
    family_name = 'volume_H1_DGLegendre'
    is_surface = False


    def __init__(self, name, dtype, shape, region, space="H1",
                 poly_space_base="dglegendre",  approx_order=0):
        """

        :param name:
        :param dtype:
        :param shape:
        :param region:
        :param space:
        :param poly_space_base: use Legendre base
        :param approx_order: 0 for finite volume methods
        """
        shape = parse_shape(shape, region.domain.shape.dim)
        if not self._check_region(region):
            raise ValueError('unsuitable region for field %s! (%s)' %
                             (name, region.name))
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)

        self.domain = region.domain
        self.region = region  # TODO is use of ragion ok?
        self.approx_order = approx_order
        # self._setup_geometry()
        # self._setup_kind()
        self._setup_shape()
        self.space = space
        self.poly_space_base = poly_space_base
        # TODO put LegendrePolySpace to table in PolySpace any_from_args, or use only Legendre for DG?
        self.poly_space = LegendrePolySpace("1_2_H1_dglegendre_1", region.domain.geom_els["1_2"], approx_order)
        # poly_space = PolySpace.any_from_args("legendre", region.domain.geom_els["1_2"], base="legendre", order=approx_order)

        self._setup_global_base()


        # wrapper for convinient integration
        self._integrate = lambda fun: Integral("dg_fi", order=approx_order).integrate(fun, order=approx_order)


    def _check_region(self, region):
        # TODO what are the requirements on region?
        return True


    def _init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        return None

    def _setup_global_base(self):
        self._init_econn()
        self.n_vertex_dof = 0  # in DG we will propably never need vertex DOFs
        self.n_edge_dof = 0 # use facets DOFS for AFS methods
        self.n_face_dof = 0 # use facet DOF for AFS methods

        self.n_element_dof = self.region.get_n_cells(self.is_surface) * (self.approx_order + 1) # is that right?
        self.element_remap = self.region.cells  # what is remap used for?
        self.element_dofs = self.region.cells  # TODO this should be kinda main region

        self.n_nod = self.n_vertex_dof + self.n_edge_dof + self.n_face_dof + self.n_element_dof

    def _setup_shape(self):
        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape

    def get_dofs_in_region(self, region, merge=True):
        """
        Return indices of DOFs that belong to the given region and group.
        """

        # node_desc = self.node_desc # TODO what is node_desc for?


        # DG probably never needs facet DOFs, but AFS...?
        # edofs = nm.empty((0,), dtype=nm.int32)
        # if node_desc.edge is not None:
        #     edofs = self._get_facet_dofs(region.edges,
        #                                  self.edge_remap,
        #                                  self.edge_dofs)
        # dofs.append(edofs)
        #
        # fdofs = nm.empty((0,), dtype=nm.int32)
        # if node_desc.face is not None:
        #     fdofs = self._get_facet_dofs(region.faces,
        #                                  self.face_remap,
        #                                  self.face_dofs)
        # dofs.append(fdofs)

        dofs = []
        eldofs = nm.empty((0,), dtype=nm.int32)
        if region.has_cells(): # TODO use "and (node_desc.element is not None)"
            els = self.element_remap[region.cells]
            eldofs = self.element_dofs[els[els >= 0]].ravel()
        dofs.append(eldofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs


    def create_mapping(self, region, integral, integration):
        # TODO create a new reference mapping, maybe steal this from FE
        raise NotImplemented

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        # TODO used to setup DOFs directly, basically sampleIC from TSs

        mesh = self.domain.mesh  # TODO use remap and self.element_dofs to get indicies!
        sic = nm.zeros((2, self.mesh.n_el, 1), dtype=nm.float64)

        c = (mesh.coors[1:] + mesh.coors[:-1]) / 2  # center
        s = (mesh.coors[1:] - mesh.coors[:-1]) / 2  # scale
        sic[0, :] = self._integrate(lambda t: fun(c + t * s)) / 2
        sic[1, :] = 3 * self._integrate(lambda t: t * fun(c + t * s)) / 2

        vals = sic  # TODO should vals correspond to shape of variable? i.e. n_nod
        nods = c # TODO what should nods be? Mapping from element_dofs? element_remap?

        return nods, vals

# _get_facets
# create_basis_context
# create_eval_mesh
# create_mesh
# create_output
# get_data_shape
# get_dofs_in_region
# get_econn
# get_true_order
# is_higher_order
# setup_extra_data
# is_surface

missing =[ '__dir__',
'__eq__',
'__ge__',
'__gt__',
'__init_subclass__',
'__le__',
'__lt__',
'__ne__',
'_check_region',
'_create_interpolant',
'_eval_basis_transform',
'_get_facet_dofs',
'_init_econn',
'_set_approx_order',
'_setup_bubble_dofs',
'_setup_edge_dofs',
'_setup_esurface',
'_setup_face_dofs',
'_setup_facet_dofs',
'_setup_facet_orientations',
'_setup_geometry',
'_setup_global_base',
'_setup_shape',
'_setup_vertex_dofs',
'_substitute_dofs',
'average_qp_to_vertices',
'basis_transform',
'bf',
'bubble_dofs',
'bubble_remap',
'clear_qp_base',
'coors',
'create_basis_context',
'create_bqp',
'create_mesh',
'create_output',
'domain',
'econn',
'econn0',
'edge_dofs',
'edge_remap',
'efaces',
'extend_dofs',
'face_dofs',
'face_remap',
'force_bubble',
'gel',
'get_base',
'get_connectivity',
'get_coor',
'get_data_shape',
'get_dofs_in_region',
'get_econn',
'get_evaluate_cache',
'get_output_approx_order',
'get_qp',
'get_true_order',
'get_vertices',
'interp_to_qp',
'interp_v_vals_to_n_vals',
'is_higher_order',
'is_surface',
'linearize',
'mappings',
'mappings0',
'n_bubble_dof',
'n_components',
'n_edge_dof',
'n_face_dof',
'n_nod',
'n_vertex_dof',
'node_desc',
'ori',
'point_data',
'qp_coors',
'remove_extra_dofs',
'restore_dofs',
'restore_substituted',
'set_basis_transform',
'set_coors',
'setup_coors',
'setup_extra_data',
'setup_point_data',
'setup_surface_data',
'stored_subs',
'substitute_dofs',
'surface_data',
'unused_dofs',
'val_shape',
'vertex_remap',
'vertex_remap_i' ]

fem_field_dir = set(['__add__',
 '__class__',
 '__delattr__',
 '__dict__',
 '__dir__',
 '__doc__',
 '__eq__',
 '__format__',
 '__ge__',
 '__getattribute__',
 '__gt__',
 '__hash__',
 '__iadd__',
 '__init__',
 '__init_subclass__',
 '__le__',
 '__lt__',
 '__module__',
 '__ne__',
 '__new__',
 '__reduce__',
 '__reduce_ex__',
 '__repr__',
 '__setattr__',
 '__sizeof__',
 '__str__',
 '__subclasshook__',
 '__weakref__',
 '_all',
 '_check_region',
 '_create_interpolant',
 '_eval_basis_transform',
 '_format_sequence',
 '_get_facet_dofs',
 '_init_econn',
 '_set_approx_order',
 '_setup_bubble_dofs',
 '_setup_edge_dofs',
 '_setup_esurface',
 '_setup_face_dofs',
 '_setup_facet_dofs',
 '_setup_facet_orientations',
 '_setup_geometry',
 '_setup_global_base',
 '_setup_kind',
 '_setup_shape',
 '_setup_vertex_dofs',
 '_str',
 '_substitute_dofs',
 'approx_order',
 'average_qp_to_vertices',
 'basis_transform',
 'bf',
 'bubble_dofs',
 'bubble_remap',
 'clear_mappings',
 'clear_qp_base',
 'coors',
 'copy',
 'create_basis_context',
 'create_bqp',
 'create_eval_mesh',
 'create_mapping',
 'create_mesh',
 'create_output',
 'domain',
 'dtype',
 'econn',
 'econn0',
 'edge_dofs',
 'edge_remap',
 'efaces',
 'evaluate_at',
 'extend_dofs',
 'face_dofs',
 'face_remap',
 'family_name',
 'force_bubble',
 'from_args',
 'from_conf',
 'gel',
 'get',
 'get_base',
 'get_connectivity',
 'get_coor',
 'get_data_shape',
 'get_dofs_in_region',
 'get_econn',
 'get_evaluate_cache',
 'get_mapping',
 'get_output_approx_order',
 'get_qp',
 'get_true_order',
 'get_vertices',
 'interp_to_qp',
 'interp_v_vals_to_n_vals',
 'is_higher_order',
 'is_surface',
 'linearize',
 'mappings',
 'mappings0',
 'n_bubble_dof',
 'n_components',
 'n_edge_dof',
 'n_face_dof',
 'n_nod',
 'n_vertex_dof',
 'name',
 'node_desc',
 'ori',
 'point_data',
 'poly_space',
 'poly_space_base',
 'qp_coors',
 'region',
 'remove_extra_dofs',
 'restore_dofs',
 'restore_substituted',
 'save_mappings',
 'set_basis_transform',
 'set_coors',
 'set_default',
 'set_dofs',
 'setup_coors',
 'setup_extra_data',
 'setup_point_data',
 'setup_surface_data',
 'shape',
 'space',
 'stored_subs',
 'str_all',
 'str_class',
 'substitute_dofs',
 'surface_data',
 'to_dict',
 'unused_dofs',
 'update',
 'val_shape',
 'vertex_remap',
 'vertex_remap_i'])


if __name__ == '__main__':
    from toolz import *

    X1 = 0.
    XN1 = 1.
    n_nod = 100
    n_el = n_nod - 1
    coors = nm.linspace(X1, XN1, n_nod).reshape((n_nod, 1))
    conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
    mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
    descs = ['1_2']
    mesh = Mesh.from_data('advection_1d', coors, None,
                          [conn], [mat_ids], descs)

    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')

    dgfield = DGField("dgfu", nm.float64, 'scalar', omega, poly_space_base=None, approx_order=1)

    dg_field_dir = set(dir(dgfield))
    # print("FEM - DG: {}".format(fem_field_dir.difference(dg_field_dir)))
    # print("DG - FEM: {}".format(dg_field_dir.difference(fem_field_dir)))

    u = FieldVariable('u', 'unknown', dgfield, history=1)
    v = FieldVariable('v', 'test', dgfield, primary_var_name='u')

