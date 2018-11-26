"""
Fields for Discontinous Galerkin
"""
import numpy as nm
import six


# sfepy imports
from sfepy.base.base import assert_, basestr, Struct
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete import Integral, FieldVariable
from six.moves import range
from sfepy.discrete.fem import Mesh, Field
from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.discrete.fem.mappings import VolumeMapping


# local imports
from dg_basis import LegendrePolySpace

class DGField(Field):
    # TODO inherit from FEField?
    # - that would require adding element dofs to it
    # TODO rename to DGVolumeField? does DGSurfaceeField make sense? Use abstrakt class and subclasses?
    family_name = 'volume_H1_DGLegendre'
    is_surface = False

    def __init__(self, name, dtype, shape, region, space="H1",
                 poly_space_base="dglegendre",  approx_order=0):
        """
        Creates DG Field, with Legendre poly space and integral corresponding to
        approx_order + 1.
        :param name:
        :param dtype:
        :param shape: shape of the problem?
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

        # geometry
        self.domain = region.domain
        self.region = region
        self.approx_order = approx_order

        self._setup_geometry()
        self._setup_shape()
        self._setup_all_dofs()

        # approximation space
        self.space = space
        self.poly_space_base = poly_space_base
        # TODO put LegendrePolySpace to table in PolySpace any_from_args, or use only Legendre for DG?
        self.poly_space = LegendrePolySpace("1_2_H1_dglegendre_1", self.gel, approx_order)
        # poly_space = PolySpace.any_from_args("legendre", region.domain.geom_els["1_2"], base="legendre", order=approx_order)

        # integral
        self.clear_qp_base()
        self.integral = Integral("dg_fi", order=approx_order+1)

        # mapping
        self.mappings = {}
        self.mappings0 = {}

        self.ori = None
        self.basis_transform = None



    def _check_region(self, region):
        # TODO what are the requirements on region?
        return True

    def _init_econn(self):
        """
        Initialize the extended DOF connectivity. What is this supposed to do?
        """

        return None

    def _setup_all_dofs(self):
        """
        Sets up all the differet kinds of DOFs, for DG only bubble DOFs
        originaly called _setup_global_base
        """
        self._init_econn()
        self.n_vertex_dof = 0  # in DG we will propably never need vertex DOFs
        self.n_edge_dof = 0 # use facets DOFS for AFS methods
        self.n_face_dof = 0 # use facet DOF for AFS methods

        (self.n_bubble_dof,
        self.bubble_remap,
        self.bubble_dofs) = self._setup_bubble_dofs()

        self.n_nod = self.n_vertex_dof + self.n_edge_dof + self.n_face_dof + self.n_bubble_dof

    def _setup_bubble_dofs(self):
        """
        Creates DOF information for  so called element, cell or bubble DOFs - the only DOFs used in DG
        n_dof is set as n_cells * order
        remap is setup to map (order) DOFs to each cell
        dofs is ???
        :return:
        """
        self.n_cell = self.region.get_n_cells(self.is_surface)
        n_dof = self.n_cell * (self.approx_order + 1)  # is that right?
        remap = nm.arange(n_dof).reshape((self.n_cell, self.approx_order + 1))  # what is remap used for?
        dofs = nm.arange(n_dof)

        self.econn = dofs[:, None]
        # TODO setup elements connectivity (it should be simple)

        return n_dof, remap, dofs

    def _setup_shape(self):
        """
        What is shape used for and what it really means.
        Does it represent shape of the problem?
        :return:
        """
        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape

    def _setup_geometry(self):
        """
        Setup the field region geometry.
        Somehow pulls the highet dimension geometry from self.region
        """
        # from VolumeField
        cmesh = self.domain.cmesh
        for key, gel in six.iteritems(self.domain.geom_els):
            ct = cmesh.cell_types
            if (ct[self.region.cells] == cmesh.key_to_index[gel.name]).all():
                self.gel = gel
                break

        else:
            raise ValueError('region %s of field %s contains multiple'
                             ' reference geometries!'
                             % (self.region.name, self.name))

    def clear_qp_base(self):
        """
        Remove cached quadrature points and base functions.
        Used in __init__ to set empty qp_coors and bf.
        """
        self.qp_coors = {}
        self.bf = {}

    def setup_extra_data(self, geometry, info, is_trace):
        """
        This is called in create_adof_conns(conn_info, var_indx=None, active_only=True, verbose=True)
        for each variable but has no effect.
        :param geometry:
        :param info:
        :param is_trace:
        :return:
        """
        # TODO placeholder, what is this used for?

        dct = info.dc_type.type # TODO check DOF connectivity type

        self.info = info
        self.is_trace = is_trace

    def get_dofs_in_region(self, region, merge=True):
        """
        Return indices of DOFs that belong to the given region and group.
        NOT really tested, called only with the ragion being the "main" region
        of the problem, i.e. self.region

        :param region:
        :param merge: merge dof tuple into one numpy array
        :return:
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

        # based on get_dofs_in_region from FEField
        dofs = []
        eldofs = nm.empty((0,), dtype=nm.int32)
        if region.has_cells(): # TODO use "and (node_desc.bubble is not None)"
            els = nm.ravel(self.bubble_remap[region.cells])
            eldofs = self.bubble_dofs[els[els >= 0]].ravel()
        dofs.append(eldofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def get_data_shape(self, integral, integration='volume', region_name=None):
        """
        Returns data shape for term, right now it is
        (n_nod, n_qp, self.gel.dim, 1)
        which results in matrix n_nod x n_nod, hovewer in
        FEM it is
        (shape.n_cell, n_qp, dim, self.econn.shape[1])
        how does this translates to domension of the matrix?
        :param integral: integral used
        :param integration:
        :param region_name: not used
        :return:
        """

        if integration in ('volume'):
            # from FEField.get_data_shape()
            _, weights = integral.get_qp(self.gel.name)
            n_qp = weights.shape[0]

            # from FEField data_shape = (shape.n_cell, n_qp, dim, self.econn.shape[1])
            data_shape = (self.n_nod, n_qp, self.gel.dim, 1)

            # TODO last is econn.shape[1]
            # we do not have connectivity in DG but need this enyway, I think

        else:
            # TODO what bout other integrations? do they make sense for DG?
            raise NotImplementedError('unsupported integration! (%s)'
                                      % integration)

        return data_shape

    def get_econn(self, conn_type, region, is_trace=False, integration=None):
        """
        getter for econn
        :param conn_type:
        :param region:
        :param is_trace:
        :param integration: 'volume' is only supported value
        :return:
        """

        ct = conn_type.type if isinstance(conn_type, Struct) else conn_type

        if ct == 'volume':
            if region.name == self.region.name:
                conn = self.econn
            else:
                raise ValueError("Bad region for the field")
        else:
            raise ValueError('unknown connectivity type! (%s)' % ct)

        return conn

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, extend=True, fill_value=None,
                      linearization=None):
        """
        Convert the DOFs corresponding to the field to a dictionary of
        output data usable by Mesh.write().
        # TODO how shloud the dictionary look
        Parameters
        ----------
        dofs : array, shape (n_nod, n_component)
            The array of DOFs reshaped so that each column corresponds
            to one component.
        var_name : str
            The variable name corresponding to `dofs`.
        dof_names : tuple of str
            The names of DOF components.
        key : str, optional
            The key to be used in the output dictionary instead of the
            variable name.
        extend : bool
            Extend the DOF values to cover the whole domain.
        fill_value : float or complex
           The value used to fill the missing DOF values if `extend` is True.
        linearization : Struct or None
            The linearization configuration for higher order approximations.

        Returns
        -------
        out : dict
            The output dictionary.
        """
        res = {}
        for i in range(self.approx_order + 1):
            res["u{}".format(i)] = Struct(mode="cell",
                              data=dofs[self.n_cell * i : self.n_cell*(i+1) ,:, None, None])

                              # nm.hstack((dofs[:self.n_cell],
                              #                 dofs[self.n_cell:])).reshape(self.n_cell, 1,1, 2)
                              #
        return res

    def create_mapping(self, region, integral, integration, return_mapping=True):
        """
        Creates and returns mapping
        :param region:
        :param integral:
        :param integration: 'volume' is so far only accepted option
        :return:
        """

        domain = self.domain
        coors = domain.get_mesh_coors(actual=True)
        dconn = domain.get_conn()
        # from FEField
        if integration == 'volume':
            # TODO qp = self.get_qp('v', integral)
            qp = self.integral.get_qp(self.gel.name)
            iels = region.get_cells()

            geo_ps = self.gel.poly_space
            ps = self.poly_space
            bf = self.get_base('v', 0, integral, iels=iels)

            conn = nm.take(dconn, iels.astype(nm.int32), axis=0)
            mapping = VolumeMapping(coors, conn, poly_space=geo_ps)
            vg = mapping.get_mapping(qp[0], qp[1], poly_space=ps,
                                     ori=self.ori,
                                     transform=self.basis_transform)

            out = vg
        else:
            raise ValueError('unknown integration geometry type: %s'
                             % integration)

        if out is not None:
            # Store the integral used.
            out.integral = integral
            out.qp = qp
            out.ps = ps
            # Update base.
            out.bf[:] = bf

        if return_mapping:
            out = (out, mapping)

        return out

    def get_base(self, key, derivative, integral, iels=None,
                 from_geometry=False, base_only=True):
        """
        Return values of base functions at quadrature points of given integral
        :param key: 'v' - volume, 's' - surface
        :param derivative:
        :param integral:
        :param iels:
        :param from_geometry:
        :param base_only:
        :return:
        """
        # from FEField
        qp = integral.get_qp(self.gel.name)

        if from_geometry:
            ps = self.gel.poly_space

        else:
            ps = self.poly_space

        _key = key if not from_geometry else 'g' + key
        bf_key = (integral.order, _key, derivative)

        if bf_key not in self.bf:
            if (iels is not None) and (self.ori is not None):
                ori = self.ori[iels]

            else:
                ori = self.ori

            self.bf[bf_key] = ps.eval_base(qp[0], diff=derivative, ori=ori,
                                           transform=self.basis_transform)

        if base_only:
            return self.bf[bf_key]
        else:
            return self.bf[bf_key], qp.weights


    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Compute projection of fun into the basis, alternatevely set DOFs directly to provided
        value or values
        :param fun: callable, scallar or array corresponding to dofs
        :param region: region to set DOFs on
        :param dpn: number of dofs per element
        :param warn: not used
        :return: nods, vals
        """

        if region is None:
            region = self.region

        aux = self.get_dofs_in_region(region)
        nods = nm.unique(nm.hstack(aux))

        if nm.isscalar(fun):
            vals = nm.repeat([fun], nods.shape[0] * dpn)

        elif isinstance(fun, nm.ndarray):
            assert_(len(fun) == dpn)
            vals = nm.repeat(fun, nods.shape[0])

        elif callable(fun):

            mesh = region.domain.mesh  # TODO use remap and self.bubble_dofs to get indicies?

            qp, weights = self.integral.get_qp(self.gel.name)
            qp, weights = qp.T, weights[:, None].T # transpose for array expansion
            qp = 2*qp - 1  # in DG we use [-1,1] reference element, change that to [0, 1]?
            weights = 2 * weights  # weights need to be tranformed as well

            def mapping1D(x):
                # TODO this should be somehow part of the mapping
                c = (mesh.coors[1:] + mesh.coors[:-1]) / 2  # center
                s = (mesh.coors[1:] - mesh.coors[:-1]) / 2  # scale
                return c + x * s
            coors = mapping1D(qp)

            # sic = nm.zeros((2, mesh.n_el, 1), dtype=nm.float64)
            # sic[0, :] = nm.sum(weights * fun(coors), axis=1)[:,  None] / 2
            # sic[1, :] = 3 * nm.sum(weights * qp * fun(coors), axis=1)[:,  None] / 2

            # TODO higher order approx seems off
            # base_vals_coors = self.poly_space.eval_base(coors)
            base_vals_qp = self.poly_space.eval_base(qp)       # -2  0   -1  0
            base_vals_qp = nm.swapaxes(nm.swapaxes(base_vals_qp, -2, 0), -1, 0)

            # left hand, so far only orthogonal basis
            lhs_diag = nm.sum(weights * base_vals_qp ** 2, axis=2)
            # right hand
            rhs_vec = nm.sum(weights * base_vals_qp * fun(coors), axis=2)

            vals = rhs_vec / lhs_diag
            # self.plot_1D_dofs(self.domain.mesh.coors, (vals,), fun)
            vals = nm.append(vals[0, :], vals[1, :])
        return nods, vals

    @staticmethod
    def plot_1D_dofs(coors, dofss, fun):
        # TODO move to visualizer
        import matplotlib.pyplot as plt
        X = (coors[1:] + coors[:-1]) / 2
        plt.figure("DOFs for function fun")
        for ii, dofs in enumerate(dofss):
            for i in range(dofs.shape[0]):
                c0 = plt.plot(X, dofs[i, :], label="fun-{}dof-{}".format(ii, i), marker=".", ls="")[0].get_color()
                # # plt.plot(coors, .1*alones(n_nod), marker=".", ls="")
                plt.step(coors[1:], dofs[i, :], color=c0)
                # plt.plot(coors[1:], sic[1, :], label="IC-1", color=c1)

        xs = nm.linspace(0,1, 500)[:, None]
        plt.legend()
        plt.plot(xs, fun(xs), label="fun-ex")
        plt.show()

    def interp_to_qp(self, dofs):
        """
        Interpolate DOFs into quadrature points.

        The quadrature order is given by the field approximation order.

        Parameters
        ----------
        dofs : array
            The array of DOF values of shape `(n_nod, n_component)`.

        Returns
        -------
        data_qp : array
            The values interpolated into the quadrature points.
        integral : Integral
            The corresponding integral defining the quadrature points.
        """
        # TODO this seems useful
        raise NotImplementedError

# _get_facets
# create_basis_context
# create_eval_mesh
# create_mesh
# create_output
# get_true_order
# is_higher_order

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

    from sfepy.discrete.fem import FEDomain
    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')

    dgfield = DGField("dgfu", nm.float64, 'scalar', omega, poly_space_base=None, approx_order=1)

    dg_field_dir = set(dir(dgfield))
    # print("FEM - DG: {}".format(fem_field_dir.difference(dg_field_dir)))
    # print("DG - FEM: {}".format(dg_field_dir.difference(fem_field_dir)))

    u = FieldVariable('u', 'unknown', dgfield, history=1)
    v = FieldVariable('v', 'test', dgfield, primary_var_name='u')

