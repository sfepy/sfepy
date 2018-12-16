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
from dg_basis import LegendrePolySpace, LegendreSimplexPolySpace, LegendreTensorProductPolySpace


def get_unraveler(n_el_nod, n_cell):

    def unravel(u):
        ur = nm.zeros((n_cell, n_el_nod, 1))
        for i in range(n_el_nod):
            ur[:, i] = u[n_cell * i: n_cell * (i + 1), None]
        return ur

    return unravel


def get_raveler(n_el_nod, n_cell):

    def ravel(u):
        ur = nm.zeros((n_cell * n_el_nod, 1))
        for i in range(n_el_nod):
            ur[n_cell * i: n_cell * (i + 1)] = u[:, i]
        return ur

    return ravel


class DGField(Field):
    family_name = 'volume_H1_DGLegendre'
    is_surface = False

    def __init__(self, name, dtype, shape, region, space="H1",
        poly_space_base=None,  approx_order=0, integral=None):
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
        self._setup_geometry()
        self._setup_connectivity()
        self.n_el_facets = self.dim + 1 if self.gel.is_simplex else 2 ** self.dim


        # approximation space
        self.space = space
        self.poly_space_base = poly_space_base
        # TODO put LegendrePolySpace into table in PolySpace any_from_args, or use only Legendre for DG?
        if poly_space_base is not None:
            self.poly_space = poly_space_base("H1_dglegendre", self.gel, approx_order)
        else:
            self.poly_space = LegendreTensorProductPolySpace("1_2_H1_dglegendre", self.gel, approx_order)
        # poly_space = PolySpace.any_from_args("legendre", self.gel, base="legendre", order=approx_order)

        # DOFs
        self.approx_order = approx_order
        self._setup_shape()
        self._setup_all_dofs()

        self.ravel_sol = get_raveler(self.n_el_nod, self.n_cell)
        self.unravel_sol = get_unraveler(self.n_el_nod, self.n_cell)

        # boundary DOFS TODO temporary
        self.boundary_val = 0.0

        # integral
        self.clear_qp_base()
        if integral is None:
            self.integral = Integral("dg_fi", order=2*approx_order)
        else:
            self.integral = integral

        self.ori = None
        self.basis_transform = None

        # mapping
        self.mappings = {}
        self.mapping = self.create_mapping(self.region, self.integral, "volume", return_mapping=True)[1]
        self.mappings0 = {}

    def _check_region(self, region):
        # TODO what are the requirements on region?
        return True

    def _setup_all_dofs(self):
        """
        Sets up all the differet kinds of DOFs, for DG only bubble DOFs
        originaly called _setup_global_base
        """
        self.n_el_nod = self.poly_space.n_nod
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
        n_dof is set as n_cells * n_el_nod
        remap is setup to map (order) DOFs to each cell
        dofs is ???
        :return:
        """
        self.n_cell = self.region.get_n_cells(self.is_surface)
        n_dof = self.n_cell * self.n_el_nod
        dofs = nm.ones((self.n_cell, self.n_el_nod), dtype=nm.int32)
        for i in range(self.n_el_nod):
            dofs[:, i] = nm.arange(self.n_cell*i, self.n_cell*(i+1), dtype=nm.int32)

        remap = nm.arange(self.n_cell)
        self.econn = dofs

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
        self.dim = cmesh.dim
        for key, gel in six.iteritems(self.domain.geom_els):
            ct = cmesh.cell_types
            if (ct[self.region.cells] == cmesh.key_to_index[gel.name]).all():
                self.gel = gel
                break

            else:
                raise ValueError('region %s of field %s contains multiple'
                                 ' reference geometries!'
                                 % (self.region.name, self.name))

    def _setup_connectivity(self):
        """
        Forces self.domain.mesh to build neccesary conectivities
        so the are available in self.get_nbrhd_dofs
        :return:
        """
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim, self.dim)
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim, self.dim - 1)


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
        # placeholder, what is this used for?

        dct = info.dc_type.type

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

        dofs = []
        eldofs = nm.empty((0,), dtype=nm.int32)
        if region.has_cells():
            els = nm.ravel(self.bubble_remap[region.cells])
            eldofs = self.bubble_dofs[els[els >= 0]]
        dofs.append(eldofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def get_nbrhd_dofs(self, region, variable):
        n_el_nod = self.n_el_nod
        n_cell = self.n_cell
        dim = self.dim
        gel = self.gel
        n_el_facets = dim + 1 if gel.is_simplex else 2 ** dim

        nb_dofs = -1 * nm.ones((n_cell,n_el_facets, n_el_nod, 1))
        dofs = self.unravel_sol(variable.data[0])

        neighbours = self.get_cell_nb_per_facet(region)
        nb_normals = self.get_cell_normals_per_facet(region, neighbours)

        ghost_nbrs = nm.where(neighbours < 0)


        # TODO treat boundary conditions
        if dim == 1:  # periodic boundary conditions in 1D
             neighbours[0, 0] = -1
             neighbours[-1, 1] = 0
        nb_dofs[:] = nm.take(dofs, neighbours, axis=0)
        # nb_dofs[ghost_nbrs] = self.boundary_val

        return nb_dofs, nb_normals



    def get_cell_nb_per_facet(self, region):
        """
        Retruns array of cell neighbbours sharing facet
        :param region:
        # TODO why should be region passed?
        :return:
        """
        n_cell = self.n_cell
        dim = self.dim
        gel = self.gel
        n_el_facets = dim + 1 if gel.is_simplex else 2 ** dim

        cmesh = region.domain.mesh.cmesh
        cells = region.cells

        neighbours = nm.zeros((n_cell, n_el_facets), dtype=nm.int32)

        c2fi, c2fo = cmesh.get_incident(dim - 1, cells, dim, ret_offsets=True)
        for ic, o1 in enumerate(c2fo[:-1]):  # loop over cells
            o2 = c2fo[ic + 1]

            c2ci, c2co = cmesh.get_incident(dim, c2fi[o1:o2], dim - 1,
                                            ret_offsets=True)  # get neighbours per facet of the cell
            nbrs = []
            for ifa, of1 in enumerate(c2co[:-1]):  # loop over facets
                of2 = c2co[ifa + 1]
                if of2 == (of1 + 1):  # facet has only one cell
                    # Surface facet.
                    nbrs.append(-1)  # c2ci[of1])  # append the cell, itself
                else:
                    if c2ci[of1] == cells[ic]:  # do not append the cell itself
                        nbrs.append(c2ci[of2 - 1])
                    else:
                        nbrs.append(c2ci[of1])
            neighbours[ic] = nbrs

        return neighbours


    def get_cell_normals_per_facet(self, region, neighbours):
        n_cell = self.n_cell
        dim = self.dim
        gel = self.gel
        n_el_facets = dim + 1 if gel.is_simplex else 2 ** dim

        cmesh = region.domain.mesh.cmesh
        cells = region.cells

        normals = cmesh.get_facet_normals()
        if dim == 1:
            normals[:, 0] = nm.tile([-1, 1], int(normals.shape[0]/2))
        normals_out = nm.zeros((n_cell, n_el_facets, dim))

        c2f = cmesh.get_conn(dim, dim - 1)
        for ic, o1 in enumerate(c2f.offsets[:-1]):
            o2 = c2f.offsets[ic + 1]
            for ifal, ifa in enumerate(c2f.indices[o1:o2]):
                normals_out[ic, ifal] = normals[o1 + ifal]

        return normals_out


    def get_data_shape(self, integral, integration='volume', region_name=None):
        """
        Returns data shape
        (n_nod, n_qp, self.gel.dim, self.n_el_nod)

        :param integral: integral used
        :param integration:
        :param region_name: not used
        :return:
        """

        if integration in ('volume'):
            # from FEField.get_data_shape()
            _, weights = integral.get_qp(self.gel.name)
            n_qp = weights.shape[0]

            data_shape = (self.n_cell, n_qp, self.gel.dim, self.n_el_nod)
            # econn.shape[1] == n_el_nod i.e. number nod in element

        else:
            # what bout other integrations? do they make sense for DG?
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

        Puts DOFs into vairables u0 ... un, where n = approx_order and marks them for writing
        as cell data.

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
            raise ValueError('unsupported integration geometry type: %s'
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

            qp, weights = self.integral.get_qp(self.gel.name)
            weights = weights.reshape(nm.shape(weights) + tuple(nm.ones(nm.ndim(qp), dtype=nm.int32))) # reshape weights to match qp
            coors = self.mapping.get_physical_qps(qp)

            # sic = nm.zeros((2, mesh.n_el, 1), dtype=nm.float64)
            # sic[0, :] = nm.sum(weights * fun(coors), axis=1)[:,  None] / 2
            # sic[1, :] = 3 * nm.sum(weights * qp * fun(coors), axis=1)[:,  None] / 2

            base_vals_qp = self.poly_space.eval_base(qp)
            # base_vals_qp = nm.swapaxes(nm.swapaxes(base_vals_qp, -2, 0), -1, 0)

            # left hand, so far only orthogonal basis
            lhs_diag = nm.sum(weights * base_vals_qp**2, axis=0)
            # for legendre base it is exactly: 1 / (2 * nm.arange(self.n_el_nod) + 1)

            # right hand TODO this is hot fix, check this in 2D and 3D case
            rhs_vec = nm.sum((weights * base_vals_qp)[:, 0, :] * fun(coors), axis=1)

            vals = (rhs_vec / lhs_diag).T

            from my_utils.visualizer import plot_1D_legendre_dofs, reconstruct_legendre_dofs
            import matplotlib.pyplot as plt
            plot_1D_legendre_dofs(self.domain.mesh.coors, (vals,), fun)
            ww, xx = reconstruct_legendre_dofs(self.domain.mesh.coors, 1, vals[..., None, None])
            plt.plot(xx ,ww[:, 0])

        return nods, vals


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

