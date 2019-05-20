"""
Fields for Discontinous Galerkin
"""
import numpy as nm
from numpy.lib.stride_tricks import as_strided
import six

# sfepy imports
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete import Integral, FieldVariable
from six.moves import range
from sfepy.discrete.fem import Mesh, Field
from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.discrete.fem.mappings import VolumeMapping
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)

# local imports
from sfepy.discrete.dg.dg_basis import LegendrePolySpace, LegendreSimplexPolySpace, LegendreTensorProductPolySpace


def get_unraveler(n_el_nod, n_cell):
    def unravel(u):
        ustride1 = u.strides[0]
        ur = as_strided(u, shape=(n_cell, n_el_nod, 1),
                        strides=(n_el_nod * ustride1, ustride1, ustride1), writeable=False)
        # FIXME writeable is not valid option for Python 2
        return ur

    return unravel


def get_raveler(n_el_nod, n_cell):
    def ravel(u):
        # ustride1 = u.strides[0]
        # ur = as_strided(u, shape=(n_el_nod*n_cell, 1),
        #                     strides=(n_cell*ustride1, ustride1))
        ur = nm.ravel(u)[:, None]
        return ur

    return ravel


def get_cell_facet_gel_name(cell_gel_name):
    """
    Returns name pf the facet geometry for given cell geometry
    :param cell_gel_name: name of the cell geometry
    :return: name of the facet geometry
    """
    if cell_gel_name == "1_2":
        return "0_1"
    elif cell_gel_name == "2_3" or cell_gel_name == "2_4":
        return "1_2"
    elif cell_gel_name == "3_4":
        return "2_3"
    elif cell_gel_name == "3_8":
        return "2_4"
    else:
        raise ValueError('unknown geometry type! {}'.format(cell_gel_name))


def get_gel(region):
    """
    :param region: sfepy region
    :return: base geometry element of the region
    """
    cmesh = region.domain.cmesh
    for key, gel in six.iteritems(region.domain.geom_els):
        ct = cmesh.cell_types
        if (ct[region.cells] == cmesh.key_to_index[gel.name]).all():
            return gel
        else:
            raise ValueError('Region {} contains multiple'
                             ' reference geometries!'.format(region))


class DGField(Field):
    family_name = 'volume_DG_legendre_discontinuous'
    is_surface = False

    def __init__(self, name, dtype, shape, region, space="H1",
                 poly_space_base=None, approx_order=0, integral=None):
        """

         Creates DGField, with Legendre poly space and default integral corresponding to
        2*approx_order.

        :param name:
        :param dtype:
        :param shape:  'vector', 'scalar'or something else
        :param region:
        :param space: default "H1"
        :param poly_space_base: optionally force polyspace
        :param approx_order: 0 for FVM
        :param integral: if None integral of order 2*approx_order is created
        """
        shape = parse_shape(shape, region.domain.shape.dim)
        # if not self._check_region(region):
        #     raise ValueError('unsuitable region for field %s! (%s)' %
        #                      (name, region.name))
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)

        self.approx_order = approx_order[0] if isinstance(approx_order, tuple) else approx_order

        # geometry
        self.domain = region.domain
        self.region = region
        self.dim = region.dim
        self._setup_geometry()
        self._setup_connectivity()
        # FIXME treat domains embeded to higher dimensional spaces
        self.n_el_facets = self.dim + 1 if self.gel.is_simplex else 2 ** self.dim

        # approximation space
        self.space = space
        self.poly_space_base = poly_space_base
        if poly_space_base is not None:
            self.poly_space = poly_space_base(self.gel.name + "H1_what?",
                                              self.gel, self.approx_order)
        elif self.gel.name in ["1_2", "2_4", "3_8"]:
            self.poly_space = LegendreTensorProductPolySpace(self.gel.name + "_DG_legendre",
                                                             self.gel, self.approx_order)
        else:
            self.poly_space = LegendreSimplexPolySpace(self.gel.name + "_DG_legendre",
                                                       self.gel, self.approx_order)

        # TODO put LegendrePolySpace into table in PolySpace any_from_args, or use only Legendre for DG?
        # poly_space = PolySpace.any_from_args("legendre", self.gel, base="legendre", order=approx_order)

        # DOFs
        self._setup_shape()
        self._setup_all_dofs()

        self.ravel_sol = get_raveler(self.n_el_nod, self.n_cell)
        self.unravel_sol = get_unraveler(self.n_el_nod, self.n_cell)

        # integral
        self.clear_qp_base()
        self.clear_facet_qp_base()
        if integral is None:
            self.integral = Integral("dg_fi", order=2 * self.approx_order)
        else:
            self.integral = integral

        self.ori = None
        self.basis_transform = None

        # mapping
        self.mappings = {}
        self.mapping = self.create_mapping(self.region, self.integral, "volume", return_mapping=True)[1]
        self.mappings0 = {}

        # neighbour facet mapping and data
        self.clear_facet_neighbour_idx_cache()
        self.clear_normals_cache()
        self.clear_facet_vols_cache()
        self.boundary_facet_local_idx = {}


    def _setup_all_dofs(self):
        """
        Sets up all the differet kinds of DOFs, for DG only bubble DOFs
        originaly called _setup_global_base
        """
        self.n_el_nod = self.poly_space.n_nod
        self.n_vertex_dof = 0  # in DG we will propably never need vertex DOFs
        self.n_edge_dof = 0  # use facets DOFS for AFS methods
        self.n_face_dof = 0  # use facet DOF for AFS methods

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
        # for i in range(self.n_el_nod):
        #     dofs[:, i] = nm.arange(self.n_cell*i, self.n_cell*(i+1), dtype=nm.int32)

        dofs = nm.arange(n_dof, dtype=nm.int32).reshape(self.n_cell, self.n_el_nod)
        remap = nm.arange(self.n_cell)
        self.econn = dofs
        self.dofs2cells = nm.repeat(nm.arange(self.n_cell), self.n_el_nod)

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
        self.gel = get_gel(self.region)

    def _setup_connectivity(self):
        """
        Forces self.domain.mesh to build neccesary conectivities
        so the are available in self.get_nbrhd_dofs
        :return:
        """
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim, self.dim)
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim - 1, self.dim)
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim, self.dim - 1)

    def clear_qp_base(self):
        """
        Remove cached quadrature points and base functions.
        Used in __init__ to set empty qp_coors and bf.
        """
        self.qp_coors = {}
        self.bf = {}

    def get_qp(self, key, integral):
        """
        Get quadrature points and weights corresponding to the given key
        and integral. The key is 'v' or 's#', where # is the number of
        face vertices.
        """
        qpkey = (integral.order, key)

        if qpkey not in self.qp_coors:
            if (key[0] == 's') and not self.is_surface:
                dim = self.gel.dim - 1
                n_fp = self.gel.surface_facet.n_vertex
                geometry = '%d_%d' % (dim, n_fp)

            else:
                geometry = self.gel.name

            vals, weights = integral.get_qp(geometry)
            self.qp_coors[qpkey] = Struct(vals=vals, weights=weights)

        return self.qp_coors[qpkey]

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

    def get_coor(self, nods=None):
        """
        Returns coors for matching nodes, uses trick to deceive
        EPBC implementation in sfepy
        :param nods: if None use all nodes
        :return:
        """

        if nods is None:
            nods = self.bubble_dofs

        cells = self.dofs2cells[nods]
        coors = self.domain.mesh.cmesh.get_centroids(self.dim)[cells]
        eps = min(self.domain.cmesh.get_volumes(self.dim)) / (self.n_el_nod + 2)
        if self.dim == 1:
            extended_coors = nm.zeros(nm.shape(coors)[:-1] + (2,))
            extended_coors[:, 0] = coors[:, 0]
            coors = extended_coors
        # shift centroid coors to lie within cells but be different for each dof
        # TODO for simplex meshes these coors fail to match
        coors += eps * nm.repeat(nm.arange(self.n_el_nod), len(nm.unique(cells)))[:, None]
        return coors

    def clear_facet_qp_base(self):
        self.facet_bf = {}
        self.facet_qp = None
        self.facet_whs = None

    def _transform_qps_to_facets(self, qps, geo_name):
        """
        Transforms points given in qps to all facets of the reference element
        with geometry geo_name.
        :param qps:
        :param geo_name:
        :return: tqps is of shape shape(qps) + (n_el_facets, geo dim)
        """
        if geo_name == "1_2":
            tqps = nm.zeros(nm.shape(qps) + (2, 1,))
            tqps[..., 0, 0] = 0.
            tqps[..., 1, 0] = 1.
        elif geo_name == "2_3":
            tqps = nm.zeros(nm.shape(qps) + (3, 2,))
            # 0.
            tqps[..., 0, 0] = qps  # x = 0 + t
            tqps[..., 0, 1] = 0.  # y = 0
            # 1.
            tqps[..., 1, 0] = 1 - qps  # x = 1 - t
            tqps[..., 1, 1] = qps  # y = t
            # 2.
            tqps[..., 2, 0] = 0  # x = 0
            tqps[..., 2, 1] = 1 - qps  # y = 1 - t
        elif geo_name == "2_4":
            tqps = nm.zeros(nm.shape(qps) + (4, 2,))
            # 0.
            tqps[..., 0, 0] = qps  # x = t
            tqps[..., 0, 1] = 0.  # y = 0
            # 1.
            tqps[..., 1, 0] = 1  # x = 1
            tqps[..., 1, 1] = qps  # y = t
            # 2.
            tqps[..., 2, 0] = 1 - qps  # x = 1 -t
            tqps[..., 2, 1] = 1  # y = 1
            # 3.
            tqps[..., 3, 0] = 0  # x = 0
            tqps[..., 3, 1] = 1 - qps  # y = 1 - t
        elif geo_name == "3_4":
            tqps = nm.zeros(nm.shape(qps) + (4, 3,))
            raise NotImplementedError("Geometry {} not supported, yet".format(geo_name))
        elif geo_name == "3_8":
            tqps = nm.zeros(nm.shape(qps) + (8, 3,))
            raise NotImplementedError("Geometry {} not supported, yet".format(geo_name))
        else:
            raise NotImplementedError("Geometry {} not supported, yet".format(geo_name))
        return tqps

    def get_facet_qp(self):
        """
        Returns dim - 1 quadrature points on all facets of the reference element in array of shape
        (n_qp, n_el_facets, dim)
        :return: qp, weights - need to be transformed to actual facets!
        """

        if self.dim == 1:
            facet_qps = self._transform_qps_to_facets(nm.zeros((1, 1)), "1_2")
            weights = nm.ones((1, 1, 1))
        else:
            qps, weights = self.integral.get_qp(get_cell_facet_gel_name(self.gel.name))
            weights = weights[None, :, None]
            facet_qps = self._transform_qps_to_facets(qps, self.gel.name)

        # from postprocess.plot_facets import plot_geometry
        # from postprocess.plot_quadrature import plot_weighted_points
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # plot_geometry(ax, self.gel)
        # facet_qps_flat = nm.concatenate([facet_qps[..., i,:] for i in range(self.n_el_facets)])
        # facet_weights_flat = nm.concatenate([weights] * self.n_el_facets)[:, 0]
        # ax.scatter(facet_qps_flat[..., 0, 0], facet_qps_flat[..., 0, 1])
        # plt.show()

        return facet_qps, weights

    def get_facet_base(self, derivative=False, base_only=False):
        """
        Returns values of base in facets quadrature points, data shape is a bit crazy right now
        currently (number of qps, 1, n_el_facets, 1, n_el_nod). This is because base eval preserves qp shape and ads
        dimension of the value - in case of derivative this will be (dim,) * derivative order and all basis values
        i.e. n_el_nod values
        :param derivative:
        :param base_only:
        :return:
        """
        if derivative:
            diff = int(derivative)
        else:
            diff = 0

        if diff in self.facet_bf:
            facet_bf = self.facet_bf[diff]
            whs = self.facet_whs
        else:
            qps, whs = self.get_facet_qp()
            ps = self.poly_space
            self.facet_qp = qps
            self.facet_whs = whs
            if derivative:
                facet_bf = nm.zeros((1,) + nm.shape(qps)[:-1] + (self.dim,) * diff + (self.n_el_nod,))
            else:
                facet_bf = nm.zeros(nm.shape(qps)[:-1] + (1, self.n_el_nod,))

            for i in range(self.n_el_facets):
                facet_bf[..., i, :, :] = ps.eval_base(qps[..., i, :], diff=diff,
                                                      transform=self.basis_transform)
            self.facet_bf[diff] = facet_bf

        if base_only:
            return facet_bf
        else:
            return facet_bf, whs

    def clear_facet_neighbour_idx_cache(self, region=None):
        """
        If region is None clear all!
        :param region:
        :return:
        """
        if region is None:
            self.facet_neighbour_index = {}
        else:
            self.facet_neighbour_index.remove(region.name)

    def get_facet_neighbor_idx(self, region, eq_map):
        """
        Returns index of cell neighbours sharing facet, along with local index
        of the facet within neighbour also treats periodic boundary conditions i.e.,
        plugs correct neighbours for cell on periodic boundary. Where there are no neighbours
        specified puts -1.

        Cashes neighbour index in self.facet_neighbours

        :param region:
        :param eq_map: eq_map from state variable containing information on EPBC
        :return: shape is (n_cell, n_el_facet, 2), first value in last axis is index of the neighbouring cell
        the second is index of the facet this nb. cell in said nb. cell
        """
        if region.name in self.facet_neighbour_index:
            return self.facet_neighbour_index[region.name]

        dim, n_cell, n_el_facets = self.get_region_info(region)

        cmesh = region.domain.mesh.cmesh
        cells = region.cells

        facet_neighbours = nm.zeros((n_cell, n_el_facets, 2), dtype=nm.int32)

        c2fi, c2fo = cmesh.get_incident(dim - 1, cells, dim, ret_offsets=True)

        for ic, o1 in enumerate(c2fo[:-1]):  # loop over cells
            o2 = c2fo[ic + 1]

            c2ci, c2co = cmesh.get_incident(dim, c2fi[o1:o2], dim - 1,
                                            ret_offsets=True)  # get neighbours per facet of the cell
            ii = cmesh.get_local_ids(c2fi[o1:o2], dim - 1, c2ci, c2co, dim)
            fis = nm.c_[c2ci, ii]

            nbrs = []
            for ifa, of1 in enumerate(c2co[:-1]):  # loop over facets
                of2 = c2co[ifa + 1]
                if of2 == (of1 + 1):  # facet has only one cell
                    # Surface facet.
                    nbrs.append([-1, -1])  # c2ci[of1])  # append no neighbours
                else:
                    if c2ci[of1] == cells[ic]:  # do not append the cell itself
                        nbrs.append(fis[of2 - 1])
                    else:
                        nbrs.append(fis[of1])
            facet_neighbours[ic, :, :] = nbrs

        # treat classical FEM EPBCs - we need to correct neighbours
        if eq_map.n_epbc > 0:
            # set neighbours of periodic cells to one another
            mcells = nm.unique(self.dofs2cells[eq_map.master])
            scells = nm.unique(self.dofs2cells[eq_map.slave])
            mcells_facets = nm.array(nm.where(facet_neighbours[mcells] == -1))[1, 0]  # facets of mcells
            scells_facets = nm.array(nm.where(facet_neighbours[scells] == -1))[1, 0]  # facets of scells
            [1, 0]  # above, first we need second axis to get axis on which facet indices are stored,
            # second we drop axis with neighbour local facet index, for multiple s/mcells this will have to be something like
            # 1 + 2*nm.arange(len(mcells)) - to skip double entries for -1 tags in neighbours and  neighbour local facet idx

            facet_neighbours[mcells, mcells_facets, 0] = scells  # set neighbours of mcells to scells
            facet_neighbours[
                mcells, mcells_facets, 1] = scells_facets  # set neighbour facets to facets of scell missing neighbour
            # TODO how to distinguish EBC and EPBC? - we do not need to, EBC over write EPBC, we only need to fix shapes

            facet_neighbours[scells, scells_facets, 0] = mcells  # set neighbours of scells to mcells
            facet_neighbours[
                scells, scells_facets, 1] = mcells_facets  # set neighbour facets to facets of mcell missing neighbour

        # treat DG EPBC - these are definetly prefered
        for master_bc2bfi, slave_bc2bfi in eq_map.dg_epbc:
            # set neighbours of periodic cells to one another
            facet_neighbours[master_bc2bfi[:, 0], master_bc2bfi[:, 1], 0] = slave_bc2bfi[:, 0]
            facet_neighbours[slave_bc2bfi[:, 0], slave_bc2bfi[:, 1], 0] = master_bc2bfi[:, 0]

            # set neigbours facets
            facet_neighbours[slave_bc2bfi[:, 0], slave_bc2bfi[:, 1], 1] = master_bc2bfi[:, 1]
            facet_neighbours[master_bc2bfi[:, 0], master_bc2bfi[:, 1], 1] = slave_bc2bfi[:, 1]

            # cache results
        self.facet_neighbour_index[region.name] = facet_neighbours

        return facet_neighbours

    def get_region_info(self, region):
        """
        Extracts information about region needed in various methods of DGField
        :param region:
        :return: dim, n_cell, n_el_facets
        """
        if not region.has_cells():
            raise ValueError("Region {} has no cells".format(region.name))
        n_cell = region.get_n_cells()
        dim = region.dim
        gel = get_gel(region)
        n_el_facets = dim + 1 if gel.is_simplex else 2 ** dim
        return dim, n_cell, n_el_facets

    def get_both_facet_qp_vals(self, state, region, derivative=None, reduce_nod=True):
        """
        Computes values of the variable represented by dofs in
        quadrature points located at facets, returns both values -
        inner and outer, along with weights.
        :param state: state variable containing BC info
        :param region:
        :param derivative: compute derivative if truthy, compute n-th derivative if number
        :return: inner_facet_values (n_cell, n_el_facets, n_qp), outer facet values (n_cell, n_el_facets, n_qp), weights
                 if derivative is True: inner_facet_values (n_cell, n_el_facets, dim, n_qp),
                                        outer facet values (n_cell, n_el_facets, dim, n_qp)

        """

        if derivative:
            diff = int(derivative)
        else:
            diff = 0
        unreduce_nod = int(not reduce_nod)

        inner_base_vals, outer_base_vals, whs = self.get_both_facet_base_vals(state, region, derivative=derivative)
        dofs = self.unravel_sol(state.data[0])

        n_qp = whs.shape[-1]
        outputs_shape = (self.n_cell, self.n_el_facets) + (self.n_el_nod,) * unreduce_nod + (self.dim,) * diff + (n_qp,)

        inner_facet_vals = nm.zeros(outputs_shape)
        if unreduce_nod:
            inner_facet_vals[:] = nm.einsum('id...,idf...->ifd...', dofs, inner_base_vals)
        else:
            inner_facet_vals[:] = nm.einsum('id...,id...->i...', dofs, inner_base_vals)

        per_facet_neighbours = self.get_facet_neighbor_idx(region, state.eq_map)

        outer_facet_vals = nm.zeros(outputs_shape)
        for facet_n in range(self.n_el_facets):
            if unreduce_nod:
                outer_facet_vals[:, facet_n, :] = nm.einsum('id...,id...->id...',
                                                            dofs[per_facet_neighbours[:, facet_n, 0]],
                                                            outer_base_vals[:, :, facet_n])
            else:
                outer_facet_vals[:, facet_n, :] = nm.einsum('id...,id...->i...',
                                                            dofs[per_facet_neighbours[:, facet_n, 0]],
                                                            outer_base_vals[:, :, facet_n])

        boundary_cells = nm.array(nm.where(per_facet_neighbours[:, :, 0] < 0)).T
        # TODO check and print boundary cells without defined BCs
        for ebc, ebc_vals in zip(state.eq_map.dg_ebc, state.eq_map.dg_ebc_val):
            outer_facet_vals[ebc[:, 0], ebc[:, 1], :] = nm.einsum("id,id...->i...",
                                                                  ebc_vals,
                                                                  inner_base_vals[0, :, ebc[:, 1]])

        # flip outer_facet_vals moved to get_both_facet_base_vals
        return inner_facet_vals, outer_facet_vals, whs

    def get_both_facet_base_vals(self, state, region, derivative=None):
        """
        Returns values of the basis function in quadrature points on facets broadcasted to all
        cells inner to the element as well as outer ones along with weights for the qps broadcasted
        and transformed to cells
        to elements
        :param state: used to get EPBC info
        :param region: for connectivity
        :param derivative: if u need derivative
        :return: inner and outer base vals, shape: (n_cell, n_el_nod, n_el_facet, n_qp) or
                        (n_cell, n_el_nod, n_el_facet, dim, n_qp) whe derivative is True or 1
                        whs, shape: (n_cell, n_el_facet, n_qp)
        """
        if derivative:
            diff = int(derivative)
        else:
            diff = 0

        facet_bf, whs = self.get_facet_base(derivative=derivative)
        n_qp = nm.shape(whs)[1]
        facet_vols = self.get_facet_vols(region)
        whs = facet_vols * whs[None, :, :, 0]

        base_shape = (self.n_cell, self.n_el_nod, self.n_el_facets) + (self.dim,) * diff + (n_qp,)
        inner_facet_base_vals = nm.zeros(base_shape)
        outer_facet_base_vals = nm.zeros(base_shape)

        if derivative:
            inner_facet_base_vals[:] = facet_bf[0, :, 0, :, :, :].swapaxes(-2, -3).T
        else:
            inner_facet_base_vals[:] = facet_bf[:, 0, :, 0, :].T

        per_facet_neighbours = self.get_facet_neighbor_idx(region, state.eq_map)
        # numpy prepends shape resulting from multiple indexing before remaining shape
        if derivative:
            outer_facet_base_vals[:] = inner_facet_base_vals[0, :, per_facet_neighbours[:, :, 1]].swapaxes(-3, -4)
        else:
            outer_facet_base_vals[:] = inner_facet_base_vals[0, :, per_facet_neighbours[:, :, 1]].swapaxes(-2, -3)

        return inner_facet_base_vals, outer_facet_base_vals[..., ::-1], whs

    def clear_normals_cache(self, region=None):
        if region is None:
            self.normals_cache = {}
        else:
            if isinstance(region, str):
                self.normals_cache.remove(region)
            else:
                self.normals_cache.remove(region.name)

    def get_cell_normals_per_facet(self, region):
        """
        Caches resuts, use clear_normals_cache to clear the cache.

        :param region:
        :return: normals of facets in array of shape (n_cell, n_el_facets, dim)
        """

        if region.name in self.normals_cache:
            return self.normals_cache[region.name]

        dim, n_cell, n_el_facets = self.get_region_info(region)

        cmesh = region.domain.mesh.cmesh
        cells = region.cells

        normals = cmesh.get_facet_normals()
        if dim == 1:
            normals[:, 0] = nm.tile([-1, 1], int(normals.shape[0] / 2))
        normals_out = nm.zeros((n_cell, n_el_facets, dim))

        c2f = cmesh.get_conn(dim, dim - 1)
        for ic, o1 in enumerate(c2f.offsets[:-1]):
            o2 = c2f.offsets[ic + 1]
            for ifal, ifa in enumerate(c2f.indices[o1:o2]):
                normals_out[ic, ifal] = normals[o1 + ifal]

        self.normals_cache[region.name] = normals_out

        return normals_out

    def clear_facet_vols_cache(self, region=None):
        if region is None:
            self.facet_vols_cache = {}
        else:
            if isinstance(region, str):
                self.facet_vols_cache.remove(region)
            else:
                self.facet_vols_cache.remove(region.name)

    def get_facet_vols(self, region):
        """

        Caches results, use clear_facet_vols_cache to clear the cach

        :param region:
        :return: volumes of the facets by cells shape is (n_cell, n_el_facets, 1)
        """

        if region.name in self.facet_vols_cache:
            return self.facet_vols_cache[region.name]

        dim, n_cell, n_el_facets = self.get_region_info(region)

        cmesh = region.domain.mesh.cmesh
        cells = region.cells

        if dim == 1:
            vols = nm.ones((cmesh.num[0], 1))
            vols[:, 0] = nm.tile([1, 1], int(vols.shape[0] / 2))
        else:
            vols = cmesh.get_volumes(dim - 1)[:, None]

        vols_out = nm.zeros((n_cell, n_el_facets, 1))

        c2f = cmesh.get_conn(dim, dim - 1)

        for ic, o1 in enumerate(c2f.offsets[:-1]):
            o2 = c2f.offsets[ic + 1]
            for ifal, ifa in enumerate(c2f.indices[o1:o2]):
                vols_out[ic, ifal] = vols[ifa]

        self.facet_vols_cache[region.name] = vols_out

        return vols_out

    def get_data_shape(self, integral, integration='volume', region_name=None):
        """
        Returns data shape
        (n_nod, n_qp, self.gel.dim, self.n_el_nod)

        :param integral: integral used
        :param integration:
        :param region_name: not used
        :return:
        """

        if integration in ('volume',):
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

        Used in BC treatement

        :param region:
        :param merge: merge dof tuple into one numpy array
        :return:
        """

        dofs = []
        if region.has_cells():  # main region or its part
            els = nm.ravel(self.bubble_remap[region.cells])
            eldofs = self.bubble_dofs[els[els >= 0]]
            dofs.append(eldofs)
        else:
            # return indicies of cells adjacent to boundary facets
            dim = self.dim
            cmesh = region.domain.mesh.cmesh
            bc_cells = cmesh.get_incident(dim, region.facets, dim - 1)
            bc_dofs = self.bubble_dofs[bc_cells]
            dofs.append(bc_dofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def get_facet_boundary_index(self, region):
        """

        Caches results in self.boundary_facet_local_idx

        :param region: surface region defining BCs
        :return: index of cells on boundary along with corresponding facets
        """

        if region.name in self.boundary_facet_local_idx:
            return self.boundary_facet_local_idx[region.name]

        bc2bfi = region.get_facet_indices()
        self.boundary_facet_local_idx[region.name] = bc2bfi

        return bc2bfi



    def create_mapping(self, region, integral, integration, return_mapping=True):
        """
        Creates and returns mapping
        :param return_mapping: default True
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
            qp = self.get_qp('v', integral)
            # qp = self.integral.get_qp(self.gel.name)
            iels = region.get_cells()

            geo_ps = self.gel.poly_space
            ps = self.poly_space
            bf = self.get_base('v', 0, integral, iels=iels)

            conn = nm.take(dconn, iels.astype(nm.int32), axis=0)
            mapping = VolumeMapping(coors, conn, poly_space=geo_ps)
            vg = mapping.get_mapping(qp.vals, qp.weights, poly_space=ps,
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

    def get_nbrhd_dofs(self, region, variable):
        """
        Puts -1 where cells has no neighbour
        :param region:
        :param variable:
        :return: (n_cell, n_el_facets, n_el_nod, 1)
        """

        n_el_nod = self.n_el_nod
        n_cell = self.n_cell
        dim = self.dim
        gel = self.gel
        n_el_facets = dim + 1 if gel.is_simplex else 2 ** dim

        nb_dofs = -1 * nm.ones((n_cell, n_el_facets, n_el_nod, 1))

        dofs = self.unravel_sol(variable.data[0])

        neighbours = self.get_facet_neighbor_idx(region)[..., 0]
        nb_normals = self.get_cell_normals_per_facet(region)

        ghost_nbrs = nm.where(neighbours < 0)

        return nb_dofs, nb_normals

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Compute projection of fun into the basis, alternatively set DOFs directly to provided
        value or values either in main volume region or in boundary region
        :param fun: callable, scallar or array corresponding to dofs
        :param region: region to set DOFs on
        :param dpn: number of dofs per element
        :param warn: not used
        :return: nods, vals
        """

        if region is None:
            region = self.region
            return self.set_cell_dofs(fun, region, dpn, warn)
        elif region.has_cells():
            return self.set_cell_dofs(fun, region, dpn, warn)
        elif region.kind_tdim == self.dim - 1:
            nods, vals = self.set_facet_dofs(fun, region, dpn, warn)
            return nods, vals

    def set_cell_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Compute projection of fun onto the basis, in main region, alternatively
        set DOFs directly to provided value or values
        :param fun: callable, scallar or array corresponding to dofs
        :param region: region to set DOFs on
        :param dpn: number of dofs per element
        :param warn: not used
        :return: nods, vals
        """

        aux = self.get_dofs_in_region(region)
        nods = nm.unique(nm.hstack(aux))

        if nm.isscalar(fun):
            vals = nm.zeros(aux.shape)
            vals[:, 0] = fun
            vals = nm.hstack(vals)

        elif isinstance(fun, nm.ndarray):
            # useful for testing, allows to pass complete array of dofs as IC
            if nm.shape(fun) == nm.shape(nods):
                vals = fun

        elif callable(fun):

            qp, weights = self.integral.get_qp(self.gel.name)
            coors = self.mapping.get_physical_qps(qp)

            base_vals_qp = self.poly_space.eval_base(qp)[:, 0, :]
            # this drops redundant axis that is returned by eval_base due to consistency with derivatives

            # left hand, so far only orthogonal basis
            # for legendre base this can be calculated exactly
            # in 1D it is: 1 / (2 * nm.arange(self.n_el_nod) + 1)
            lhs_diag = nm.einsum("q,q...->...", weights, base_vals_qp ** 2)

            rhs_vec = nm.einsum("q,q...,iq...->i...", weights, base_vals_qp, fun(coors))

            vals = (rhs_vec / lhs_diag)

            # plot for 1D
            # from my_utils.visualizer import plot_1D_legendre_dofs, reconstruct_legendre_dofs
            # import matplotlib.pyplot as plt
            # plot_1D_legendre_dofs(self.domain.mesh.coors, (vals,), fun)
            # ww, xx = reconstruct_legendre_dofs(self.domain.mesh.coors, 1, vals.T[..., None, None])
            # plt.plot(xx, ww[:, 0], label="reconstructed dofs")
            # plt.show()

        return nods, vals

    def set_facet_dofs(self, fun, region, dpn, warn):
        """
        Compute projection of fun onto the basis, in main region, alternatively
        set DOFs directly to provided value or values
        :param fun: callable, scallar or array corresponding to dofs
        :param region: region to set DOFs on
        :param dpn: number of dofs per element
        :param warn: not used
        :return: nods, vals
        """

        aux = self.get_dofs_in_region(region)
        nods = nm.unique(nm.hstack(aux))

        if nm.isscalar(fun):
            vals = nm.zeros(aux.shape)
            vals[:, 0] = fun
            vals = nm.hstack(vals)

        elif isinstance(fun, nm.ndarray):
            assert_(len(fun) == dpn)
            vals = nm.zeros(aux.shape)
            vals[:, 0] = nm.repeat(fun, vals.shape[0])

        elif callable(fun):
            vals = nm.zeros(aux.shape)
            # set zero DOF to value fun, set other DOFs to zero
            # FIXME only temporary to test BCs
            vals[:, 0] = fun(1)

            # get facets QPs
            # get facets weights
            # get facet basis vals
            # get coors
            # solve for boundary cell DOFs

        return nods, vals

    def get_nodal_values(self, dofs, region, ref_nodes=None):
        """
        Computes nodal representation of the DOFs
        :param dofs:
        :param region: will we use this?
        :param ref_nodes: defaults to proper set of nodes to get best interpolant properties
        :return:
        """
        if ref_nodes is None:
            # TODO we need exactly as many ref_nodes as we have basis functions, maybe poly_space should provide those?
            ref_nodes = self.get_qp('v', Integral("I", order=self.approx_order + 1)).vals
            # ref_nodes = nm.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=nm.float64)
        base_vals_node = self.poly_space.eval_base(ref_nodes)[:, 0, :]
        dofs = self.unravel_sol(dofs[:, 0])

        nodal_vals = nm.sum(dofs * base_vals_node.T, axis=1)
        nodes = self.mapping.get_physical_qps(ref_nodes)

        # import matplotlib.pyplot as plt
        # plt.plot(nodes[:, 0], nodal_vals)
        # plt.show()

        return nodes, nodal_vals

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, extend=True, fill_value=None,
                      linearization=None):
        """
        Convert the DOFs corresponding to the field to a dictionary of
        output data usable by Mesh.write().

        Puts DOFs into vairables u0 ... un, where n = approx_order and marks them for writing
        as cell data.

        Also get node values and adds them to dictionary as cell_nodes

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
        udofs = self.unravel_sol(dofs)

        for i in range(self.n_el_nod):
            res["u_modal{}".format(i)] = Struct(mode="cell",
                                                data=udofs[:, i, None, None])

        unravel = get_unraveler(self.n_el_nod, self.n_cell)
        res["u_modal_cell_nodes"] = Struct(mode="cell_nodes",
                                           data=unravel(dofs)[..., 0],
                                           interpolation_scheme=self.poly_space.get_interpol_scheme())
        # TODO somehow choose nodal vs modal output
        # cell_nodes, nodal_dofs = self.get_nodal_values(dofs, None, None)
        # res["u_nodal"] = Struct(mode="cell_nodes", data=nodal_dofs)
        return res
