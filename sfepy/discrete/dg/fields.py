# -*- coding: utf-8 -*-
"""
Fields for Discontinous Galerkin method
"""
import numpy as nm
from numpy.lib.stride_tricks import as_strided
import six
from six.moves import range

# sfepy imports
from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)
from sfepy.discrete.fem.fields_base import FEField
from sfepy.discrete import Integral, FieldVariable
from sfepy.discrete.fem.mappings import VolumeMapping
from sfepy.discrete.common.fields import parse_shape

# local imports
from sfepy.discrete.dg.dg_poly_spaces import LegendreSimplexPolySpace
from sfepy.discrete.dg.dg_poly_spaces import LegendreTensorProductPolySpace


def get_unraveler(n_el_nod, n_cell):
    """
    Returns function for unraveling i.e. unpacking dof data from
    serialized array from shape (n_el_nod*n_cell, 1) to (n_cell, n_el_nod, 1).

    The unraveler returns non-writeable view into the input array.

    :param n_el_nod:
    :param n_el_nod, n_cell: expected dimensions of dofs array
    :return:
    """
    def unravel(u):
        """
        Returns non-writeable view into the input array reshaped (n*m, 1)
        to (m, n, 1) .
        :param u:
        :return:
        """
        ustride1 = u.strides[0]
        ur = as_strided(u,
                        shape=(n_cell, n_el_nod, 1),
                        strides=(n_el_nod * ustride1, ustride1, ustride1),
                        writeable=False)
        return ur

    return unravel


def get_raveler(n_el_nod, n_cell):
    """
    Returns function for raveling i.e. packing dof data from
    two dimensional array of shape (n_cell, n_el_nod, 1) to (n_el_nod*n_cell, 1)

    The raveler returns view into the input array.

    :param n_el_nod:
    :param n_el_nod, n_cell: expected dimensions of dofs array
    """
    def ravel(u):
        """
        Returns view into the input array reshaped from (m, n, 1) to (n*m, 1)
        to (m, n, 1) .
        :param u:
        :return:
        """

        # ustride1 = u.strides[0]
        # ur = as_strided(u, shape=(n_el_nod*n_cell, 1),
        #                     strides=(n_cell*ustride1, ustride1))
        ur = nm.ravel(u)[:, None]
        # possibly use according to
        #  https://docs.scipy.org/doc/numpy/reference/generated/numpy.ravel.html
        # ur = u.reshape(-1)
        return ur

    return ravel


# mapping between geometry element types
# and their facets types
# TODO move to sfepy/discrete/fem/geometry_element.py?
cell_facet_gel_name = {
    "1_2": "0_1",
    "2_3": "1_2",
    "2_4": "1_2",
    "3_4": "2_3",
    "3_8": "2_4"
}


def get_gel(region):
    """
    region : sfepy.discrete.common.region.Region
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


class DGField(FEField):
    """
    Class for usage with DG terms, provides functionality for Discontinous
    Galerkin method like neighbour look up, projection to discontinuous basis
    and correct DOF treatment.
    """
    family_name = 'volume_DG_legendre_discontinuous'
    is_surface = False

    def __init__(self, name, dtype, shape, region, space="H1",
                 poly_space_base="legendre", approx_order=1, integral=None):
        """
        Creates DGField, with Legendre polyspace and default integral
        corresponding to 2 * approx_order.

        :param name:
        :param dtype:
        :param shape:  'vector', 'scalar' or something else
        :param region : sfepy.discrete.common.region.Region
        :param space: default "H1"
        :param poly_space_base: optionally force polyspace
        :param approx_order: 0 for FVM, default 1
        :param integral: if None integral of order 2*approx_order is created
        """
        shape = parse_shape(shape, region.domain.shape.dim)
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)

        if isinstance(approx_order, tuple):
            self.approx_order = approx_order[0]
        else:
            self.approx_order = approx_order

        # geometry
        self.domain = region.domain
        self.region = region
        self.dim = region.tdim
        self._setup_geometry()
        self._setup_connectivity()
        # TODO treat domains embedded into higher dimensional spaces?
        self.n_el_facets = self.dim + 1 if self.gel.is_simplex else 2**self.dim

        # approximation space
        self.space = space
        self.poly_space_base = poly_space_base
        self.force_bubble = False
        if self.gel.name in ["1_2", "2_4", "3_8"]:
            self.extended = True
            self.poly_space = LegendreTensorProductPolySpace(
                self.gel.name + "_DG_legendre",
                self.gel, self.approx_order)
        else:
            self.poly_space = LegendreSimplexPolySpace(
                self.gel.name + "_DG_legendre",
                self.gel, self.approx_order)

        # TODO make LegendrePolySpace work through PolySpace any_from_args, or
        #  use only Legendre for DG?
        # poly_space = PolySpace.any_from_args("legendre", self.gel,
        #                                   base="legendre", order=approx_order)
        self._create_interpolant()

        # DOFs
        self._setup_shape()
        self._setup_all_dofs()

        self.ravel_sol = get_raveler(self.n_el_nod, self.n_cell)
        self.unravel_sol = get_unraveler(self.n_el_nod, self.n_cell)

        # integral
        self.clear_qp_base()
        self.clear_facet_qp_base()
        if integral is None:
            self.integral = Integral("dg_fi", order=3 * self.approx_order)
        else:
            self.integral = integral

        self.ori = None
        self.basis_transform = None

        # mapping
        self.mappings = {}
        self.mapping = self.create_mapping(self.region, self.integral, "volume",
                                           return_mapping=True)[1]
        self.mappings0 = {}

        # neighbour facet mapping and data caches
        # TODO use lru cache or different method?
        self.clear_facet_neighbour_idx_cache()
        self.clear_normals_cache()
        self.clear_facet_vols_cache()
        self.boundary_facet_local_idx = {}

    def _create_interpolant(self):
        name = '%s_%s_%s_%d%s' % (self.gel.name, self.space,
                                  self.poly_space_base, self.approx_order,
                                  'B' * self.force_bubble)
        ps = PolySpace.any_from_args(name, self.gel, self.approx_order,
                                     base=self.poly_space_base,
                                     force_bubble=self.force_bubble)
        self.poly_space = ps

    def _setup_all_dofs(self):
        """
        Sets up all the differet kinds of DOFs, for DG only bubble DOFs
        """
        self.n_el_nod = self.poly_space.n_nod
        self.n_vertex_dof = 0  # in DG we will propably never need vertex DOFs
        self.n_edge_dof = 0  # use facets DOFS for AFS methods
        self.n_face_dof = 0  # use facet DOF for AFS methods

        (self.n_bubble_dof,
         self.bubble_remap,
         self.bubble_dofs) = self._setup_bubble_dofs()

        self.n_nod = self.n_vertex_dof + self.n_edge_dof \
                     + self.n_face_dof + self.n_bubble_dof

    def _setup_bubble_dofs(self):
        """
        Creates DOF information for so called element, cell or bubble DOFs
        - the only DOFs used in DG
        n_dof = n_cells * n_el_nod
        remap optional remapping between cells
        dofs is mapping between dofs and cells

        :return: n_dof, remap, dofs
        """
        self.n_cell = self.region.get_n_cells(self.is_surface)
        n_dof = self.n_cell * self.n_el_nod

        dofs = nm.arange(n_dof, dtype=nm.int32)\
                   .reshape(self.n_cell, self.n_el_nod)
        remap = nm.arange(self.n_cell)
        self.econn = dofs
        self.dofs2cells = nm.repeat(nm.arange(self.n_cell), self.n_el_nod)

        return n_dof, remap, dofs

    def _setup_shape(self):
        """
        What is shape used for and what it really means.
        Does it represent shape of the problem?
        """
        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape

    def _setup_geometry(self):
        """
        Setup the field region geometry.
        Somehow
        """
        # get_gel extracts the highest dimension geometry from self.region
        self.gel = get_gel(self.region)

    def _setup_connectivity(self):
        """
        Forces self.domain.mesh to build necessary conductivities
        so they are available in self.get_nbrhd_dofs
        """
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim, self.dim)
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim - 1, self.dim)
        self.region.domain.mesh.cmesh.setup_connectivity(self.dim, self.dim - 1)

    def get_coor(self, nods=None):
        """
        Returns coors for matching nodes
        # TODO revise DG_EPBC and EPBC matching
        :param nods: if None use all nodes
        :return: coors on surface
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
        # use coors of facet QPs?
        coors += eps * nm.repeat(nm.arange(self.n_el_nod),
                                 len(nm.unique(cells)))[:, None]
        return coors

    def clear_facet_qp_base(self):
        """
        Clears facet_qp_base cache
        """
        self.facet_bf = {}
        self.facet_qp = None
        self.facet_whs = None

    def _transform_qps_to_facets(self, qps, geo_name):
        """
        Transforms points given in qps to all facets of the reference element
        with geometry geo_name.
        :param qps: qps corresponding to facet dimension to be transformed
        :param geo_name: element type
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
            # tqps = nm.zeros(nm.shape(qps) + (4, 3,))
            raise NotImplementedError("Geometry {} not supported, yet"
                                      .format(geo_name))
        elif geo_name == "3_8":
            # tqps = nm.zeros(nm.shape(qps) + (8, 3,))
            raise NotImplementedError("Geometry {} not supported, yet"
                                      .format(geo_name))
        else:
            raise NotImplementedError("Geometry {} not supported, yet"
                                      .format(geo_name))
        return tqps

    def get_facet_qp(self):
        """
        Returns quadrature points on all facets of the reference element in
        array of shape (n_qp, 1 , n_el_facets, dim)

        Returns
        -------
        qps : array
            quadrature points
        weights : array
            Still needs to be transformed to actual facets!
        """

        if self.dim == 1:
            facet_qps = self._transform_qps_to_facets(nm.zeros((1, 1)), "1_2")
            weights = nm.ones((1, 1, 1))
        else:
            qps, weights = self.integral.get_qp(cell_facet_gel_name[self.gel.name])
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
        Returns values of base in facets quadrature points, data shape is a bit
        crazy right now:
            (number of qps, 1, n_el_facets, 1, n_el_nod)
        end for derivatine:
            (1, number of qps, (dim,) * derivative, n_el_facets, 1, n_el_nod)

        Parameters
        ----------
        derivative: truthy of integer
        base_only: do not return weights

        Returns
        -------
        facet_bf : array
            values of basis functions in facet qps
        weights : array, optionally
            weights of qps
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

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
            If None clear all.
        """
        if region is None:
            self.facet_neighbour_index = {}
        else:
            self.facet_neighbour_index.pop(region.name)

    def get_facet_neighbor_idx(self, region=None, eq_map=None):
        """
        Returns index of cell neighbours sharing facet, along with local index
        of the facet within neighbour, also treats periodic boundary conditions
        i.e. plugs correct neighbours for cell on periodic boundary.
        Where there are no neighbours specified puts -1  instead of neighbour
        and facet id

        Cashes neighbour index in self.facet_neighbours

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
            Main region, must contain cells.
        eq_map :
            eq_map from state variable containing information on
            EPBC and DG EPBC.

        Returns
        -------
        facet_neighbours : array
             Shape is
                (n_cell, n_el_facet, 2),
             first value is index of the neighbouring cell,
             the second is index of the facet in said nb. cell.
        """
        if region is None or eq_map is None:
            # HOTFIX enabling limiter to obtain connectivity data without
            # knowing eq_map or region
            if self.region.name in self.facet_neighbour_index:
                return self.facet_neighbour_index[self.region.name]
            else:
                raise ValueError("No facet neighbour mapping for main region {}".
                                 format(self.region.name) +
                                 "cached yet, call with region and eq_map first")

        if region.name in self.facet_neighbour_index:
            return self.facet_neighbour_index[region.name]

        dim, n_cell, n_el_facets = self.get_region_info(region)

        cmesh = region.domain.mesh.cmesh
        cells = region.cells

        facet_neighbours = nm.zeros((n_cell, n_el_facets, 2), dtype=nm.int32)

        c2fi, c2fo = cmesh.get_incident(dim - 1, cells, dim, ret_offsets=True)

        for ic, o1 in enumerate(c2fo[:-1]):  # loop over cells
            o2 = c2fo[ic + 1]

            # get neighbours per facet of the cell
            c2ci, c2co = cmesh.get_incident(dim, c2fi[o1:o2], dim - 1,
                                            ret_offsets=True)
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

        facet_neighbours = \
            self._set_fem_periodic_facet_neighbours(facet_neighbours, eq_map)

        facet_neighbours = \
            self._set_dg_periodic_facet_neighbours(facet_neighbours, eq_map)

        # cache results
        self.facet_neighbour_index[region.name] = facet_neighbours

        return facet_neighbours

    def _set_dg_periodic_facet_neighbours(self, facet_neighbours, eq_map):
        """

        Parameters
        ----------
        facet_neighbours : ndarray
            Shape is
                (n_cell, n_el_facet, 2),
            first value is index of the neighbouring cell
            the second is index of the facet in said nb. cell.

        eq_map :
            must contain dg_ep_bc a List with pairs of slave and master boundary
            cell boundary facet mapping

        Returns
        -------
        facet_neighbours : array
            Updated incidence array.

        """
        # if eq_map.
        # treat DG EPBC - these are definitely preferred
        if eq_map.n_dg_epbc > 0 and self.gel.name not in ["1_2", "2_4", "3_6"]:
            raise ValueError(
                "Periodic boundary conditions not supported " +
                "for geometry {} elements.".format(self.gel.name))

        dg_epbc = eq_map.dg_epbc

        for master_bc2bfi, slave_bc2bfi in dg_epbc:
            # set neighbours of periodic cells to one another
            facet_neighbours[master_bc2bfi[:, 0], master_bc2bfi[:, 1], 0] = \
                slave_bc2bfi[:, 0]
            facet_neighbours[slave_bc2bfi[:, 0], slave_bc2bfi[:, 1], 0] = \
                master_bc2bfi[:, 0]

            # set neighbours facets
            facet_neighbours[slave_bc2bfi[:, 0], slave_bc2bfi[:, 1], 1] = \
                master_bc2bfi[:, 1]
            facet_neighbours[master_bc2bfi[:, 0], master_bc2bfi[:, 1], 1] =\
                slave_bc2bfi[:, 1]

        return facet_neighbours

    def _set_fem_periodic_facet_neighbours(self, facet_neighbours, eq_map):
        """
        TODO maybe remove after DG EPBC revision in self.get_coor
        Parameters
        ----------
        facet_neighbours : ndarray
            Shape is (n_cell, n_el_facet, 2), first value is index of the
            neighbouring cell the second is index of the facet in said nb. cell.

        eq_map :
            eq_map from state variable containing information on
            EPBC and DG EPBC.

        Returns
        -------
        facet_neighbours : ndarray
            Updated incidence array.
        """

        # treat classical FEM EPBCs - we need to correct neighbours
        if eq_map.n_epbc > 0:
            # set neighbours of periodic cells to one another
            mcells = nm.unique(self.dofs2cells[eq_map.master])
            scells = nm.unique(self.dofs2cells[eq_map.slave])
            mcells_facets = nm.array(
                nm.where(facet_neighbours[mcells] == -1))[1, 0]  # facets mcells
            scells_facets = nm.array(
                nm.where(facet_neighbours[scells] == -1))[1, 0]  # facets scells
            # [1, 0]  above, first we need second axis to get axis on which
            # facet indices are stored, second we drop axis with neighbour
            # local facet index,
            #
            # for multiple s/mcells this will have to be
            # something like 1 + 2*nm.arange(len(mcells)) - to skip double
            # entries for -1 tags in neighbours and  neighbour local facet idx

            # set neighbours of mcells to scells
            facet_neighbours[mcells, mcells_facets, 0] = scells
            # set neighbour facets to facets of scell missing neighbour
            facet_neighbours[
                mcells, mcells_facets, 1] = scells_facets
            # we do not need to distinguish EBC and EPBC cells, EBC overwrite
            # EPBC, we only need to fix shapes

            # set neighbours of scells to mcells
            facet_neighbours[scells, scells_facets, 0] = mcells
            # set neighbour facets to facets of mcell missing neighbour0
            facet_neighbours[
                scells, scells_facets, 1] = mcells_facets

        return facet_neighbours

    @staticmethod
    def get_region_info(region):
        """
        Extracts information about region needed in various methods of DGField

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region

        Returns
        -------
            dim, n_cell, n_el_facets
        """
        if not region.has_cells():
            raise ValueError("Region {} has no cells".format(region.name))
        n_cell = region.get_n_cells()
        dim = region.tdim
        gel = get_gel(region)
        n_el_facets = dim + 1 if gel.is_simplex else 2 ** dim
        return dim, n_cell, n_el_facets

    def get_both_facet_state_vals(self, state, region,
                                  derivative=None, reduce_nod=True):
        """
        Computes values of the variable represented by dofs in
        quadrature points located at facets, returns both values -
        inner and outer, along with weights.
        Parameters
        ---------
        state: state variable containing BC info
        region : sfepy.discrete.common.region.Region
        derivative: compute derivative if truthy,
                           compute n-th derivative if a number
        reduce_nod: if False DOES NOT sum nodes into values at QPs

        Returns
        -------
        inner_facet_values (n_cell, n_el_facets, n_qp),
                 outer facet values (n_cell, n_el_facets, n_qp),
                 weights,
                 if derivative is True:
                    inner_facet_values (n_cell, n_el_facets, dim, n_qp),
                    outer_facet values (n_cell, n_el_facets, dim, n_qp)

        """
        if derivative:
            diff = int(derivative)
        else:
            diff = 0
        unreduce_nod = int(not reduce_nod)

        inner_base_vals, outer_base_vals, whs = \
            self.get_both_facet_base_vals(state, region, derivative=derivative)
        dofs = self.unravel_sol(state.data[0])

        n_qp = whs.shape[-1]
        outputs_shape = (self.n_cell, self.n_el_facets) + \
                        (self.n_el_nod,) * unreduce_nod + \
                        (self.dim,) * diff + \
                        (n_qp,)

        inner_facet_vals = nm.zeros(outputs_shape)
        if unreduce_nod:
            inner_facet_vals[:] = nm.einsum('id...,idf...->ifd...',
                                            dofs, inner_base_vals)
        else:
            inner_facet_vals[:] = nm.einsum('id...,id...->i...',
                                            dofs, inner_base_vals)

        per_facet_neighbours = self.get_facet_neighbor_idx(region, state.eq_map)

        outer_facet_vals = nm.zeros(outputs_shape)
        for facet_n in range(self.n_el_facets):
            if unreduce_nod:
                outer_facet_vals[:, facet_n, :] = \
                    nm.einsum('id...,id...->id...',
                              dofs[per_facet_neighbours[:, facet_n, 0]],
                              outer_base_vals[:, :, facet_n])
            else:
                outer_facet_vals[:, facet_n, :] = \
                    nm.einsum('id...,id...->i...',
                              dofs[per_facet_neighbours[:, facet_n, 0]],
                              outer_base_vals[:, :, facet_n])

        boundary_cells = nm.array(nm.where(per_facet_neighbours[:, :, 0] < 0)).T
        outer_facet_vals[boundary_cells[:, 0], boundary_cells[:, 1]] = 0.0
        # TODO detect and print boundary cells without defined BCs
        for ebc, ebc_vals in zip(state.eq_map.dg_ebc.get(diff, []),
                                 state.eq_map.dg_ebc_val.get(diff, [])):
            if unreduce_nod:
                raise NotImplementedError("Unreduced DOFs are not available " +
                                          "for boundary outer facets")
                outer_facet_vals[ebc[:, 0], ebc[:, 1], :] = \
                    nm.einsum("id,id...->id...",
                              ebc_vals, inner_base_vals[0, :, ebc[:, 1]])
            else:
                # FIXME contains quick fix flipping qp order to accomodate for
                #  opposite facet orientation of neighbours
                # this is partially taken care of in get_both_facet_base_vals,
                # but needs to be repeated here
                outer_facet_vals[ebc[:, 0], ebc[:, 1], :] = ebc_vals[:, ::-1]

        # flip outer_facet_vals moved to get_both_facet_base_vals
        return inner_facet_vals, outer_facet_vals, whs

    def get_both_facet_base_vals(self, state, region, derivative=None):
        """
        Returns values of the basis function in quadrature points on facets
        broadcasted to all cells inner to the element as well as outer ones
        along with weights for the qps broadcasted and transformed to elements.

        Contains quick fix to flip facet QPs for right integration order.


        Parameters
        ----------
        state: used to get EPBC info
        region : sfepy.discrete.common.region.Region for connectivity
        derivative: if u need derivative

        Returns
        -------
        outer_facet_base_vals:
        inner_facet_base_vals:
                 shape (n_cell, n_el_nod, n_el_facet, n_qp) or
                       (n_cell, n_el_nod, n_el_facet, dim, n_qp)
                 when derivative is True or 1
        whs: shape (n_cell, n_el_facet, n_qp)
        """
        if derivative:
            diff = int(derivative)
        else:
            diff = 0

        facet_bf, whs = self.get_facet_base(derivative=derivative)
        n_qp = nm.shape(whs)[1]
        facet_vols = self.get_facet_vols(region)
        whs = facet_vols * whs[None, :, :, 0]

        base_shape = (self.n_cell, self.n_el_nod, self.n_el_facets) + \
                     (self.dim,) * diff + \
                     (n_qp,)
        inner_facet_base_vals = nm.zeros(base_shape)
        outer_facet_base_vals = nm.zeros(base_shape)

        if derivative:
            inner_facet_base_vals[:] = facet_bf[0, :, 0, :, :, :].swapaxes(-2, -3).T
        else:
            inner_facet_base_vals[:] = facet_bf[:, 0, :, 0, :].T

        per_facet_neighbours = self.get_facet_neighbor_idx(region, state.eq_map)

        # numpy prepends shape resulting from multiple
        # indexing before remaining shape
        if derivative:
            outer_facet_base_vals[:] = inner_facet_base_vals[0, :, per_facet_neighbours[:, :, 1]].swapaxes(-3, -4)
        else:
            outer_facet_base_vals[:] = inner_facet_base_vals[0, :, per_facet_neighbours[:, :, 1]].swapaxes(-2, -3)

        # FIXME quick fix to flip facet QPs for right integration order
        return inner_facet_base_vals, outer_facet_base_vals[..., ::-1], whs

    def clear_normals_cache(self, region=None):
        """
        Clears normals cache for given region or all regions.

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
            region to clear cache or None to clear all
        """
        if region is None:
            self.normals_cache = {}
        else:
            if isinstance(region, str):
                self.normals_cache.pop(region)
            else:
                self.normals_cache.pop(region.name)

    def get_cell_normals_per_facet(self, region):
        """
        Caches results, use clear_normals_cache to clear the cache.

        Parameters
        ----------
        region: sfepy.discrete.common.region.Region
            Main region, must contain cells.
        Returns
        -------
        normals: array
            normals of facets in array of shape (n_cell, n_el_facets, dim)
        """

        if region.name in self.normals_cache:
            return self.normals_cache[region.name]

        dim, n_cell, n_el_facets = self.get_region_info(region)

        cmesh = region.domain.mesh.cmesh
        normals = cmesh.get_facet_normals()
        normals_out = nm.zeros((n_cell, n_el_facets, dim))

        c2f = cmesh.get_conn(dim, dim - 1)
        for ic, o1 in enumerate(c2f.offsets[:-1]):
            o2 = c2f.offsets[ic + 1]
            for ifal, ifa in enumerate(c2f.indices[o1:o2]):
                normals_out[ic, ifal] = normals[o1 + ifal]

        self.normals_cache[region.name] = normals_out

        return normals_out

    def clear_facet_vols_cache(self, region=None):
        """
        Clears facet volume cache for given region or all regions.

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
            region to clear cache or None to clear all
        """
        if region is None:
            self.facet_vols_cache = {}
        else:
            if isinstance(region, str):
                self.facet_vols_cache.pop(region)
            else:
                self.facet_vols_cache.pop(region.name)

    def get_facet_vols(self, region):
        """
        Caches results, use clear_facet_vols_cache to clear the cache

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region

        Returns
        -------
        vols_out: array
            volumes of the facets by cells shape (n_cell, n_el_facets, 1)
        """

        if region.name in self.facet_vols_cache:
            return self.facet_vols_cache[region.name]

        dim, n_cell, n_el_facets = self.get_region_info(region)

        cmesh = region.domain.mesh.cmesh

        if dim == 1:
            vols = nm.ones((cmesh.num[0], 1))
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

        Parameters
        ----------
        integral: integral used
        integration:
            'volume' is only supported value
        region_name: not used

        Returns
        -------
        data_shape : tuple
        """

        if integration in ('volume',):
            # from FEField.get_data_shape()
            _, weights = integral.get_qp(self.gel.name)
            n_qp = weights.shape[0]

            data_shape = (self.n_cell, n_qp, self.gel.dim, self.n_el_nod)
            # econn.shape[1] == n_el_nod i.e. number nod in element
        else:
            raise NotImplementedError('unsupported integration! (%s)'
                                      % integration)

        return data_shape

    def get_econn(self, conn_type, region, is_trace=False, integration=None):
        """
        Getter for econn
        Parameters
        ----------
        conn_type: string or Struct
            'volume' is only supported
        region : sfepy.discrete.common.region.Region
        is_trace: ignored
        integration: ignored

        Returns
        -------
        connectivity
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
        This is called in create_adof_conns(conn_info, var_indx=None,
                                                active_only=True, verbose=True)
        for each variable but has no effect.
        """
        # placeholder, what is this used for?

        # dct = info.dc_type.type

        self.info = info
        self.is_trace = is_trace

    def get_dofs_in_region(self, region, merge=True):
        """
        Return indices of DOFs that belong to the given region.

        Not Used in BC treatment

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
        merge: bool
            merge dof tuple into one numpy array, default True

        Returns
        -------
        dofs : array
        """

        dofs = []
        if region.has_cells():  # main region or its part
            els = nm.ravel(self.bubble_remap[region.cells])
            eldofs = self.bubble_dofs[els[els >= 0]]
            dofs.append(eldofs)
        else:
            # return indices of cells adjacent to boundary facets
            dim = self.dim
            cmesh = region.domain.mesh.cmesh
            bc_cells = cmesh.get_incident(dim, region.facets, dim - 1)
            bc_dofs = self.bubble_dofs[bc_cells]
            dofs.append(bc_dofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def get_bc_facet_idx(self, region):
        """

        Caches results in self.boundary_facet_local_idx

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
            surface region defining BCs

        Returns
        -------
        bc2bfi : array
            index of cells on boundary along with corresponding facets
        """

        if region.name in self.boundary_facet_local_idx:
            return self.boundary_facet_local_idx[region.name]

        bc2bfi = region.get_facet_indices()
        self.boundary_facet_local_idx[region.name] = bc2bfi

        return bc2bfi

    def create_mapping(self, region, integral, integration,
                       return_mapping=True):
        """
        Creates and returns mapping

        Parameters
        ----------
        region : sfepy.discrete.common.region.Region
        integral:
        integration: 'volume' is only accepted option
        return_mapping: default True

        Returns
        -------
        mapping
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

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Compute projection of fun into the basis, alternatively set DOFs
        directly to provided value or values either in main volume region
        or in boundary region.

        Parameters
        ----------
        fun: callable, scalar or array corresponding to dofs
        region : sfepy.discrete.common.region.Region
            region to set DOFs on
        dpn: number of dofs per element
        warn: not used
        Returns
        -------
        nods : array
        vals : array
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

        Parameters
        ----------
        fun: callable, scallar or array corresponding to dofs
        region : sfepy.discrete.common.region.Region
        region to set DOFs on
        dpn: number of dofs per element
        warn: not used

        Returns
        -------
        nods
        vals
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
            # this drops redundant axis that is returned by eval_base due to
            # consistency with derivatives

            # left hand, so far only orthogonal basis
            # for legendre base this can be calculated exactly
            # in 1D it is: 1 / (2 * nm.arange(self.n_el_nod) + 1)
            lhs_diag = nm.einsum("q,q...->...", weights, base_vals_qp ** 2)

            rhs_vec = nm.einsum("q,q...,iq...->i...",
                                weights, base_vals_qp, fun(coors))

            vals = (rhs_vec / lhs_diag)

            # plot for 1D
            # from my_utils.visualizer import plot1D_legendre_dofs, reconstruct
            # _legendre_dofs
            # import matplotlib.pyplot as plt
            # plot1D_legendre_dofs(self.domain.mesh.coors, (vals,), fun)
            # ww, xx = reconstruct_legendre_dofs(self.domain.mesh.coors, 1,
            # vals.T[..., None, None])
            # plt.plot(xx, ww[:, 0], label="reconstructed dofs")
            # plt.show()

        return nods, vals

    def set_facet_dofs(self, fun, region, dpn, warn):
        """
        Compute projection of fun onto the basis on facets, alternatively
        set DOFs directly to provided value or values

        Parameters
        ----------
        fun: callable, scalar or array corresponding to dofs
        region : sfepy.discrete.common.region.Region
        region to set DOFs on
        dpn: number of dofs per element
        warn: not used

        Returns
        -------
        nods
        vals
        """
        raise NotImplementedError(
            "Setting facet DOFs is not supported with DGField, " +
            "use values at qp directly. " +
            "This is usually result of using ebc instead of dgebc")

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
            # get facets QPs
            qp, weights = self.get_facet_qp()
            weights = weights[0, :, 0]
            qp = qp[:, 0, :, :]
            # get facets weights ?

            # get coors
            bc2bfi = self.get_bc_facet_idx(region)
            coors = self.mapping.get_physical_qps(qp)

            # get_physical_qps returns data in strange format, swapping
            # some axis and flipping qps order
            bcoors = coors[bc2bfi[:, 1], ::-1, bc2bfi[:, 0], :]

            # get facet basis vals
            base_vals_qp = self.poly_space.eval_base(qp)[:, 0, 0, :]

            # solve for boundary cell DOFs
            bc_val = fun(bcoors)
            # this returns singular matrix - projection on the boundary should
            # be into facet dim space
            #lhs = nm.einsum("q,qd,qc->dc", weights, base_vals_qp, base_vals_qp)
            # inv_lhs = nm.linalg.inv(lhs)
            # rhs_vec = nm.einsum("q,q...,iq...->i...",
            #                       weights, base_vals_qp, bc_val)

        return nods, vals

    def get_bc_facet_values(self, fun, region, ret_coors=False, diff=0):
        """
        Returns values of fun in facet QPs of the region

        Parameter
        ---------
        diff: derivative 0 or 1 supported
        fun: Function value or values to set qps values to
        region : sfepy.discrete.common.region.Region
            boundary region
        ret_coors: default False,
                Return physical coors of qps in shape (n_cell, n_qp, dim).

        Returns
        -------
        vals
            In shape (n_cell,) + (self.dim,) * diff + (n_qp,)
        """
        if region.has_cells():
            raise NotImplementedError(
                "Region {} has cells and can't be used as boundary region".
                format(region))

        # get facets QPs
        qp, weights = self.get_facet_qp()
        weights = weights[0, :, 0]
        qp = qp[:, 0, :, :]
        n_qp = qp.shape[0]
        # get facets weights ?

        # get physical coors
        bc2bfi = self.get_bc_facet_idx(region)
        n_cell = bc2bfi.shape[0]
        coors = self.mapping.get_physical_qps(qp)

        # get_physical_qps returns data in strange format,
        # swapping some axis and flipping qps order
        # to get coors in shape (n_facet, n_qp, n_cell, dim)
        if len(coors.shape) == 3:
            coors = coors[:, None, :, :]  # add axis for qps when it is missing
            coors = coors.swapaxes(0, 2)
        bcoors = coors[bc2bfi[:, 1], ::-1, bc2bfi[:, 0], :]
        diff_shape = (self.dim,) * diff
        output_shape = (n_cell,) + diff_shape + (n_qp,)
        vals = nm.zeros(output_shape)
        # we do not need last axis of coors, values are scalars

        if nm.isscalar(fun):
            if sum(diff_shape) > 1:
                output("Warning: Setting gradient of shape {} in region {} " +
                       "with scalar value {}"
                              .format(diff_shape, region.name, fun))
            vals[:] = fun

        elif isinstance(fun, nm.ndarray):
            try:
                vals[:] = fun[:, None]
            except ValueError:
                raise ValueError("Provided values of shape {} could not be used"
                                 + " to set BC qps of shape {} in region {}"
                    .format(fun.shape, vals.shape, region.name))

        elif callable(fun):
            # get boundary values
            vals[:] = fun(bcoors)

        if ret_coors:
            return bcoors, vals
        return vals

    def get_nodal_values(self, dofs, region, ref_nodes=None):
        """
        Computes nodal representation of the DOFs

        Parameter
        ---------
        dofs : array
            dofs to transform to nodes
        region : ignored
        ref_nodes:
            reference node to use instead of default qps

        Returns
        -------
        nodes
        nodal_vals
        """
        if ref_nodes is None:
            # TODO poly_space should provide special nodes
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
        Converts the DOFs corresponding to the field to a dictionary of
        output data usable by Mesh.write().

        Puts DOFs into vairables u0 ... un, where n = approx_order and marks
        them for writing as cell data.

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
        # TODO revise output dictionary structure and naming
        res = {}
        udofs = self.unravel_sol(dofs)

        if self.dim == 1:
            for i in range(self.n_el_nod):
                res["u_modal{}".format(i)] = Struct(mode="cell",
                                                    data=udofs[:, i, None, None])
        else:
            unravel = get_unraveler(self.n_el_nod, self.n_cell)
            res["u_modal_cell_nodes"] = Struct(mode="cell_nodes",
                                               data=unravel(dofs)[..., 0],
                                               interpolation_scheme=self.poly_space.get_interpol_scheme())
        return res