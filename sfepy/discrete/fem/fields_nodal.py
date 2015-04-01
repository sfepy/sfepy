"""
Notes
-----

Important attributes of continuous (order > 0) :class:`Field` and
:class:`SurfaceField` instances:

- `vertex_remap` : `econn[:, :n_vertex] = vertex_remap[conn]`
- `vertex_remap_i` : `conn = vertex_remap_i[econn[:, :n_vertex]]`

where `conn` is the mesh vertex connectivity, `econn` is the
region-local field connectivity.
"""
import time
import numpy as nm

from sfepy.base.base import output, assert_
import fea
from sfepy.discrete.fem.utils import prepare_remap
from sfepy.discrete.common.dof_info import expand_nodes_to_dofs
from sfepy.discrete.fem.global_interp import get_ref_coors
from sfepy.discrete.fem.facets import get_facet_dof_permutations
from sfepy.discrete.fem.fields_base import (FEField, VolumeField, SurfaceField,
                                            H1Mixin)
from sfepy.discrete.fem.extmods.bases import evaluate_in_rc

class H1NodalMixin(H1Mixin):

    def _setup_facet_orientations(self):
        order = self.approx_order
        self.node_desc = self.interp.describe_nodes()

        edge_nodes = self.node_desc.edge_nodes
        if edge_nodes is not None:
            n_fp = self.gel.edges.shape[1]
            self.edge_dof_perms = get_facet_dof_permutations(n_fp, order)

        face_nodes = self.node_desc.face_nodes
        if face_nodes is not None:
            n_fp = self.gel.faces.shape[1]
            self.face_dof_perms = get_facet_dof_permutations(n_fp, order)

    def _setup_edge_dofs(self):
        """
        Setup edge DOF connectivity.
        """
        if self.node_desc.edge is None:
            return 0, None, None

        return self._setup_facet_dofs(1, self.node_desc.edge,
                                      self.edge_dof_perms,
                                      self.n_vertex_dof)

    def _setup_face_dofs(self):
        """
        Setup face DOF connectivity.
        """
        if self.node_desc.face is None:
            return 0, None, None

        return self._setup_facet_dofs(self.domain.shape.tdim - 1,
                                      self.node_desc.face,
                                      self.face_dof_perms,
                                      self.n_vertex_dof + self.n_edge_dof)

    def _setup_facet_dofs(self, dim, facet_desc, facet_perms, offset):
        """
        Helper function to setup facet DOF connectivity, works for both
        edges and faces.
        """
        facet_desc = nm.array(facet_desc)
        n_dof_per_facet = facet_desc.shape[1]

        cmesh = self.domain.cmesh

        facets = self.region.entities[dim]
        ii = nm.arange(facets.shape[0], dtype=nm.int32)
        all_dofs = offset + expand_nodes_to_dofs(ii, n_dof_per_facet)

        # Prepare global facet id remapping to field-local numbering.
        remap = prepare_remap(facets, cmesh.num[dim])

        cconn = self.region.domain.cmesh.get_conn(self.region.tdim, dim)
        offs = cconn.offsets

        n_f = self.gel.edges.shape[0] if dim == 1 else self.gel.faces.shape[0]

        oris = cmesh.get_orientations(dim)

        gcells = self.region.get_cells()
        n_el = gcells.shape[0]

        indices = cconn.indices[offs[gcells[0]]:offs[gcells[-1]+1]]
        facets_of_cells = remap[indices]

        ori = oris[offs[gcells[0]]:offs[gcells[-1]+1]]
        perms = facet_perms[ori]

        # Define global facet dof numbers.
        gdofs = offset + expand_nodes_to_dofs(facets_of_cells,
                                              n_dof_per_facet)

        # Elements of facets.
        iel = nm.arange(n_el, dtype=nm.int32).repeat(n_f)
        ies = nm.tile(nm.arange(n_f, dtype=nm.int32), n_el)

        # DOF columns in econn for each facet.
        iep = facet_desc[ies]

        iaux = nm.arange(gdofs.shape[0], dtype=nm.int32)
        self.ap.econn[iel[:, None], iep] = gdofs[iaux[:, None], perms]

        n_dof = n_dof_per_facet * facets.shape[0]
        assert_(n_dof == nm.prod(all_dofs.shape))

        return n_dof, all_dofs, remap

    def _setup_bubble_dofs(self):
        """
        Setup bubble DOF connectivity.
        """
        if self.node_desc.bubble is None:
            return 0, None, None

        offset = self.n_vertex_dof + self.n_edge_dof + self.n_face_dof
        n_dof_per_cell = self.node_desc.bubble.shape[0]

        ap = self.ap
        ii = self.region.get_cells()
        n_cell = ii.shape[0]
        n_dof = n_dof_per_cell * n_cell

        all_dofs = nm.arange(offset, offset + n_dof, dtype=nm.int32)
        all_dofs.shape = (n_cell, n_dof_per_cell)
        iep = self.node_desc.bubble[0]
        ap.econn[:,iep:] = all_dofs

        return n_dof, all_dofs

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Set the values of DOFs in a given region using a function of space
        coordinates or value `fun`.
        """
        if region is None:
            region = self.region

        if dpn is None:
            dpn = self.n_components

        aux = self.get_dofs_in_region(region, clean=True, warn=warn)
        nods = nm.unique(nm.hstack(aux))

        if callable(fun):
            vals = fun(self.get_coor(nods))

        elif nm.isscalar(fun):
            vals = nm.repeat([fun], nods.shape[0] * dpn)

        elif isinstance(fun, nm.ndarray):
            assert_(len(fun) == dpn)
            vals = nm.repeat(fun, nods.shape[0])

        else:
            raise ValueError('unknown function/value type! (%s)' % type(fun))

        return nods, vals

    def evaluate_at(self, coors, source_vals, strategy='kdtree',
                    close_limit=0.1, cache=None, ret_cells=False,
                    ret_status=False, ret_ref_coors=False, verbose=False):
        """
        Evaluate source DOF values corresponding to the field in the given
        coordinates using the field interpolation.

        Parameters
        ----------
        coors : array
            The coordinates the source values should be interpolated into.
        source_vals : array
            The source DOF values corresponding to the field.
        strategy : str, optional
            The strategy for finding the elements that contain the
            coordinates. Only 'kdtree' is supported for the moment.
        close_limit : float, optional
            The maximum limit distance of a point from the closest
            element allowed for extrapolation.
        cache : Struct, optional
            To speed up a sequence of evaluations, the field mesh, the inverse
            connectivity of the field mesh and the KDTree instance can
            be cached as `cache.mesh`, `cache.offsets`, `cache.iconn` and
            `cache.kdtree`. Optionally, the cache can also contain the
            reference element coordinates as `cache.ref_coors`,
            `cache.cells` and `cache.status`, if the evaluation occurs
            in the same coordinates repeatedly. In that case the KDTree
            related data are ignored.
        ret_cells : bool, optional
            If True, return also the cell indices the coordinates are in.
        ret_status : bool, optional
            If True, return also the status for each point: 0 is
            success, 1 is extrapolation within `close_limit`, 2 is
            extrapolation outside `close_limit`, 3 is failure.
        ret_ref_coors : bool, optional
            If True, return also the found reference element coordinates.
        verbose : bool
            If False, reduce verbosity.

        Returns
        -------
        vals : array
            The interpolated values.
        cells : array
            The cell indices, if `ret_cells` or `ret_status` are True.
        status : array
            The status, if `ret_status` is True.
        """
        output('evaluating in %d points...' % coors.shape[0], verbose=verbose)

        ref_coors, cells, status = get_ref_coors(self, coors,
                                                 strategy=strategy,
                                                 close_limit=close_limit,
                                                 cache=cache,
                                                 verbose=verbose)

        tt = time.clock()
        vertex_coorss, nodess, orders, mtx_is = [], [], [], []
        conns = []
        for ap in self.aps.itervalues():
            ps = ap.interp.poly_spaces['v']

            vertex_coorss.append(ps.geometry.coors)
            nodess.append(ps.nodes)
            mtx_is.append(ps.get_mtx_i())

            orders.append(ps.order)
            conns.append(ap.econn)

        orders = nm.array(orders, dtype=nm.int32)

        # Interpolate to the reference coordinates.
        vals = nm.empty((coors.shape[0], source_vals.shape[1]),
                        dtype=source_vals.dtype)

        evaluate_in_rc(vals, ref_coors, cells, status, source_vals,
                       conns, vertex_coorss, nodess, orders, mtx_is,
                       1e-15)
        output('interpolation: %f s' % (time.clock()-tt),verbose=verbose)

        output('...done',verbose=verbose)

        if ret_ref_coors:
            return vals, ref_coors, cells, status

        elif ret_status:
            return vals, cells, status

        elif ret_cells:
            return vals, cells

        else:
            return vals

class H1NodalVolumeField(H1NodalMixin, VolumeField):
    family_name = 'volume_H1_lagrange'

    def interp_v_vals_to_n_vals(self, vec):
        """
        Interpolate a function defined by vertex DOF values using the FE
        geometry base (P1 or Q1) into the extra nodes, i.e. define the
        extra DOF values.
        """
        if not self.node_desc.has_extra_nodes():
            enod_vol_val = vec.copy()

        else:
            dim = vec.shape[1]
            enod_vol_val = nm.zeros((self.n_nod, dim), dtype=nm.float64)
            for ig, ap in self.aps.iteritems():
                group = self.domain.groups[ig]
                econn = ap.econn

                coors = ap.interp.poly_spaces['v'].node_coors

                ginterp = ap.interp.gel.interp
                bf = ginterp.poly_spaces['v'].eval_base(coors)
                bf = bf[:,0,:].copy()

                evec = nm.dot(bf, vec[group.conn])
                enod_vol_val[econn] = nm.swapaxes(evec, 0, 1)

        return enod_vol_val

class H1DiscontinuousField(H1NodalMixin, VolumeField):
    family_name = 'volume_H1_lagrange_discontinuous'

    def _setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.DiscontinuousApproximation(name, self.interp,
                                                self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def _setup_global_base( self ):
        """
        Setup global DOF/base function indices and connectivity of the field.
        """
        self._setup_facet_orientations()

        self._init_econn()

        n_dof = 0
        all_dofs = {}
        remaps = {}
        for ig, ap in self.aps.iteritems():
            ii = self.region.get_cells(ig)
            nd = nm.prod(ap.econn.shape)

            group = self.domain.groups[ig]
            remaps[ig] = prepare_remap(ii, group.shape.n_el)

            aux = nm.arange(n_dof, n_dof + nd, dtype=nm.int32)
            aux.shape = ap.econn.shape

            ap.econn[:] = aux
            all_dofs[ig] = aux

            n_dof += nd

        self.n_nod = n_dof

        self.n_bubble_dof = n_dof
        self.bubble_dofs = all_dofs
        self.bubble_remaps = remaps

        self.n_vertex_dof = self.n_edge_dof = self.n_face_dof = 0

        self._setup_esurface()

    def extend_dofs(self, dofs, fill_value=None):
        """
        Extend DOFs to the whole domain using the `fill_value`, or the
        smallest value in `dofs` if `fill_value` is None.
        """
        if self.approx_order != 0:
            dofs = self.average_to_vertices(dofs)

        new_dofs = FEField.extend_dofs(self, dofs)

        return new_dofs

    def remove_extra_dofs(self, dofs):
        """
        Remove DOFs defined in higher order nodes (order > 1).
        """
        if self.approx_order != 0:
            dofs = self.average_to_vertices(dofs)

        new_dofs = FEField.remove_extra_dofs(self, dofs)

        return new_dofs

    def average_to_vertices(self, dofs):
        """
        Average DOFs of the discontinuous field into the field region
        vertices.
        """
        data_qp, integral = self.interp_to_qp(dofs)
        vertex_dofs = self.average_qp_to_vertices(data_qp, integral)

        return vertex_dofs

class H1NodalSurfaceField(H1NodalMixin, SurfaceField):
    """
    A field defined on a surface region.
    """
    family_name = 'surface_H1_lagrange'

    def interp_v_vals_to_n_vals(self, vec):
        """
        Interpolate a function defined by vertex DOF values using the FE
        surface geometry base (P1 or Q1) into the extra nodes, i.e. define the
        extra DOF values.
        """
        if not self.node_desc.has_extra_nodes():
            enod_vol_val = vec.copy()

        else:
            msg = 'surface nodal fields do not support higher order nodes yet!'
            raise NotImplementedError(msg)

        return enod_vol_val
