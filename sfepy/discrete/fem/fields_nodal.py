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
import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.discrete.integrals import Integral
from sfepy.discrete.fem.utils import prepare_remap
from sfepy.discrete.common.dof_info import expand_nodes_to_dofs
from sfepy.discrete.common.mappings import get_physical_qps
from sfepy.discrete.fem.facets import get_facet_dof_permutations
from sfepy.discrete.fem.fields_base import FEField, H1Mixin

class GlobalNodalLikeBasis(Struct):

    def _setup_facet_orientations(self):
        order = self.approx_order
        self.node_desc = self.poly_space.describe_nodes()

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

        cmesh = self.cmesh

        facets = self.region.entities[dim]
        ii = nm.arange(facets.shape[0], dtype=nm.int32)
        all_dofs = offset + expand_nodes_to_dofs(ii, n_dof_per_facet)

        # Prepare global facet id remapping to field-local numbering.
        remap = prepare_remap(facets, cmesh.num[dim])

        cconn = cmesh.get_conn(self.region.tdim, dim)
        offs = cconn.offsets

        n_f = self.gel.edges.shape[0] if dim == 1 else self.gel.faces.shape[0]

        oris = cmesh.get_orientations(dim)

        gcells = self.region.get_cells()
        n_el = gcells.shape[0]

        # Elements of facets.
        iel = nm.arange(n_el, dtype=nm.int32).repeat(n_f)
        ies = nm.tile(nm.arange(n_f, dtype=nm.int32), n_el)

        aux = offs[gcells][:, None] + ies.reshape((n_el, n_f))

        indices = cconn.indices[aux]
        facets_of_cells = remap[indices].ravel()

        ori = oris[aux].ravel()
        perms = facet_perms[ori]

        # Define global facet dof numbers.
        gdofs = offset + expand_nodes_to_dofs(facets_of_cells,
                                              n_dof_per_facet)

        # DOF columns in econn for each facet.
        iep = facet_desc[ies]

        iaux = nm.arange(gdofs.shape[0], dtype=nm.int32)
        self.econn[iel[:, None], iep] = gdofs[iaux[:, None], perms]

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

        ii = self.region.get_cells()
        remap = prepare_remap(ii, self.cmesh.n_el)

        n_cell = ii.shape[0]
        n_dof = n_dof_per_cell * n_cell

        all_dofs = nm.arange(offset, offset + n_dof, dtype=nm.int32)
        all_dofs.shape = (n_cell, n_dof_per_cell)
        iep = self.node_desc.bubble[0]
        self.econn[:,iep:] = all_dofs

        return n_dof, all_dofs, remap

    def get_surface_basis(self, region):
        """
        Get basis for projections to region's facets.

        Notes
        -----
        Cannot be uses for all fields because IGA does not support surface
        mappings.
        """
        order = self.approx_order

        integral = Integral('i', order=2*order)
        geo, mapping = self.get_mapping(region, integral, 'surface')
        pqps = get_physical_qps(region, integral)
        qps = pqps.values.reshape(pqps.shape)

        bfs = nm.broadcast_to(
            geo.bf[..., 0, :],
            (qps.shape[0], qps.shape[1], geo.bf.shape[3]),
        )

        return qps, bfs, geo.det[..., 0]

class H1NodalMixin(H1Mixin, GlobalNodalLikeBasis):

    def _substitute_dofs(self, subs):
        """
        Perform facet DOF substitutions according to `subs`.

        Modifies `self.econn` in-place.
        """
        if self.gel.name == '2_4':
            ef = self.efaces

            for ii, sub in enumerate(subs):
                # 2_4 edges always in opposite orientation.
                ee = ef[sub[1]].copy()
                ee[0], ee[1] = ee[1], ee[0] # Swap vertex DOFs.
                ee[2:] = ee[-1:1:-1] # Swap edge DOFs.

                master = self.econn[sub[0], ee]
                self.econn[sub[2], ef[sub[3]]] = master
                self.econn[sub[4], ef[sub[5]]] = master

        elif self.gel.name == '3_8':
            def _sort4(p):
                key = 0

                if (p[0] < p[1]): key += 1
                if (p[0] < p[2]): key += 2
                if (p[1] < p[2]): key += 4
                if (p[0] < p[3]): key += 8
                if (p[1] < p[3]): key += 16
                if (p[2] < p[3]): key += 32

                return key

            if subs[0] is not None:
                ef = self.efaces
                epf = self.gel.get_edges_per_face()
                nde = self.node_desc.edge
                ndf = self.node_desc.face
                gedges = self.gel.edges
                gfaces = self.gel.faces

                for ii, sub in enumerate(subs[0]):
                    master = self.econn[sub[0]]
                    fmaster = master[ef[sub[1]]]
                    lmaster = fmaster.tolist()

                    for ic in range(4):
                        ia, ib = 2 + 2 * ic, 2 + 2 * ic + 1
                        cell = self.econn[sub[ia]]

                        # Corner vertex is always the first for faces 0, 1, 2.
                        iv = cell[ef[sub[ib]][0]]
                        i0 = lmaster.index(iv)
                        for ik in range(4):
                            cell[ef[sub[ib]][ik]] = fmaster[:4][i0 - ik]

                        # Treat edge DOFs.
                        if nde is not None:
                            sedges = epf[sub[ib]]
                            medges = epf[sub[1]]
                            for ie, sedge in enumerate(sedges):
                                iim =  i0 - 1 - ie
                                ies = nde[sedge]
                                medge = medges[iim]
                                iem = nde[medge]

                                vm = master[gedges[medge]][0]
                                vs = cell[gedges[sedge]][0]
                                if vm == vs:
                                    cell[ies] = self.econn[sub[0], iem]

                                else:
                                    cell[ies] = self.econn[sub[0], iem[::-1]]

                        # Treat face DOFs.
                        if ndf is not None:
                            new_ori = _sort4(cell[gfaces[sub[ib]]])
                            smaster = nm.sort(master[ndf[sub[1]]])
                            aux = self.face_dof_perms[new_ori]
                            cell[ndf[sub[ib]]] = smaster[aux]

            if subs[1] is not None:
                ef = self.eedges
                for ii, sub in enumerate(subs[1]):
                    master = self.econn[sub[0]]

                    me = master[gedges[sub[1]]]
                    for ic in range(2):
                        ia, ib = 2 + 2 * ic, 2 + 2 * ic + 1
                        cell = self.econn[sub[ia]]
                        ce = cell[gedges[sub[ib]]]

                        if (me[0] == ce[0]) or (me[1] == ce[1]):
                            cell[ef[sub[ib]]] = master[ef[sub[1]]]

                        else:
                            ee = ef[sub[1]].copy()
                            ee[0], ee[1] = ee[1], ee[0] # Swap vertex DOFs.
                            ee[2:] = ee[-1:1:-1] # Swap edge DOFs.
                            cell[ef[sub[ib]]] = master[ee]

        else:
            raise ValueError('unsupported reference element type! (%s)'
                             % self.gel.name)

    def _eval_basis_transform(self, subs):
        """
        """
        from sfepy.discrete import Integral
        from sfepy.discrete.fem import Mesh, FEDomain, Field

        transform = nm.tile(nm.eye(self.econn.shape[1]),
                            (self.econn.shape[0], 1, 1))
        if subs is None:
            return transform

        gel = self.gel
        ao = self.approx_order

        conn = [gel.conn]
        mesh = Mesh.from_data('a', gel.coors, None, [conn], [nm.array([0])],
                              [gel.name])
        cdomain = FEDomain('d', mesh)
        comega = cdomain.create_region('Omega', 'all')
        rcfield = Field.from_args('rc', self.dtype, 1, comega, approx_order=ao)

        fdomain = cdomain.refine()
        fomega = fdomain.create_region('Omega', 'all')
        rffield = Field.from_args('rf', self.dtype, 1, fomega, approx_order=ao)

        def assign_transform(transform, bf, subs, ef):
            if not len(subs): return

            n_sub = (subs.shape[1] - 2) // 2

            for ii, sub in enumerate(subs):
                for ij in range(n_sub):
                    ik = 2 * (ij + 1)

                    fface = ef[sub[ik+1]]

                    mtx = transform[sub[ik]]
                    ix, iy = nm.meshgrid(fface, fface)

                    cbf = bf[iy, 0, ix]

                    mtx[ix, iy] = cbf

        fcoors = rffield.get_coor()

        coors = fcoors[rffield.econn[0]]
        integral = Integral('i', coors=coors, weights=nm.ones_like(coors[:, 0]))

        rcfield.clear_qp_basis()
        bf = rcfield.eval_basis('v', False, integral)

        if gel.name == '2_4':
            fsubs = subs
            esubs = None

            assign_transform(transform, bf, fsubs, rffield.efaces)

        else:
            fsubs = subs[0]
            esubs = subs[1]

            assign_transform(transform, bf, fsubs, rffield.efaces)
            if esubs is not None:
                assign_transform(transform, bf, esubs, rffield.eedges)

        assert_((nm.abs(transform.sum(1) - 1.0) < 1e-15).all())

        return transform

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Set the values of DOFs in a given `region` using a function of space
        coordinates or value `fun`.
        """
        if region is None:
            region = self.region

        if dpn is None:
            dpn = self.n_components

        aux = self.get_dofs_in_region(region)
        nods = nm.unique(aux)

        if callable(fun):
            coors = self.get_coor(nods)
            vals = nm.asarray(fun(coors))
            if (vals.ndim > 1) and (vals.shape != (len(coors), dpn)):
                raise ValueError('The projected function return value should be'
                                 ' (n_point, dpn) == %s, instead of %s!'
                                 % ((len(coors), dpn), vals.shape))

        elif nm.isscalar(fun):
            vals = nm.repeat([fun], nods.shape[0] * dpn)

        elif isinstance(fun, nm.ndarray):
            try:
                assert_(len(fun) == dpn)

            except (TypeError, ValueError):
                msg = ('wrong array value shape for setting'
                       ' DOFs of "%s" field!'
                       ' (shape %s should be %s)'
                       % (self.name, fun.shape, (dpn,)))
                raise ValueError(msg)

            vals = nm.repeat(fun, nods.shape[0])

        else:
            raise ValueError('unknown function/value type! (%s)' % type(fun))

        vals.shape = (len(nods), -1)

        return nods, vals

    def create_basis_context(self):
        """
        Create the context required for evaluating the field basis.
        """
        ps = self.poly_space
        gps = self.gel.poly_space

        mesh = self.create_mesh(extra_nodes=False)

        ctx = ps.create_context(None, 0, 1e-15, 100, 1e-8,
                                tdim=mesh.cmesh.tdim)
        geo_ctx = gps.create_context(mesh.cmesh, 0, 1e-15, 100, 1e-8)

        ctx.geo_ctx = geo_ctx

        return ctx

class H1NodalVolumeField(H1NodalMixin, FEField):
    """
    Lagrange basis nodal approximation.
    """
    family_name = 'volume_H1_lagrange'

    def interp_v_vals_to_n_vals(self, vec):
        """
        Interpolate a function defined by vertex DOF values using the FE
        geometry basis (P1 or Q1) into the extra nodes, i.e. define the
        extra DOF values.
        """
        if not self.node_desc.has_extra_nodes():
            enod_vol_val = vec.copy()

        else:
            dim = vec.shape[1]
            enod_vol_val = nm.zeros((self.n_nod, dim), dtype=nm.float64)

            coors = self.poly_space.node_coors

            bf = self.gel.poly_space.eval_basis(coors)
            bf = bf[:,0,:].copy()

            conn = self.econn[:, :self.gel.n_vertex]

            evec = nm.dot(bf, vec[conn])
            enod_vol_val[self.econn] = nm.swapaxes(evec, 0, 1)

        return enod_vol_val

class H1SNodalVolumeField(H1NodalVolumeField):
    """
    Lagrange basis nodal serendipity approximation with order <= 3.
    """
    family_name = 'volume_H1_serendipity'

    def create_basis_context(self):
        """
        Create the context required for evaluating the field basis.
        """
        # Hack for tests to pass - the reference coordinates are determined
        # from vertices only - we can use the Lagrange basis context for the
        # moment. The true context for Field.evaluate_at() is not implemented.
        gps = self.gel.poly_space
        mesh = self.create_mesh(extra_nodes=False)

        ctx = geo_ctx = gps.create_context(self.cmesh, 0, 1e-15, 100, 1e-8)
        ctx.geo_ctx = geo_ctx

        return ctx

class H1SEMVolumeField(H1NodalVolumeField):
    """
    Spectral element method approximation.

    Uses the Lagrange basis with Legendre-Gauss-Lobatto nodes and quadrature.
    """
    family_name = 'volume_H1_sem'

    def create_basis_context(self):
        """
        Create the context required for evaluating the field basis.
        """
        # Hack for tests to pass - the reference coordinates are determined
        # from vertices only - we can use the Lagrange basis context for the
        # moment. The true context for Field.evaluate_at() is not implemented.
        gps = self.gel.poly_space
        mesh = self.create_mesh(extra_nodes=False)

        ctx = geo_ctx = gps.create_context(self.cmesh, 0, 1e-15, 100, 1e-8)
        ctx.geo_ctx = geo_ctx

        return ctx

class H1DiscontinuousField(H1NodalMixin, FEField):
    """
    The C0 constant-per-cell approximation.
    """
    family_name = 'volume_H1_lagrange_discontinuous'

    def _setup_global_basis(self):
        """
        Setup global DOF/basis function indices and connectivity of the field.
        """
        self._setup_facet_orientations()

        self._init_econn()

        ii = self.region.get_cells()
        self.bubble_remap = prepare_remap(ii, self.cmesh.n_el)

        n_dof = nm.prod(self.econn.shape)
        all_dofs = nm.arange(n_dof, dtype=nm.int32)
        all_dofs.shape = self.econn.shape

        self.econn[:] = all_dofs

        self.n_nod = n_dof

        self.n_bubble_dof = n_dof
        self.bubble_dofs = all_dofs

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

class H1NodalSurfaceField(H1NodalMixin, FEField):
    """
    A field defined on a surface region.
    """
    family_name = 'surface_H1_lagrange'

    def interp_v_vals_to_n_vals(self, vec):
        """
        Interpolate a function defined by vertex DOF values using the FE
        surface geometry basis (P1 or Q1) into the extra nodes, i.e. define the
        extra DOF values.
        """
        if not self.node_desc.has_extra_nodes():
            enod_vol_val = vec.copy()

        else:
            msg = 'surface nodal fields do not support higher order nodes yet!'
            raise NotImplementedError(msg)

        return enod_vol_val

class H1SNodalSurfaceField(H1NodalSurfaceField):
    family_name = 'surface_H1_serendipity'

class H1SEMSurfaceField(H1NodalSurfaceField):
    family_name = 'surface_H1_sem'
