import numpy as nm

from sfepy.base.base import assert_
from sfepy.discrete.fem.utils import prepare_remap, prepare_translate
from sfepy.discrete.common.dof_info import expand_nodes_to_dofs
from sfepy.discrete.fem.fields_base import FEField, H1Mixin

class H1HierarchicVolumeField(H1Mixin, FEField):
    """
    Hierarchical basis approximation with Lobatto polynomials.
    """
    family_name = 'volume_H1_lobatto'

    def _init_econn(self):
        """
        Initialize the extended DOF connectivity and facet orientation array.
        """
        FEField._init_econn(self)

        self.ori = nm.zeros_like(self.econn)

    def _setup_facet_orientations(self):
        self.node_desc = self.poly_space.describe_nodes()

    def _setup_edge_dofs(self):
        """
        Setup edge DOF connectivity.
        """
        if self.node_desc.edge is None:
            return 0, None, None

        return self._setup_facet_dofs(1,
                                      self.node_desc.edge,
                                      self.n_vertex_dof)

    def _setup_face_dofs(self):
        """
        Setup face DOF connectivity.
        """
        if self.node_desc.face is None:
            return 0, None, None

        return self._setup_facet_dofs(self.domain.shape.tdim - 1,
                                      self.node_desc.face,
                                      self.n_vertex_dof + self.n_edge_dof)

    def _setup_facet_dofs(self, dim, facet_desc, offset):
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
        n_fp = 2 if dim == 1 else self.gel.surface_facet.n_vertex

        oris = cmesh.get_orientations(dim)

        gcells = self.region.get_cells()
        n_el = gcells.shape[0]

        # Elements of facets.
        iel = nm.arange(n_el, dtype=nm.int32).repeat(n_f)
        ies = nm.tile(nm.arange(n_f, dtype=nm.int32), n_el)

        aux = offs[gcells][:, None] + ies.reshape((n_el, n_f))

        indices = cconn.indices[aux]
        facets_of_cells = remap[indices].ravel()

        # Define global facet dof numbers.
        gdofs = offset + expand_nodes_to_dofs(facets_of_cells,
                                              n_dof_per_facet)

        # DOF columns in econn for each facet (repeating same values for
        # each element.
        iep = facet_desc[ies]

        self.econn[iel[:, None], iep] = gdofs

        ori = oris[aux].ravel()

        if (n_fp == 2) and (self.gel.name in ['2_4', '3_8']):
            tp_edges = self.gel.edges
            ecs = self.gel.coors[tp_edges]
            # True = positive, False = negative edge orientation w.r.t.
            # reference tensor product axes.
            tp_edge_ori = (nm.diff(ecs, axis=1).sum(axis=2) > 0).squeeze()
            aux = nm.tile(tp_edge_ori, n_el)
            ori = nm.where(aux, ori, 1 - ori)

        if n_fp == 2: # Edges.
            # ori == 1 means the basis has to be multiplied by -1.
            ps = self.poly_space
            orders = ps.node_orders
            eori = nm.repeat(ori[:, None], n_dof_per_facet, 1)
            eoo = orders[iep] % 2 # Odd orders.
            self.ori[iel[:, None], iep] = eori * eoo

        elif n_fp == 3: # Triangular faces.
            raise NotImplementedError

        else: # Quadrilateral faces.
            # ori encoding in 3 bits:
            # 0: axis swap, 1: axis 1 sign, 2: axis 2 sign
            # 0 = + or False, 1 = - or True
            # 63 -> 000 = 0
            #  0 -> 001 = 1
            # 30 -> 010 = 2
            # 33 -> 011 = 3
            # 11 -> 100 = 4
            #  7 -> 101 = 5
            # 52 -> 110 = 6
            # 56 -> 111 = 7
            # Special cases:
            # Both orders same and even -> 000
            # Both orders same and odd -> 0??
            # Bits 1, 2 are multiplied by (swapped) axial order % 2.
            new = nm.repeat(nm.arange(8, dtype=nm.int32), 3)
            translate = prepare_translate([31, 59, 63,
                                           0, 1, 4,
                                           22, 30, 62,
                                           32, 33, 41,
                                           11, 15, 43,
                                           3, 6, 7,
                                           20, 52, 60,
                                           48, 56, 57], new)
            ori = translate[ori]
            eori = nm.repeat(ori[:, None], n_dof_per_facet, 1)

            ps = self.poly_space
            orders = ps.face_axes_nodes[iep - ps.face_indx[0]]
            eoo = orders % 2
            eoo0, eoo1 = eoo[..., 0], eoo[..., 1]

            i0 = nm.where(eori < 4)
            i1 = nm.where(eori >= 4)

            eori[i0] = nm.bitwise_and(eori[i0], 2*eoo0[i0] + 5)
            eori[i0] = nm.bitwise_and(eori[i0], eoo1[i0] + 6)

            eori[i1] = nm.bitwise_and(eori[i1], eoo0[i1] + 6)
            eori[i1] = nm.bitwise_and(eori[i1], 2*eoo1[i1] + 5)

            self.ori[iel[:, None], iep] = eori

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

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Set the values of DOFs in a given `region` using a function of space
        coordinates or value `fun`.
        """
        if region is None:
            region = self.region

        if dpn is None:
            dpn = self.n_components

        # Hack - use only vertex DOFs.
        gnods = self.get_dofs_in_region(region, merge=False)
        nods = nm.concatenate(gnods)
        n_dof = dpn * nods.shape[0]

        if nm.isscalar(fun):
            vals = nm.zeros(n_dof, dtype=nm.dtype(type(fun)))
            vals[:gnods[0].shape[0] * dpn] = fun

        elif callable(fun):
            coors = self.get_coor(gnods[0])
            vv = nm.asarray(fun(coors))
            if (vv.ndim > 1) and (vv.shape != (len(coors), dpn)):
                raise ValueError('The projected function return value should be'
                                 ' (n_point, dpn) == %s, instead of %s!'
                                 % ((len(coors), dpn), vv.shape))

            vals = nm.zeros(n_dof, dtype=vv.dtype)
            vals[:gnods[0].shape[0] * dpn] = vv.ravel()

        else:
            raise ValueError('unknown function/value type! (%s)' % type(fun))

        nods, indx = nm.unique(nods, return_index=True)
        ii = (nm.tile(dpn * indx, dpn)
              + nm.tile(nm.arange(dpn, dtype=nm.int32), indx.shape[0]))
        vals = vals[ii]

        vals.shape = (len(nods), -1)

        return nods, vals

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
