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

from sfepy.base.base import output, iter_dict_of_lists, get_default, assert_
from sfepy.base.base import Struct
import fea
from sfepy.fem.mesh import Mesh, make_inverse_connectivity
from sfepy.fem.utils import extend_cell_data, prepare_remap, invert_remap
from sfepy.fem.fe_surface import FESurface
from sfepy.fem.dof_info import expand_nodes_to_dofs
from sfepy.fem.integrals import Integral
from sfepy.fem.extmods.fem import evaluate_at

def parse_approx_order(approx_order):
    """
    Parse the uniform approximation order value (str or int).
    """
    ao_msg = 'unsupported approximation order! (%s)'
    force_bubble = False
    discontinuous = False

    try:
        ao = int(approx_order)
    except ValueError:
        mode = approx_order[-1].lower()
        if mode == 'b':
            ao = int(approx_order[:-1])
            force_bubble = True

        elif mode == 'd':
            ao = int(approx_order[:-1])
            discontinuous = True

        else:
            raise ValueError(ao_msg % approx_order)

    if ao < 0:
        raise ValueError(ao_msg % approx_order)

    elif ao == 0:
        discontinuous = True

    return ao, force_bubble, discontinuous

def create_dof_conn(conn, dpn):
    """Given element a node connectivity, create the dof connectivity."""
    if dpn == 1:
        dc = conn.copy()
    else:
        n_el, n_ep = conn.shape
        n_ed = n_ep * dpn
        dc = nm.empty( (n_el, n_ed), dtype = conn.dtype )
        for ic in range( n_ed ):
            inod = ic / dpn
            idof = ic % dpn
##                    iloc = ic
            iloc = n_ep * idof + inod # Hack: For DBD order.
            dc[:,iloc] = dpn * conn[:,inod] + idof

    return dc

def fields_from_conf(conf, regions):
    fields = {}
    for key, val in conf.iteritems():
        field = Field.from_conf(val, regions)
        fields[field.name] = field

    return fields

def setup_extra_data(conn_info):
    """
    Setup extra data required for non-volume integration.
    """
    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        ## print key, ii
        ## print info
        for var in info.all_vars:
            field = var.get_field()
            field.setup_extra_data(info.ps_tg, info, info.is_trace)

def setup_dof_conns(conn_info, dof_conns=None,
                    make_virtual=False, verbose=True):
    """
    Dof connectivity key:
        (field.name, var.n_components, region.name, type, ig)
    """
    if verbose:
        output('setting up dof connectivities...')
        tt = time.clock()

    dof_conns = get_default(dof_conns, {})

    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        ## print key, ii
        ## print info

        if info.primary is not None:
            var = info.primary
            field = var.get_field()
            field.setup_extra_data(info.ps_tg, info, info.is_trace)
            field.setup_dof_conns(dof_conns, var.n_components,
                                  info.dc_type, info.get_region())

        if info.has_virtual and not info.is_trace:
            # This is needed regardless make_virtual.
            var = info.virtual
            field = var.get_field()
            field.setup_extra_data(info.v_tg, info, False)
            field.setup_dof_conns(dof_conns, var.n_components,
                                  info.dc_type,
                                  info.get_region(can_trace=False))

    ## print dof_conns
    ## pause()
    if verbose:
        output('...done in %.2f s' % (time.clock() - tt))

    return dof_conns

##
# 14.07.2006, c
class Field( Struct ):
    """
    Finite element field.

    Notes
    -----

    - Region can span over several groups -> different Aproximation
      instances
    - interps and hence node_descs are per region (must have single
      geometry!)
    - no two interps can be in a same group -> no two aps (with
      different regions) can be in a same group -> aps can be uniquely
      indexed with ig
    """

    @staticmethod
    def from_conf(conf, regions):
        """To refactor... very hackish now."""
        space = conf.get_default_attr('space', 'H1')
        poly_space_base = conf.get_default_attr('poly_space_base', 'lagrange')

        approx_order = parse_approx_order(conf.approx_order)
        ao, force_bubble, discontinuous = approx_order

        if isinstance(conf.region, tuple):
            region_name, kind = conf.region
            region = regions[region_name]
            if kind == 'surface':
                obj = SurfaceField(conf.name, conf.dtype, conf.shape, region,
                                   space=space,
                                   poly_space_base=poly_space_base,
                                   approx_order=approx_order[:2])

            else:
                raise ValueError('unknown field kind! (%s)', kind)

        else:
            if discontinuous:
                cls = DiscontinuousField

            else:
                cls = Field

            obj = cls(conf.name, conf.dtype, conf.shape, regions[conf.region],
                      space=space,
                      poly_space_base=poly_space_base,
                      approx_order=approx_order[:2])

        return obj

    def __init__(self, name, dtype, shape, region,
                 space='H1', poly_space_base='lagrange', approx_order=1):
        """Create a Field.

        Parameters
        ----------
        name : str
            Object name.
        dtype : numpy.dtype
            Field data type: float64 or complex128.
        shape : int/tuple/str
            Field shape: 1 or (1,) or 'scalar', space dimension (2, or
            (2,) or 3 or (3,)) or 'vector'. The field shape determines
            the shape of the FE base functions and can be different from
            a FieldVariable instance shape. (TODO)

        region : Region
            The region where the field is defined.
        space : str
            The function space name.
        poly_space_base : str
            The name of polynomial space base.
        approx_order : int/str
            FE approximation order, e.g. 0, 1, 2, '1B' (1 with bubble).

        Notes
        -----
        Assumes one cell type for the whole region!
        """
        if isinstance(shape, str):
            try:
                shape = {'scalar' : (1,),
                         'vector' : (region.domain.shape.dim,)}[shape]
            except KeyError:
                raise ValueError('unsupported field shape! (%s)', shape)

        elif isinstance(shape, int):
            shape = (shape,)

        Struct.__init__(self,
                        name = name,
                        dtype = dtype,
                        shape = shape,
                        region = region,
                        space = space,
                        poly_space_base = poly_space_base)
        self.domain = self.region.domain
        self.igs = self.region.igs

        self.clear_dof_conns()

        self.set_approx_order(approx_order)
        self.setup_geometry()

        self.create_interpolant()
        self.setup_approximations()
        self.setup_global_base()
        self.setup_coors()
        self.clear_mappings(clear_all=True)

    def set_approx_order(self, approx_order):
        """
        Set a uniform approximation order.
        """
        if isinstance(approx_order, tuple):
            self.approx_order = approx_order[0]
            self.force_bubble = approx_order[1]

        else:
            self.approx_order = approx_order
            self.force_bubble = False

    def setup_geometry(self):
        """
        Setup the field region geometry.
        """
        self.gel = self.domain.groups[self.region.igs[0]].gel

    def create_interpolant(self):
        name = '%s_%d%s' % (self.gel.name, self.approx_order,
                            'B' * self.force_bubble)
        self.interp = fea.Interpolant(name, self.gel, self.approx_order,
                                      self.force_bubble)

    def setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.Approximation(name, self.interp, self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def get_true_order(self):
        """
        Get the true approximation order depending on the reference
        element geometry.

        For example, for P1 (linear) approximation the true order is 1,
        while for Q1 (bilinear) approximation in 2D the true order is 2.
        """
        gel = self.gel
        if (gel.dim + 1) == gel.n_vertex:
            order = self.approx_order

        else:
            order = gel.dim * self.approx_order

        return order

    def setup_global_base(self):
        """
        Setup global DOF/base function indices and connectivity of the field.
        """
        self.setup_facet_orientations()

        self.init_econn()

        self.n_vertex_dof, self.vertex_remap = self.setup_vertex_dofs()
        self.vertex_remap_i = invert_remap(self.vertex_remap)

        aux = self.setup_edge_dofs()
        self.n_edge_dof, self.edge_dofs, self.edge_remap = aux

        aux = self.setup_face_dofs()
        self.n_face_dof, self.face_dofs, self.face_remap = aux

        aux = self.setup_bubble_dofs()
        self.n_bubble_dof, self.bubble_dofs, self.bubble_remaps = aux

        self.n_nod = self.n_vertex_dof + self.n_edge_dof \
                     + self.n_face_dof + self.n_bubble_dof

        self.setup_esurface()

    def setup_facet_orientations(self):
        self.node_desc = self.interp.describe_nodes()

        edge_nodes = self.node_desc.edge_nodes
        if edge_nodes is not None:
            ed = self.domain.ed
            self.edge_dof_perms = ed.get_facet_dof_permutations(edge_nodes)

        face_nodes = self.node_desc.face_nodes
        if face_nodes is not None:
            fa = self.domain.fa
            self.face_dof_perms = fa.get_facet_dof_permutations(face_nodes)

    def init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        for ig, ap in self.aps.iteritems():
            n_ep = ap.n_ep['v']
            n_cell = self.region.get_n_cells(ig)
            ap.econn = nm.zeros((n_cell, n_ep), nm.int32)

    def setup_vertex_dofs(self):
        """
        Setup vertex DOF connectivity.
        """
        if self.node_desc.vertex is None:
            return 0, None

        region = self.region

        remap = prepare_remap(region.all_vertices, region.n_v_max)
        n_dof = region.all_vertices.shape[0]

        ##
        # Remap vertex node connectivity to field-local numbering.
        for ig, ap in self.aps.iteritems():
            group = self.domain.groups[ig]
            offset = group.shape.n_ep
            cells = region.get_cells(ig)
            ap.econn[:,:offset] = remap[group.conn[cells]]

        return n_dof, remap

    def _setup_facet_dofs(self, facets, facet_desc, facet_perms,
                          get_facets, offset):
        """
        Helper function to setup facet DOF connectivity, works for both
        edges and faces.
        """
        facet_desc = nm.array(facet_desc)
        n_dof_per_facet = facet_desc.shape[1]

        # Prepare global facet id remapping to field-local numbering.
        uids = []
        for ig, ap in self.aps.iteritems():
            ii = get_facets(ig)
            uid_i = facets.uid_i[ii]

            uids.append(uid_i)

        uids = nm.unique(nm.concatenate(uids))
        n_uid = uids.shape[0]
        lids = nm.arange(n_uid, dtype=nm.int32)
        remap = prepare_remap(uids, facets.n_unique)

        all_dofs = offset + expand_nodes_to_dofs(lids, n_dof_per_facet)

        for ig, ap in self.aps.iteritems():
            ori = facets.oris[ig]
            perms = facet_perms[ig][ori]

            ii = get_facets(ig)
            g_uid = facets.uid_i[ii]
            uid = remap[g_uid]

            # Define global facet dof numbers.
            gdofs = offset + expand_nodes_to_dofs(uid, n_dof_per_facet)

            # Elements of facets.
            iel = facets.indices[ii, 1]

            ies = facets.indices[ii, 2]
            # DOF columns in econn for each facet.
            iep = facet_desc[ies]

            iaux = nm.arange(gdofs.shape[0], dtype=nm.int32)
            ap.econn[iel[:, None], iep] = gdofs[iaux[:, None], perms]

        n_dof = n_dof_per_facet * n_uid
        assert_(n_dof == nm.prod(all_dofs.shape))

        return n_dof, all_dofs, remap

    def setup_edge_dofs(self):
        """
        Setup edge DOF connectivity.
        """
        if self.node_desc.edge is None:
            return 0, None, None

        return self._setup_facet_dofs(self.domain.ed,
                                      self.node_desc.edge,
                                      self.edge_dof_perms,
                                      self.region.get_edges,
                                      self.n_vertex_dof)

    def setup_face_dofs(self):
        """
        Setup face DOF connectivity.
        """
        if self.node_desc.face is None:
            return 0, None, None

        return self._setup_facet_dofs(self.domain.fa,
                                      self.node_desc.face,
                                      self.face_dof_perms,
                                      self.region.get_faces,
                                      self.n_vertex_dof + self.n_edge_dof)

    def setup_bubble_dofs(self):
        """
        Setup bubble DOF connectivity.
        """
        if self.node_desc.bubble is None:
            return 0, None, None

        offset = self.n_vertex_dof + self.n_edge_dof + self.n_face_dof
        n_dof = 0
        n_dof_per_cell = self.node_desc.bubble.shape[0]
        all_dofs = {}
        remaps = {}
        for ig, ap in self.aps.iteritems():
            ii = self.region.get_cells(ig)
            n_cell = ii.shape[0]
            nd = n_dof_per_cell * n_cell

            group = self.domain.groups[ig]
            remaps[ig] = prepare_remap(ii, group.shape.n_el)

            aux = nm.arange(offset + n_dof, offset + n_dof + nd,
                            dtype=nm.int32)
            aux.shape = (n_cell, n_dof_per_cell)
            iep = self.node_desc.bubble[0]
            ap.econn[:,iep:] = aux
            all_dofs[ig] = aux

            n_dof += nd

        return n_dof, all_dofs, remaps

    def setup_esurface(self):
        """
        Setup extended surface entities (edges in 2D, faces in 3D),
        i.e. indices of surface entities into the extended connectivity.
        """
        node_desc = self.node_desc

        for ig, ap in self.aps.iteritems():
            gel = ap.interp.gel
            ap.efaces = gel.get_surface_entities().copy()

            nd = node_desc.edge
            if nd is not None:
                efs = []
                for eof in gel.get_edges_per_face():
                    efs.append(nm.concatenate([nd[ie] for ie in eof]))
                efs = nm.array(efs).squeeze()

                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
                ap.efaces = nm.hstack((ap.efaces, efs))

            efs = node_desc.face
            if efs is not None:
                efs = nm.array(efs).squeeze()

                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
                ap.efaces = nm.hstack((ap.efaces, efs))

    def setup_coors(self, coors=None):
        """
        Setup coordinates of field nodes.
        """
        mesh = self.domain.mesh
        self.coors = nm.empty((self.n_nod, mesh.dim), nm.float64)

        if coors is None:
            coors = mesh.coors

        # Mesh vertex nodes.
        if self.n_vertex_dof:
            indx = self.region.all_vertices
            self.coors[:self.n_vertex_dof] = coors[indx]

        for ig, ap in self.aps.iteritems():
            ap.eval_extra_coor(self.coors, coors)

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if geometry != None:
            geometry_flag = 'surface' in geometry
        else:
            geometry_flag = False

        if (dct == 'surface') or (geometry_flag):
            reg = info.get_region()
            # Calls reg.select_cells_of_surface(reset=False)...
            self.domain.create_surface_group(reg)
            self.setup_surface_data(reg)

        elif dct == 'edge':
            raise NotImplementedError('dof connectivity type %s' % dct)

        elif dct == 'point':
            self.setup_point_data(self, info.region)

        elif dct not in ('volume', 'scalar'):
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def setup_surface_data(self, region):
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue
            if region.name not in ap.surface_data:
                ap.setup_surface_data(region)

    def setup_point_data(self, field, region):
        # Point data only in the first group to avoid multiple
        # assembling of nodes on group boundaries.
        ap = self.aps[self.igs[0]]
        if region.name not in ap.point_data:
            ap.setup_point_data(field, region)

    def get_vertices(self):
        """
        Return indices of vertices belonging to the field region.
        """
        return self.region.all_vertices

    def get_dofs_in_region(self, region, merge=False, clean=False,
                           warn=False, igs=None):
        """
        Return indices of DOFs that belong to the given region.
        """
        if igs is None:
            igs = region.igs

        nods = []
        for ig in self.igs:
            if not ig in igs:
                nods.append(None)
                continue

            nn = self.get_dofs_in_region_group(region, ig)
            nods.append(nn)

        if merge:
            nods = [nn for nn in nods if nn is not None]
            nods = nm.unique(nm.hstack(nods))

        elif clean:
            for nn in nods[:]:
                if nn is None:
                    nods.remove(nn)
                    if warn is not None:
                        output(warn + ('%s' % region.name))

        return nods

    def _get_facet_dofs(self, facets, get_facets, remap, dofs, ig):
        ii = get_facets(ig)
        g_uid = facets.uid_i[ii]
        uid = remap[g_uid]

        return dofs[uid[uid >= 0]].ravel()

    def get_dofs_in_region_group(self, region, ig):
        """
        Return indices of DOFs that belong to the given region and group.
        """
        node_desc = self.node_desc

        dofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.vertex is not None:
            ii = region.get_vertices(ig)
            vdofs = self.vertex_remap[ii]
            dofs = nm.concatenate((dofs, vdofs[vdofs >= 0]))

        if node_desc.edge is not None:
            edofs = self._get_facet_dofs(self.domain.ed,
                                         region.get_edges,
                                         self.edge_remap,
                                         self.edge_dofs, ig)
            dofs = nm.concatenate((dofs, edofs))

        if node_desc.face is not None:
            fdofs = self._get_facet_dofs(self.domain.fa,
                                         region.get_faces,
                                         self.face_remap,
                                         self.face_dofs, ig)
            dofs = nm.concatenate((dofs, fdofs))

        if (node_desc.bubble is not None) and region.can_cells:
            ii = region.get_cells(ig)
            group_els = self.bubble_remaps[ig][ii]
            bdofs = self.bubble_dofs[ig][group_els[group_els >= 0]].ravel()
            dofs = nm.concatenate((dofs, bdofs))

        return dofs

    def extend_dofs(self, dofs, fill_value=None):
        """
        Extend DOFs to the whole domain using the `fill_value`, or the
        smallest value in `dofs` if `fill_value` is None.
        """
        if fill_value is None:
            if dofs.shape[1]: # Vector.
                fill_value = nm.amin(nm.abs(dofs))
            else: # Scalar.
                fill_value = nm.amin(dofs)

        if self.approx_order != 0:
            indx = self.get_vertices()

            n_nod = self.domain.shape.n_nod
            new_dofs = nm.empty((n_nod, dofs.shape[1]), dtype=self.dtype)
            new_dofs.fill(fill_value)
            new_dofs[indx] = dofs[:indx.size]

        else:
            new_dofs = extend_cell_data(dofs, self.domain, self.region,
                                        val=fill_value)

        return new_dofs

    def remove_extra_dofs(self, dofs):
        """
        Remove DOFs defined in higher order nodes (order > 1).
        """
        if self.approx_order != 0:
            new_dofs = dofs[:self.n_vertex_dof]

        else:
            new_dofs = dofs

        return new_dofs

    def get_output_approx_order(self):
        """
        Get the approximation order used in the output file.
        """
        return min(self.approx_order, 1)

    def clear_dof_conns(self):
        self.dof_conns = {}

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region):
        """Setup dof connectivities of various kinds as needed by terms."""
        dct = dc_type.type

        ##
        # Expand nodes into dofs.
        can_point = True
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue

            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig)
            if key in dof_conns:
                self.dof_conns[key] = dof_conns[key]

                if dct == 'point':
                    can_point = False
                continue

            if dct == 'volume':
                dc = create_dof_conn(ap.econn, dpn)
                self.dof_conns[key] = dc

            elif dct == 'surface':
                sd = ap.surface_data[region_name]
                dc = create_dof_conn(sd.econn, dpn)
                self.dof_conns[key] = dc

            elif dct == 'edge':
                raise NotImplementedError('dof connectivity type %s' % dct)

            elif dct == 'point':
                if can_point:
                    # Point data only in the first group to avoid multiple
                    # assembling of nodes on group boundaries.
                    conn = ap.point_data[region_name]
                    dc = create_dof_conn(conn, dpn)
                    self.dof_conns[key] = dc
                    can_point = False

            else:
                raise ValueError('unknown dof connectivity type! (%s)' % dct)

        dof_conns.update(self.dof_conns)

    def create_mesh(self, extra_nodes=True):
        """
        Create a mesh from the field region, optionally including the field
        extra nodes.
        """
        mesh = self.domain.mesh

        if self.approx_order != 0:
            conns, mat_ids, descs = [], [], []
            for ig, ap in self.aps.iteritems():
                group = self.domain.groups[ig]
                if extra_nodes:
                    conn = ap.econn
                else:
                    offset = group.shape.n_ep
                    conn = ap.econn[:,:offset]
                conns.append(conn)
                mat_ids.append(mesh.mat_ids[ig])
                descs.append(mesh.descs[ig])

            if extra_nodes:
                coors = self.coors

            else:
                coors = self.coors[:self.n_vertex_dof]

            mesh = Mesh.from_data(self.name, coors, None, conns,
                                  mat_ids, descs)

        return mesh

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
        integral = Integral('i', order=self.get_true_order())

        data_qp = []
        for ig, ap in self.aps.iteritems():
            bf = ap.get_base('v', False, integral)
            bf = bf[:,0,:].copy()

            vals = nm.dot(bf, dofs[ap.econn])
            vals = nm.swapaxes(vals, 0, 1)
            vals.shape = vals.shape + (1,)

            data_qp.append(vals)

        data_qp = nm.concatenate(data_qp, axis=0)

        return data_qp, integral

    def average_qp_to_vertices(self, data_qp, integral):
        """
        Average data given in quadrature points in region elements into
        region vertices.

        .. math::
           u_n = \sum_e (u_{e,avg} * volume_e) / \sum_e volume_e
               = \sum_e \int_{volume_e} u / \sum volume_e
        """
        domain = self.domain
        if domain.shape.n_el != data_qp.shape[0]:
            msg = 'incomatible shape! (%d == %d)' % (domain.shape.n_el,
                                                     data_qp.shape[0])
            raise ValueError(msg)

        n_vertex = domain.shape.n_nod
        dim = data_qp.shape[2]

        nod_vol = nm.zeros((n_vertex,), dtype=nm.float64)
        data_vertex = nm.zeros((n_vertex, dim), dtype=nm.float64)
        for ig, ap in self.aps.iteritems():
            vg = ap.describe_geometry(self, 'volume', ap.region, integral)

            volume = nm.squeeze(vg.variable(2))
            iels = ap.region.cells[ig]

            data_e = nm.zeros((volume.shape[0], 1, dim, 1), dtype=nm.float64)
            vg.integrate(data_e, data_qp[iels])

            ir = nm.arange(dim, dtype=nm.int32)

            conn = domain.groups[ig].conn
            for ii, cc in enumerate(conn):
                # Assumes unique nodes in cc!
                ind2, ind1 = nm.meshgrid(ir, cc)
                data_vertex[ind1,ind2] += data_e[iels[ii],0,:,0]
                nod_vol[cc] += volume[ii]
        data_vertex /= nod_vol[:,nm.newaxis]

        return data_vertex

    def get_coor(self, nods=None):
        """
        Get coordinates of the field nodes.

        Parameters
        ----------
        nods : array, optional
           The indices of the required nodes. If not given, the
           coordinates of all the nodes are returned.
        """
        if nods is None:
            return self.coors
        else:
            return self.coors[nods]

    def clear_mappings(self, clear_all=False):
        """
        Clear current reference mappings.
        """
        self.mappings = {}
        if clear_all:
            self.mappings0 = {}

    def save_mappings(self):
        """
        Save current reference mappings to `mappings0` attribute.
        """
        self.mappings0 = self.mappings.copy()

    def create_mapping(self, ig, region, integral, integration):
        """
        Create a new reference mapping.
        """
        ap = self.aps[ig]

        out = ap.describe_geometry(self, integration, region, integral)
        return out

    def get_mapping(self, ig, region, integral, integration):
        """
        For given region, integral and integration type, get a reference
        mapping, i.e. jacobians, element volumes and base function
        derivatives for Volume-type geometries, and jacobians, normals
        and base function derivatives for Surface-type geometries
        corresponding to the field approximation.

        The mappings are cached in the field instance in `mappings`
        attribute.  The mappings can be saved to `mappings0` using
        `Field.save_mappings`.
        """
        ap = self.aps[ig]
        # Share full group mappings.
        if region.shape[ig].n_vertex == self.domain.groups[ig].shape.n_vertex:
            region_name = ig

        else:
            region_name = region.name

        key = (integral.get_key(), region_name, ig, integration)

        # out is (geo, mapping) tuple.
        out = self.mappings.get(key, None)
        if out is None:
            out = ap.describe_geometry(self, integration, region, integral,
                                       return_mapping=True)
            self.mappings[key] = out

        return out[0]

    def describe_geometry(self, geometry_type, ig, region,
                          term_region, integral):
        """
        For give approximation, compute jacobians, element volumes and
        base function derivatives for Volume-type geometries, and
        jacobians, normals and base function derivatives for
        Surface-type geometries.

        Usually, region is term_region. Only if is_trace is True, then region
        is the mirror region and term_region is the true term region.
        """
        ap = self.aps[ig]

        geo = ap.describe_geometry(self, geometry_type, region, integral)
        return geo

    def update_geometry(self, regions, geometries):
        for geom_key, geom in geometries.iteritems():
            iname, region_name, geometry_type, ap_name = geom_key
            ap = self.aps_by_name[ap_name]
            geom = ap.describe_geometry(self, geometry_type,
                                        regions[region_name], geom.integral)
            geometries[geom_key] = geom


    def evaluate_at(self, coors, source_vals, strategy='kdtree',
                    close_limit=0.1, cache=None, ret_cells=False,
                    ret_status=False):
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
            To speed up a sequence of evaluations, the inverse
            connectivity of the field mesh and the KDTree instance can
            be cached as `cache.offsets`, `cache.iconn` and
            `cache.kdtree`.
        ret_cells : bool, optional
            If True, return also the cell indices the coordinates are in.
        ret_status : bool, optional
            If True, return also the status for each point: 0 is
            success, 1 is extrapolation within `close_limit`, 2 is
            extrapolation outside `close_limit`, 3 is failure.

        Returns
        -------
        vals : array
            The interpolated values.
        cells : array
            The cell indices, if `ret_cells` or `ret_status` are True.
        status : array
            The status, if `ret_status` is True.
        """
        mesh = self.create_mesh()
        scoors = mesh.coors

        output('interpolating from %d nodes to %d nodes...' % (scoors.shape[0],
                                                               coors.shape[0]))

        if cache is None:
            offsets, iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                                       ret_offsets=True)
        else:
            offsets, iconn = cache.offsets, cache.iconn

        if strategy == 'kdtree':
            if cache is None:
                from scipy.spatial import cKDTree as KDTree
                ## from scipy.spatial import KDTree

                tt = time.clock()
                ctree = KDTree(scoors)
                output('ctree: %f s' % (time.clock()-tt))

            else:
                ctree = cache.ctree

            tt = time.clock()

            vals = nm.empty((coors.shape[0], source_vals.shape[1]),
                            dtype=source_vals.dtype)
            cells = nm.empty((coors.shape[0], 2), dtype=nm.int32)
            status = nm.empty((coors.shape[0],), dtype=nm.int32)

            ics = ctree.query(coors)[1]
            ics = nm.asarray(ics, dtype=nm.int32)

            vertex_coorss, nodess, orders, mtx_is = [], [], [], []
            conns, conns0 = [], []
            for ap in self.aps.itervalues():
                ps = ap.interp.poly_spaces['v']
                if ps.order == 0:
                    # Use geometry element space and connectivity to locate an
                    # element a point is in.
                    ps = ap.interp.gel.interp.poly_spaces['v']
                    assert_(ps.order == 1)

                    orders.append(0) # Important!
                    conn = self.domain.groups[ap.ig].conn
                    conns.append(conn)

                else:
                    orders.append(ps.order)
                    conns.append(ap.econn)

                vertex_coorss.append(ps.geometry.coors)
                nodess.append(ps.nodes)
                mtx_is.append(ps.get_mtx_i())

                # Always the true connectivity for extracting source values.
                conns0.append(ap.econn)

            orders = nm.array(orders, dtype=nm.int32)

            evaluate_at(vals, cells, status, coors, source_vals,
                        ics, offsets, iconn,
                        scoors, conns0, conns,
                        vertex_coorss, nodess, orders, mtx_is,
                        1, close_limit, 1e-15, 100, 1e-8)

            output('interpolator: %f s' % (time.clock()-tt))

        elif strategy == 'crawl':
            raise NotImplementedError

        else:
            raise ValueError('unknown search strategy! (%s)' % strategy)

        output('...done')

        if ret_status:
            return vals, cells, status

        elif ret_cells:
            return vals, cells

        else:
            return vals

class DiscontinuousField(Field):

    def setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.DiscontinuousApproximation(name, self.interp,
                                                self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def setup_global_base( self ):
        """
        Setup global DOF/base function indices and connectivity of the field.
        """
        self.setup_facet_orientations()

        self.init_econn()

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

        self.setup_esurface()

    def setup_coors(self, coors=None):
        """
        Setup coordinates of field nodes.
        """
        mesh = self.domain.mesh
        self.coors = nm.empty((self.n_nod, mesh.dim), nm.float64)

        if coors is None:
            coors = mesh.coors

        for ig, ap in self.aps.iteritems():
            ap.eval_extra_coor(self.coors, coors)

    def extend_dofs(self, dofs, fill_value=None):
        """
        Extend DOFs to the whole domain using the `fill_value`, or the
        smallest value in `dofs` if `fill_value` is None.
        """
        if self.approx_order != 0:
            dofs = self.average_to_vertices(dofs)

        new_dofs = Field.extend_dofs(self, dofs)

        return new_dofs

    def remove_extra_dofs(self, dofs):
        """
        Remove DOFs defined in higher order nodes (order > 1).
        """
        if self.approx_order != 0:
            dofs = self.average_to_vertices(dofs)

        new_dofs = Field.remove_extra_dofs(self, dofs)

        return new_dofs

    def average_to_vertices(self, dofs):
        """
        Average DOFs of the discontinuous field into the field region
        vertices.
        """
        data_qp, integral = self.interp_to_qp(dofs)
        vertex_dofs = self.average_qp_to_vertices(data_qp, integral)

        return vertex_dofs

class SurfaceField(Field):
    """
    A field defined on a surface region.
    """

    def __init__(self, name, dtype, shape, region,
                 space='H1', poly_space_base='lagrange', approx_order=1):
        region.setup_face_indices()

        Field.__init__(self, name, dtype, shape, region,
                       space=space, poly_space_base=poly_space_base,
                       approx_order=approx_order)

    def setup_geometry(self):
        """
        Setup the field region geometry.
        """
        self.gel = self.domain.groups[self.region.igs[0]].gel.surface_facet
        if self.gel is None:
            raise ValueError('element group has no surface!')

    def create_interpolant(self):
        name = '%s_%d%s' % (self.gel.name, self.approx_order,
                            'B' * self.force_bubble)
        self.interp = fea.SurfaceInterpolant(name, self.gel, self.approx_order,
                                             self.force_bubble)

    def setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.SurfaceApproximation(name, self.interp, self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if dct != 'surface':
            msg = "dof connectivity type must be 'surface'! (%s)" % dct
            raise ValueError(msg)

        reg = info.get_region()
        reg.select_cells_of_surface(reset=False)

        for ig, ap in self.aps.iteritems():
            if ig not in reg.igs: continue

            if reg.name not in ap.surface_data:
                # Defined in setup_vertex_dofs()
                msg = 'no surface data of surface field! (%s)' % reg.name
                raise ValueError(msg)

    def init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        for ig, ap in self.aps.iteritems():
            n_ep = ap.n_ep['v']
            n_cell = self.region.get_n_cells(ig, True)
            ap.econn = nm.zeros((n_cell, n_ep), nm.int32)

    def setup_vertex_dofs(self):
        """
        Setup vertex DOF connectivity.
        """
        if self.node_desc.vertex is None:
            return 0, None

        region = self.region

        remap = prepare_remap(region.all_vertices, region.n_v_max)
        n_dof = region.all_vertices.shape[0]

        ##
        # Remap vertex node connectivity to field-local numbering.
        for ig, ap in self.aps.iteritems():
            group = self.domain.groups[ig]
            faces = group.gel.get_surface_entities()
            aux = FESurface('aux', region, faces, group.conn, ig)
            ap.econn[:,:aux.n_fp] = aux.leconn
            ap.surface_data[region.name] = aux

        return n_dof, remap

    def setup_bubble_dofs(self):
        """
        Setup bubble DOF connectivity.
        """
        return 0, None, None

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region):
        """Setup dof connectivities of various kinds as needed by terms."""
        dct = dc_type.type

        if dct != 'surface':
            msg = "dof connectivity type must be 'surface'! (%s)" % dct
            raise ValueError(msg)

        ##
        # Expand nodes into dofs.
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue

            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig)
            if key in dof_conns:
                self.dof_conns[key] = dof_conns[key]
                continue

            sd = ap.surface_data[region_name]
            dc = create_dof_conn(sd.leconn, dpn)
            self.dof_conns[key] = dc

        dof_conns.update(self.dof_conns)
