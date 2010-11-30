import numpy as nm

from sfepy.base.base import Struct, Container, OneTypeList, assert_
from sfepy.fem.mappings import VolumeMapping, SurfaceMapping
from poly_spaces import PolySpace
from fe_surface import FESurface
import extmods.meshutils as mu

def set_mesh_coors( domain, fields, geometries, coors, update_state = False ):
    domain.mesh.coors = coors.copy()
    if update_state:
        for field in fields.itervalues():
            field.setup_coors()
            field.aps.update_geometry( field, domain.regions, geometries )


##
# 04.08.2005, c
def _interp_to_faces( vertex_vals, bfs, faces ):
    dim = vertex_vals.shape[1]
    n_face = faces.shape[0]
    n_qp = bfs.shape[0]
    
    faces_vals = nm.zeros( (n_face, n_qp, dim), nm.float64 )
    for ii, face in enumerate( faces ):
        vals = vertex_vals[face,:dim]
        faces_vals[ii,:,:] = nm.dot( bfs[:,0,:], vals )

    return( faces_vals )

class Interpolant( Struct ):
    """A simple wrapper around PolySpace."""

    def __init__(self, name, gel, approx_order=1, force_bubble=False):
        self.name = name
        self.gel = gel

        self.poly_spaces = poly_spaces = {}
        poly_spaces['v'] = PolySpace.any_from_args(name, gel, approx_order,
                                                   base='lagrange',
                                                   force_bubble=force_bubble)
        gel = gel.surface_facet
        if gel is not None:
            ps = PolySpace.any_from_args(name, gel, approx_order,
                                         base='lagrange',
                                         force_bubble=False)
            skey = 's%d' % ps.n_nod
            poly_spaces[skey] = ps

    def describe_nodes( self ):
        ps = self.poly_spaces['v']
        node_desc = ps.describe_nodes()

        return node_desc

    ##
    # 16.11.2007, c
    def get_n_nodes( self ):
        nn = {}
        for key, ps in self.poly_spaces.iteritems():
            nn[key] = ps.nodes.shape[0]
        return nn

##
# 18.07.2006, c
class Approximation( Struct ):
    ##
    # 18.07.2006, c
    # 10.10.2006
    # 11.07.2007
    # 17.07.2007
    def __init__(self, name, interp, region, ig, is_surface=False):
        """interp, region are borrowed."""

        self.name = name
        self.interp = interp
        self.region = region
        self.ig = ig
        self.is_surface = is_surface
        self.surface_data = {}
        self.edge_data = {}
        self.point_data = {}
        self.qp_coors = {}
        self.bf = {}
        self.n_ep = self.interp.get_n_nodes()

    ##
    # 19.07.2006, c
    def describe_nodes( self, node_desc ):
        self.node_desc = node_desc

        self.has_extra_edge_nodes = self.has_extra_face_nodes = False
        if self.node_desc.edge is not None:
            self.has_extra_edge_nodes = True
        if self.node_desc.face is not None:
            self.has_extra_face_nodes = True

##         print self.node_desc
##        pause()

    ##
    # 03.10.2005, c
    # 04.10.2005
    # 10.10.2005
    # 26.10.2005
    # 19.07.2006
    # 02.08.2006
    # 04.08.2006
    # 13.02.2007
    # 20.02.2007
    def eval_extra_coor( self, coors, mesh ):
        """Evaluate coordinates of extra nodes."""
        node_offsets = self.node_offsets
        n_nod = nm.sum( node_offsets[1:,1] - node_offsets[1:,0] )
#        print self.name
#        print n_nod
##         pause()

        if n_nod == 0: return

        ##
        # Evaluate geometry interpolation base functions in extra nodes.
        ginterp = self.interp.gel.interp
        ps = self.interp.poly_spaces['v']

        iex = (ps.nts[:,0] > 0).nonzero()[0]
        qp_coors = ps.node_coors[iex,:]

        bf = ginterp.poly_spaces['v'].eval_base(qp_coors)
        bf = bf[:,0,:].copy()

        ##
        # Evaulate extra coordinates with 'bf'.
#        print self.econn.shape, self.sub.n_el

        econn = nm.zeros_like( self.econn[:,iex] ).copy()
        for ii in range( 1, 4 ):
            ix = nm.where( ps.nts[:,0] == ii )[0]
            if not len( ix ): continue
            
##             print ii, ix, iex[0], ix-iex[0]
##             pause()
            econn[:,ix-iex[0]] = self.econn[:,ix]

##         print self.econn
##         print econn
##         print coors.shape, nm.amax( econn )
##         pause()

        region = self.region
        group = region.domain.groups[self.ig]
        cells = region.get_cells( self.ig )
        mu.interp_vertex_data( coors, econn, mesh.coors,
                               group.conn[cells], bf, 0 )

    def get_v_data_shape(self, integral=None):
        """returns (n_el, n_qp, dim, n_ep)"""
        if integral is not None:
            bf_vg = self.get_base('v', 1, integral)
            return (self.region.shape[self.ig].n_cell,) + bf_vg.shape

        else:
            return (self.region.shape[self.ig].n_cell,
                    self.interp.gel.dim, self.n_ep['v'])

    def get_s_data_shape(self, integral, key):
        """returns (n_fa, n_qp, dim, n_fp)"""
        if not self.surface_data:
            return 0, 0, 0, 0

        sd = self.surface_data[key]
        bf_sg = self.get_base(sd.face_type, 1, integral)
        n_qp, dim, n_fp = bf_sg.shape
        assert_(n_fp == sd.n_fp)

        return sd.n_fa, n_qp, dim + 1, n_fp

    ##
    # c: 05.09.2006, r: 09.05.2008
    def setup_surface_data( self, region ):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        sd = FESurface('surface_data_%s' % region.name, region,
                       self.efaces, self.econn, self.ig)
        self.surface_data[region.name] = sd
        return sd

    ##
    # 11.07.2007, c
    def setup_point_data( self, field, region ):
        conn = field.get_dofs_in_region(region, merge=True, igs=region.igs)
##         conn = [nods]\
##                + [nm.empty( (0,), dtype = nm.int32 )]\
##                * (len( region.igs ) - 1)

        conn.shape += (1,)
        self.point_data[region.name] = conn

    def get_qp(self, key, integral):
        """
        Get quadrature points and weights corresponding to the given key
        and integral. The key is 'v' or 's#', where # is the number of
        face vertices.
        """
        qpkey = (integral.name, key)

        if not self.qp_coors.has_key(qpkey):
            interp = self.interp
            if (key[0] == 's') and not self.is_surface:
                dim = interp.gel.dim - 1
                n_fp = interp.gel.surface_facet.n_vertex
                geometry = '%d_%d' % (dim, n_fp)

            else:
                geometry = interp.gel.name

            vals, weights = integral.get_qp(geometry)
            self.qp_coors[qpkey] = Struct(vals=vals, weights=weights)

        return self.qp_coors[qpkey]

    def get_base(self, key, derivative, integral,
                 from_geometry=False, base_only=True):
        if self.is_surface:
            key = 'v'

        qp = self.get_qp(key, integral)

        if from_geometry:
            if key == 'v':
                gkey = key
            else:
                gkey = 's%d' % self.interp.gel.surface_facet.n_vertex
            ps = self.interp.gel.interp.poly_spaces[gkey]
            bf_key = (integral.name, 'g' + key, derivative)

        else:
            ps = self.interp.poly_spaces[key]
            bf_key = (integral.name, key, derivative)

        if not self.bf.has_key( bf_key ):
            self.bf[bf_key] = ps.eval_base(qp.vals, diff=derivative)

        if base_only:
            return self.bf[bf_key]
        else:
            return self.bf[bf_key], qp.weights

    def describe_geometry(self, field, gtype, region, integral=None, coors=None):
        """Compute jacobians, element volumes and base function derivatives
        for Volume-type geometries, and jacobians, normals and base function
        derivatives for Surface-type geometries.
        """
        if coors is None:
            coors = field.aps.coors

        if gtype == 'volume':
            if integral is None:
                from sfepy.fem import Integral
                dim = field.domain.shape.dim
                quad_name = 'gauss_o1_d%d' % dim
                integral = Integral('i_tmp', 'v', quad_name)

            qp = self.get_qp('v', integral)

            mapping = VolumeMapping(coors, self.econn,
                                    poly_space=self.interp.poly_spaces['v'])

            try:
                vg = mapping.get_mapping(qp.vals, qp.weights)

            except ValueError: # Constant bubble.
                domain = self.region.domain
                group = domain.groups[self.ig]
                coors = domain.get_mesh_coors()

                ps = self.interp.gel.interp.poly_spaces['v']
                mapping = VolumeMapping(coors, group.conn, poly_space=ps)

                vg = mapping.get_mapping(qp.vals, qp.weights)

            out = vg

        elif (gtype == 'surface') or (gtype == 'surface_extra'):
            sd = self.surface_data[region.name]
            qp = self.get_qp(sd.face_type, integral)

            if not self.is_surface:
                ps = self.interp.poly_spaces[sd.face_type]

            else:
                ps = self.interp.poly_spaces['v']

            econn = sd.get_connectivity(self.is_surface)
            mapping = SurfaceMapping(coors, econn, poly_space=ps)
            sg = mapping.get_mapping(qp.vals, qp.weights)

            if gtype == 'surface_extra':
                sg.alloc_extra_data( self.get_v_data_shape()[2] )

                self.create_bqp( region.name, integral )
                bf_bg = self.get_base(sd.bkey, 1, integral)
                sg.evaluate_bfbgm( bf_bg, coors, sd.fis, self.econn )

            out =  sg

        elif gtype == 'point':
            out = None

        else:
            raise ValueError('unknown geometry type: %s' % gtype)

        if out is not None:
            # Store the integral used.
            out.integral = integral

        return out

    def _create_bqp(self, skey, bf_s, weights, integral_name):
        interp = self.interp
        gel = interp.gel
        bkey = 'b%s' % skey[1:]
        bqpkey = (integral_name, bkey)
        coors, faces = gel.coors, gel.get_surface_entities()

        vals = _interp_to_faces(coors, bf_s, faces)
        self.qp_coors[bqpkey] = Struct(name = 'BQP_%s' % bkey,
                                       vals = vals, weights = weights)
        interp.poly_spaces[bkey] = interp.poly_spaces['v']
        return bkey

    def create_bqp(self, region_name, integral):
        sd = self.surface_data[region_name]
        bqpkey = (integral.name, sd.bkey)
        if not bqpkey in self.qp_coors:
            bf_s = self.get_base(sd.face_type, 0, integral,
                                 from_geometry=True)
            qp = self.get_qp(sd.face_type, integral)

            bkey = self._create_bqp(sd.face_type, bf_s, qp.weights,
                                    integral.name)
            assert_(bkey == sd.bkey)

##
# 14.07.2006, c
class Approximations( Container ):
    """- Region can span over several groups -> different Aproximation
    instances

    - interps and hence node_descs are per region (must have single
    geometry!)
    - no two interps can be in a same group -> no two aps (with different
    regions) can be in a same group -> aps can be uniquely indexed with ig"""

    def __init__(self, interp, region, is_surface=False):

        self.igs = region.igs
        self.interp = interp
        self.aps_per_group = {}
        self.region = region
        self.is_surface = is_surface

        objs = OneTypeList( Approximation )
        for ig in region.igs:
            ap = Approximation(interp.name + '_%s_ig%d' % (region.name, ig),
                               interp, region, ig, is_surface=is_surface)
            self.aps_per_group[ig] = ap
            objs.append(ap)

        self.update(objs)

    def iter_aps(self, igs=None):
        for ig, ap in self.aps_per_group.iteritems():
            if igs is not None:
                if not ig in igs: continue
            yield ig, ap

    def get_approximation(self, ig):
        """
        Returns
        -------
            Approximation instance for the given ig.
        """
        ap = self.aps_per_group[ig]

        return ap

    def setup_facet_orientations(self):
        self.node_desc = self.interp.describe_nodes()

        edge_nodes = self.node_desc.edge_nodes
        if edge_nodes is not None:
            ed = self.region.domain.ed
            self.edge_dof_perms = ed.get_facet_dof_permutations(edge_nodes)

        face_nodes = self.node_desc.face_nodes
        if face_nodes is not None:
            fa = self.region.domain.fa
            self.face_dof_perms = fa.get_facet_dof_permutations(face_nodes)

    def setup_global_base(self):
        """
        Setup global DOF/base function indices and connectivity.
        """
        self.init_econn()

        self.n_vertex_dof, remap = self.setup_vertex_dofs()
        self.n_edge_dof = self.setup_edge_dofs()
        self.n_face_dof = self.setup_face_dofs()
        self.n_bubble_dof = self.setup_bubble_dofs()

        self.n_nod = self.n_vertex_dof + self.n_edge_dof \
                     + self.n_face_dof + self.n_bubble_dof

        self.setup_esurface()

        return self.n_nod, remap

    def init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        for ig, ap in self.iter_aps():
            n_ep = ap.n_ep['v']
            n_cell = self.region.get_n_cells(ig, self.is_surface)
            ap.econn = nm.zeros((n_cell, n_ep), nm.int32)

    def setup_vertex_dofs(self):
        """
        Setup vertex DOF connectivity.
        """
        if self.node_desc.vertex is None:
            return 0

        region = self.region

        vertices = region.all_vertices
        n_dof = vertices.shape[0]
        remap = nm.empty((vertices.max() + 1,), dtype=nm.int32)
        remap.fill(-1)
        remap[vertices] = nm.arange(n_dof, dtype=nm.int32)

        ##
        # Remap vertex node connectivity to field-local numbering.
        for ig, ap in self.iter_aps():
            group = region.domain.groups[ig]
            if not self.is_surface:
                offset = group.shape.n_ep
                cells = region.get_cells(ig)
                ap.econn[:,:offset] = remap[group.conn[cells]]

            else:
                faces = group.gel.get_surface_entities()
                aux = FESurface('aux', region, faces, group.conn, ig)
                ap.econn[:,:aux.n_fp] = aux.leconn
                ap.surface_data[region.name] = aux

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
        for ig, ap in self.iter_aps():
            ii = get_facets(ig)
            uid_i = facets.uid_i[ii]

            uids.append(uid_i)

        uids = nm.unique(nm.concatenate(uids))
        n_uid = uids.shape[0]
        remap = nm.empty((uids.max() + 1,), dtype=nm.int32)
        remap.fill(-1)
        remap[uids] = nm.arange(n_uid, dtype=nm.int32)

        for ig, ap in self.iter_aps():
            ori = facets.oris[ig]
            perms = facet_perms[ig][ori]

            ii = get_facets(ig)
            g_uid = facets.uid_i[ii]
            uid = remap[g_uid]

            # Define global facet dof numbers.
            gdofs = nm.repeat(uid, n_dof_per_facet)
            gdofs.shape = (ii.shape[0], n_dof_per_facet)
            idof = nm.arange(n_dof_per_facet, dtype=nm.int32)
            gdofs = offset + n_dof_per_facet * gdofs + idof

            # Elements of facets.
            iel = facets.indices[ii, 1]

            ies = facets.indices[ii, 2]
            # DOF columns in econn for each facet.
            iep = facet_desc[ies]

            iaux = nm.arange(gdofs.shape[0], dtype=nm.int32)
            ap.econn[iel[:, None], iep] = gdofs[iaux[:, None], perms]

        n_dof = n_dof_per_facet * n_uid

        return n_dof

    def setup_edge_dofs(self):
        """
        Setup edge DOF connectivity.
        """
        if self.node_desc.edge is None:
            return 0

        return self._setup_facet_dofs(self.region.domain.ed,
                                      self.node_desc.edge,
                                      self.edge_dof_perms,
                                      self.region.get_edges,
                                      self.n_vertex_dof)

    def setup_face_dofs(self):
        """
        Setup face DOF connectivity.
        """
        if self.node_desc.face is None:
            return 0

        return self._setup_facet_dofs(self.region.domain.fa,
                                      self.node_desc.face,
                                      self.face_dof_perms,
                                      self.region.get_faces,
                                      self.n_vertex_dof + self.n_edge_dof)

    def setup_bubble_dofs(self):
        """
        Setup bubble DOF connectivity.
        """
        if self.node_desc.bubble is None:
            return 0

        offset = self.n_vertex_dof + self.n_edge_dof + self.n_face_dof
        n_dof = 0
        for ig, ap in self.iter_aps():
            n_bubble = self.node_desc.bubble.shape[0]
            n_cell = self.region.get_n_cells(ig, self.is_surface)
            aux = nm.arange(offset, offset + n_bubble * n_cell, dtype=nm.int32)
            aux.shape = (n_cell, n_bubble)
            iep = self.node_desc.bubble[0]
            ap.econn[:,iep:] = aux

            n_dof += n_bubble * n_cell

        return n_dof

    def setup_esurface(self):
        """
        Setup extended surface entities (edges in 2D, faces in 3D),
        i.e. indices of surface entities into the extended connectivity.
        """
        node_desc = self.node_desc

        for ig, ap in self.iter_aps():
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

    ##
    # c: 19.07.2006, r: 15.01.2008
    def setup_coors( self, mesh, cnt_vn ):
        """id column is set to 1."""
        noft = self.node_offset_table

        n_nod = noft[-1,-1]
        self.coors = nm.empty( (n_nod, mesh.dim), nm.float64 )
#        print n_nod

        # Mesh vertex nodes.
        inod = slice( noft[0,0], noft[0,-1] )
#        print inod
        if inod:
            indx = cnt_vn[cnt_vn >= 0]
            self.coors[inod,:] = mesh.coors[indx,:]
##         print self.coors
##         print mesh.coors
##         pause()
        for ig, ap in self.iter_aps():
            ap.eval_extra_coor( self.coors, mesh )

##         print self.coors
##         print self.coors.shape
##         pause()

    def setup_surface_data(self, region):
        for ig, ap in self.iter_aps(igs=region.igs):
            if region.name not in ap.surface_data:
                if self.is_surface:
                    msg = 'no surface data of surface field! (%s)' % region.name
                    raise ValueError(msg)

                ap.setup_surface_data(region)

    def setup_point_data(self, field, region):
        # Point data only in the first group to avoid multiple
        # assembling of nodes on group boundaries.
        for ig, ap in self.iter_aps(igs=region.igs[:1]):
            if region.name not in ap.point_data:
                ap.setup_point_data(field, region)

    def describe_geometry(self, field, geometry_type, ig, region,
                          term_region, integral):
        """
        For give approximation, compute jacobians, element volumes and
        base function derivatives for Volume-type geometries, and
        jacobians, normals and base function derivatives for
        Surface-type geometries.

        Usually, region is term_region. Only if is_trace is True, then region
        is the mirror region and term_region is the true term region.
        """
        is_trace = region is not term_region

        ap = self.aps_per_group[ig]

        geo = ap.describe_geometry(field, geometry_type, region, integral,
                                   self.coors)
        return geo

    ##
    # created:       21.12.2007
    # last revision: 21.12.2007
    def update_geometry( self, field, regions, geometries ):
        for geom_key, geom in geometries.iteritems():
            iname, geometry_type, tregion_name, ap_name = geom_key
            ap = self[ap_name]
            geom = ap.describe_geometry(field, geometry_type,
                                        regions[tregion_name],
                                        geom.integral, self.coors)
            geometries[geom_key] = geom
            
    ##
    # 03.05.2007, c
    # 17.07.2007
    def purge_surface_data( self ):
        for ap in self:
            ap.surface_data = {}
