from sfepy.base.base import *
from sfepy.fem.mappings import VolumeMapping, SurfaceMapping
from poly_spaces import PolySpace
from fe_surface import FESurface
import extmods.meshutils as mu
import extmods.geometry as gm

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

def _get_i_name( iname, integral, key ):
    if iname is None:
        if integral is None:
            print 'no integral name given for key "%s"' % key
            raise ValueError
        else:
            iname = integral.name
    return iname

class Interpolant( Struct ):
    """A simple wrapper around PolySpace."""

    # A dirty hack to satisfy the tests. The key should be given by PolySpace.
    _map = {
        '2_3_P2'  : 's3',
        '2_3_P2B' : 's3',
        '3_4_P2'  : 's6',
        '3_4_P2B' : 's6',
    }

    def __init__(self, name, gel):
        self.name = name
        self.gel = gel

        is_bubble = name[-1] == 'B'
        if is_bubble:
            approx_order = int(name[-2])

        else:
            approx_order = int(name[-1])

        if name in self._map:
            skey = self._map[name]

        else:
            if gel.surface_facet is not None:
                n_fp = gel.surface_facet.n_vertex

            else:
                n_fp = 2
            skey = 's%d' % n_fp

        key_map = {0 : 'v', 1 : skey}

        is_bubble = name[-1] == 'B'
        ## self.i_key_map = i_key_map = invert_dict(self.key_map)

        poly_spaces = {}
        for ii, gel in enumerate([gel, gel.surface_facet]):
            if gel is None: continue
            force_bubble = is_bubble and (ii == 0)

            ps = PolySpace.any_from_args(None, gel, approx_order,
                                         base='lagrange',
                                         force_bubble=force_bubble)
            key = key_map[ii]
            poly_spaces[key] = ps

        self.poly_spaces = poly_spaces

    ##
    # 02.08.2005, c
    # 22.09.2005
    # 23.09.2005
    # 30.09.2005
    # 03.10.2005
    def list_extra_node_types( self, et, ft ):
        gel = self.gel
        ps = self.poly_spaces['v']
        max_ao = nm.amax( nm.sum( ps.nodes, 1 ) )

##         print nodes
##         print 'sdsd', max_ao
##         pause()

        for ii, nt in enumerate( ps.nts ):
            if (nt[0] == 1): # Edge node.
                edge = gel.edges[nt[1]]
                tpar = float( ps.nodes[ii,edge[1]] )
                key = int( round( tpar / max_ao, 5 ) * 1e5 )
                if not et.has_key( key ): et[key] = len( et )
##                 print ii, nt
##                 print tpar
##                 print ps.nodes[ii], edge

            elif (nt[0] == 2): # Face node.
                face = gel.faces[nt[1]]
                upar = float( ps.nodes[ii,face[1]] )
                vpar = float( ps.nodes[ii,face[-1]] )
                key = [int( round( ii / max_ao, 5 ) * 1e5 )
                       for ii in [upar, vpar]]
                key = tuple( key )
                if not ft.has_key( key ): ft[key] = len( ft )
##                 print upar, vpar
##                 print ps.nodes[ii], face

        return (et, ft)


    ##
    # 1. col ... iep, 2. col ... is_extra (simplified for now)
    # c: 03.10.2005, r: 04.02.2008
    def describe_nodes( self ):

        def _describe_nodes_aux( inds ):
            objs = []
            ii = 0
            while (ii < len( inds )):
                nt = nts[inds[ii]]
                first = nt[1]
                obj = []
                while (first == nt[1]):
                    obj.append( [inds[ii], 1] )
                    ii += 1
                    if (ii >= len( inds )):
                        break
                    nt = nts[inds[ii]]
                objs.append( nm.array( obj, dtype = nm.int32 ) )

            return objs

        nts = self.poly_spaces['v'].nts

        node_desc = Struct()

        # Vertex node.
        ii = (nts[:,0] == 0).nonzero()[0]
        zz = nm.zeros( (ii.shape[0], 1), dtype = nm.int32 )
        node_desc.vertex = nm.concatenate( (ii[:,nm.newaxis], zz), 1 )

        # Edge nodes.
        inds = (nts[:,0] == 1).nonzero()[0]
        node_desc.edge = _describe_nodes_aux( inds )
        
        # Face nodes.
        inds = (nts[:,0] == 2).nonzero()[0]
        node_desc.face = _describe_nodes_aux( inds )

        # Bubble nodes.
        inds = (nts[:,0] == 3).nonzero()[0]
        node_desc.bubble = _describe_nodes_aux( inds )
        
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
    def __init__( self, name, interp, region, ig ):
        """interp, region are borrowed."""

        self.name = name
        self.interp = interp
        self.region = region
        self.ig = ig
        self.surface_data = {}
        self.edge_data = {}
        self.point_data = {}
        self.qp_coors = {}
        self.bf = {}
        self.integrals = {}
        self.n_ep = self.interp.get_n_nodes()

    ##
    # 19.07.2006, c
    def describe_nodes( self, node_desc ):
        self.node_desc = node_desc

        self.has_extra_edge_nodes = self.has_extra_face_nodes = False
        if self.node_desc.edge:
            self.has_extra_edge_nodes = True
        if self.node_desc.face:
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

    ##
    # 24.07.2006, c
    def get_v_data_shape( self, iname = None ):
        """returns (n_el, n_qp, dim, n_ep)"""
        if iname:
            bf_vg = self.get_base('v', 1, iname)
            return (self.region.shape[self.ig].n_cell,) + bf_vg.shape

        else:
            return (self.region.shape[self.ig].n_cell,
                    self.interp.gel.dim, self.n_ep['v'])

    ##
    # 05.09.2006, c
    # 20.02.2007
    # 17.07.2007
    def get_s_data_shape( self, iname, key ):
        """returns (n_fa, n_qp, dim, n_fp)"""
        if not self.surface_data:
            return 0, 0, 0, 0
        sd = self.surface_data[key]
        bf_sg = self.get_base(sd.face_type, 1, iname)
        n_qp, dim, n_fp = bf_sg.shape
        assert_( n_fp == sd.n_fp )
        
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
        conn = region.get_field_nodes( field, merge = True, igs = region.igs )
##         conn = [nods]\
##                + [nm.empty( (0,), dtype = nm.int32 )]\
##                * (len( region.igs ) - 1)

        conn.shape += (1,)
        self.point_data[region.name] = conn

    ##
    # created:       29.11.2007
    # last revision: 02.10.2008
    def get_qp( self, key, iname, integral = None ):
        """integrals are stored in self.integrals..."""
        qpkey = (iname, key)

        if not self.qp_coors.has_key( qpkey ):
            if integral is None:
                try:
                    integral = self.integrals[iname]
                except:
                    print self.name
                    print 'no integral given for key "%s"' % str( qpkey )
                    raise ValueError

            interp = self.interp
            if key[0] == 's':
                dim = interp.gel.dim - 1
                n_fp = interp.gel.surface_facet.n_vertex
                geometry = '%d_%d' % (dim, n_fp)

            else:
                geometry = interp.gel.name

            vals, weights = integral.get_qp( geometry )
            self.qp_coors[qpkey] = Struct( vals = vals, weights = weights )
##             print self.name, self.qp_coors
##             pause()
        return self.qp_coors[qpkey]

    ##
    # created:       29.11.2007
    def get_base( self, key, derivative, iname = None, integral = None,
                 from_geometry = False, base_only = True ):
        iname = _get_i_name( iname, integral, key )
        qp = self.get_qp( key, iname, integral )

        if from_geometry:
            if key == 'v':
                gkey = key
            else:
                gkey = 's%d' % self.interp.gel.surface_facet.n_vertex
            ps = self.interp.gel.interp.poly_spaces[gkey]
            bf_key = (iname, 'g' + key, derivative)

        else:
            ps = self.interp.poly_spaces[key]
            bf_key = (iname, key, derivative)

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

            qp = self.get_qp('v', integral.name, integral)

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
            qp = self.get_qp(sd.face_type, integral.name, integral)

            ps = self.interp.poly_spaces[sd.face_type]
            mapping = SurfaceMapping(coors, sd.econn, poly_space=ps)

            sg = mapping.get_mapping(qp.vals, qp.weights)

            if gtype == 'surface_extra':
                sg.alloc_extra_data( self.get_v_data_shape()[2] )

                self.create_bqp( region.name, integral )
                bf_bg = self.get_base( sd.bkey, 1, integral = integral )
                sg.evaluate_bfbgm( bf_bg, coors, sd.fis, self.econn )

            out =  sg

        elif gtype == 'point':
            out = None

        else:
            raise ValueError('unknown geometry type: %s' % gtype)

        return out

    ##
    #
    def _create_bqp( self, skey, bf_s, weights, iname ):
        interp = self.interp
        gel = interp.gel
        bkey = 'b%s' % skey[1:]
        bqpkey = (iname, bkey)
        coors, faces = gel.coors, gel.get_surface_entities()

        vals = _interp_to_faces( coors, bf_s, faces )
        self.qp_coors[bqpkey] = Struct( name = 'BQP_%s' % bkey,
                                     vals = vals, weights = weights )
##         print self.qp_coors[bqpkey]
        interp.poly_spaces[bkey] = interp.poly_spaces['v']
        return bkey

    def create_bqp(self, region_name, integral ):
        sd = self.surface_data[region_name]
        bqpkey = (integral.name, sd.bkey)
        if not bqpkey in self.qp_coors:
            bf_s = self.get_base( sd.face_type, 0, integral = integral,
                                  from_geometry = True )
            qp = self.get_qp( sd.face_type, None, integral )

            bkey = self._create_bqp( sd.face_type, bf_s, qp.weights,
                                     integral.name )
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

    def __init__(self, interp, region):

        self.igs = region.igs
        self.interp = interp
        self.aps_per_group = {}
        self.region = region

        objs = OneTypeList( Approximation )
        for ig in region.igs:
            ap = Approximation(interp.name + '_%s_ig%d' % (region.name, ig),
                               interp, region, ig )
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

    ##
    # c: 19.07.2006, r: 15.01.2008
    def describe_nodes( self ):
        self.ent = ent = {}
        self.fnt = fnt = {}

        self.node_descs = {}
        self.interp.list_extra_node_types(ent, fnt)
        self.node_desc = self.interp.describe_nodes()

##             print ent, fnt
##         pause()

        for ig, ap in self.iter_aps():
            ap.describe_nodes(self.node_desc)

    ##
    # c: 23.09.2005, r: 15.01.2008
    def setup_nodes( self ):

        self.edge_oris = {}
        self.face_oris = {}
        for _, ap in self.iter_aps():
##             print ap
            domain = ap.region.domain
            igs = ap.region.igs
            for ig in igs:
                if not self.edge_oris.has_key( ig ):
                    self.edge_oris[ig] = domain.get_orientation( ig )
                if not self.face_oris.has_key( ig ):
                    self.face_oris[ig] = None
#            pause()

        ##
        # Edge node type permutation table.
        n_ent = len( self.ent )
        entt = nm.zeros( (2, n_ent), nm.int32 )
        entt[0] = nm.arange( n_ent )
        entt[1] = nm.arange( n_ent - 1, -1, -1 )
        self.ent_table = entt

##         print self.ent
##         print self.ent_table
#        pause()

        ##
        # Face node type permutation table.

    def setup_global_base(self, is_surface=False):
        """
        efaces: indices of faces into econn.
        """
        region = self.region
        node_desc = self.node_desc

        node_offset_table = nm.zeros( (4, len( self ) + 1), dtype = nm.int32 )

        i_vertex, i_edge, i_face, i_bubble = 0, 1, 2, 3

        # Global node number.
        iseq = 0

        ##
        # Vertex nodes.
        n_v = region.n_v_max
        cnt_vn = nm.empty( (n_v,), dtype = nm.int32 )
        cnt_vn.fill( -1 )
        
        node_offset_table[i_vertex,0] = iseq
        ia = 0
        for ig, ap in self.iter_aps():
            n_ep = ap.n_ep['v']
            n_cell = region.get_n_cells(ig, is_surface)
            ap.econn = nm.zeros((n_cell, n_ep), nm.int32)

##             print ap.econn.shape
#            pause()
            if node_desc.vertex.size:
                vertices = region.get_vertices( ig )
                n_new = (nm.where( cnt_vn[vertices] == -1 )[0]).shape[0]
                cnt_vn[vertices] = vertices
#                print n_new
                iseq += n_new
            node_offset_table[i_vertex,ia+1] = iseq
            ia += 1

        ##
        # Remap vertex node connectivity to field-local numbering.
        indx = nm.arange( iseq, dtype = nm.int32 )
        remap = nm.empty( (n_v,), dtype = nm.int32 )
        remap.fill( -1 )
        remap[nm.where( cnt_vn >= 0 )[0]] = indx
##         print remap
##         pause()
##         print iseq, remap
##         pause()
        for ig, ap in self.iter_aps():
            group = region.domain.groups[ig]
            if node_desc.vertex.size:
                if not is_surface:
                    offset = group.shape.n_ep
                    cells = region.get_cells( ig )
                    ap.econn[:,:offset] = remap[group.conn[cells]]

                else:
                    faces = group.gel.get_surface_entities()
                    aux = FESurface('aux', region, faces, group.conn, ig)
                    ap.econn[:,:aux.n_fp] = aux.leconn

        ed, ned, fa, nfa = region.domain.get_neighbour_lists()
        entt = self.ent_table
        cnt_en = nm.zeros( (entt.shape[1], ed.n_unique), nm.int32 ) - 1

        ##
        # Edge nodes.
        node_offset_table[i_edge,0] = iseq
        ia = 0
        for ig, ap in self.iter_aps():
            if node_desc.edge:
                cptr0 = int( ned.pel[ned.pg[ig]] )
                ori = self.edge_oris[ig]
                iseq = mu.assign_edge_nodes( iseq, ap.econn, cnt_en, \
                                           ori, entt, ned.uid, \
                                           node_desc.edge, cptr0 )[1]
##                 print ap.econn
##                 pause()
            node_offset_table[i_edge,ia+1] = iseq
            ia += 1

        #    cnt_fn = nm.zeros( (fntt.shape[1], fa.n_unique), nm.int32 ) - 1
        node_offset_table[i_face,0] = iseq
        ia = 0
        for ig, ap in self.iter_aps():
            node_offset_table[i_face,ia+1] = iseq
            ia += 1

        #    cnt_bn = nm.zeros( (fntt.shape[1], fa.n_unique), nm.int32 ) - 1
        ##
        # Bubble nodes.
        node_offset_table[i_bubble,0] = iseq
        ia = 0
        for ig, ap in self.iter_aps():
            if (node_desc.bubble):
                n_bubble = node_desc.bubble[0].shape[0]
                n_cell = region.get_n_cells(ig, is_surface)
                aux = nm.arange( iseq, iseq + n_bubble * n_cell )
                aux.shape = (n_cell, n_bubble)
                offset = node_desc.bubble[0][0,0]
                ap.econn[:,offset:] = aux[:,:]
                iseq += n_bubble * n_cell

            node_offset_table[i_bubble,ia+1] = iseq
            ia += 1
            
##         print node_offset_table
        if node_offset_table[-1,-1] != iseq:
            raise RuntimeError
        
        self.node_offset_table = node_offset_table

        ia = 0
        for ig, ap in self.iter_aps():
            ap.node_offsets = self.node_offset_table[:,ia:ia+2]
            ia += 1
##             print ia
##             print ap.econn
##             print ap.node_offsets
#            pause()
            
        for ig, ap in self.iter_aps():
            gel = ap.interp.gel
            ap.efaces = gel.get_surface_entities().copy()

            if ap.has_extra_edge_nodes:
                nd = node_desc.edge
                efs = []
                for eof in gel.get_edges_per_face():
                    ef = [nd[ie][:,0] for ie in eof]
                    efs.append( ef )
                efs = nm.array( efs ).squeeze()
                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
#                print efs
                ap.efaces = nm.hstack( (ap.efaces, efs ) )

            if ap.has_extra_face_nodes:
                nd = node_desc.face
                efs = [ef[:,0] for ef in nd]
#                print nd
                efs = nm.array( efs ).squeeze()
                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
#                print efs
                ap.efaces = nm.hstack( (ap.efaces, efs ) )
                
        return iseq, remap, cnt_vn, cnt_en

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
            integral = ap.integrals[iname]
            geom = ap.describe_geometry(field, geometry_type,
                                        regions[tregion_name],
                                        integral, self.coors)
            geometries[geom_key] = geom
            
    ##
    # 03.05.2007, c
    # 17.07.2007
    def purge_surface_data( self ):
        for ap in self:
            ap.surface_data = {}
