from sfepy.base.base import *
from feGenerators import ap_node_generators, ap_bf_generators
from baseFunction import BaseFunction
from quadratures import collect_quadratures
import extmods.meshutils as mu
import extmods.geometry as gm

def set_mesh_coors( domain, fields, geometries, coors, update_state = False ):
    domain.mesh.coors = coors.copy()
    if update_state:
        fields.setup_coors()
        for field in fields:
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

##
# 02.08.2005, c
# 30.09.2005
# 03.10.2005
# 07.10.2005
# 19.12.2005
# 02.08.2006
# 11.04.2007
# 03.05.2007
# 04.05.2007
# 24.05.2007
def eval_bf( qp, base_fun, nodes, derivative ):
    """bf_b(g) indexing is (ifa,iqp,:,n_ep) -> can use FMF_SetCell"""
    fun, var_set = base_fun.fun, base_fun.var_set
    if (qp.vals.ndim == 2):
        if derivative == 0:
            val = fun.value( qp.vals, nodes )
        else:
            val = fun.value( qp.vals, nodes, var_set )

    else: # Boundary QP.
        sh = qp.vals.shape[0:2]
        if derivative == 0:
            val = nm.zeros( sh + (1, len( nodes )), nm.float64 )
            for ifa, bqp in enumerate( qp.vals ):
                val[ifa,:,:,:] = fun.value( bqp, nodes )
        else:
            val = nm.zeros( sh + (len( var_set ), len( nodes )), nm.float64 )
            for ifa, bqp in enumerate( qp.vals ):
                val[ifa,:,:,:] = fun.value( bqp, nodes, var_set )


    if derivative == 0:
        aux = nm.sum( val, -1 ).ravel()
        if not nm.alltrue( nm.absolute( aux - 1.0 ) < 1e-14 ):
            raise ValueError( '%s' % aux )

    return val

##
# Field approximation.
# 29.11.2004, c
# 30.11.2004
# 01.12.2004
# 04.02.2005
class Interpolant( Struct ):

    ##
    # 31.03.2005
    # 16.06.2005
    def __init__( self ):
        self.setup_done = 0

    ##
    # 09.03.2005, c
    # 31.03.2005
    # 22.06.2005
    # 19.07.2005
    # 20.07.2005
    # 02.08.2005
    # 18.07.2006
    # 23.08.2006
    # 11.04.2007
    def setup( self, gel = None ):
        if (self.setup_done): return

        if gel is not None:
            self.gel = gel

        # Input transformation and defaults.
        for key, nod in self.nodes.iteritems():
            if not nod.has_key( 'share' ):
                nod['share'] = 1

        self.desc = dict_to_struct( self.desc, flag = (1,) )
        self.nodes, self.base_funs\
                    = dict_to_struct( self.nodes, self.base_funs,
                                    flag = (0, 1) )
        name = self.desc.family
        self.i_key_map = i_key_map = invert_dict( self.key_map )
        # Nodes.
        for key, nod in self.nodes.iteritems():
            if nod.mode == 'generate':
                try:
                    gen = ap_node_generators[name]
                except:
                    print self.desc.family
                    raise NotImplementedError, name
                if not key == 'v':
                    gen = gen[1]
                else:
                    gen = gen[0]
                aux = gen( self.desc.approx_order, i_key_map[key],
                           nod, self.gel, want_bar_coors = 1 )
                # ntx := node types, a row = (type, entity (vefb)) index
                nod.nts, nod.vals, nod.bar_coors = aux
            else:
                print 'unknown mode', nod.mode
                raise NotImplementedError

        # Base functions.
        for key, bf in self.base_funs.items():
            gd = self.gel.data[i_key_map[key]]
            bf.var_set = range( gd.dim )
            
            if bf.mode == 'generate':
                try:
                    gen = ap_bf_generators[name]
                except:
                    print self.desc.family
                    raise NotImplementedError, name
                if not key == 'v':
                    gen = gen[1]
                else:
                    gen = gen[0]
                bfgen = gen( gd.coors, bf, self.nodes[key].vals, bf.var_set )
                bf.fun = BaseFunction( bfgen, self.nodes[key].vals, bf.var_set )

            else:
                print 'unknown mode', bf.mode
                raise NotImplementedError


##         for key, nod in self.nodes.iteritems():
##             print key, nod
##         pause()
            
        self.setup_done = 1



    ##
    # 02.08.2005, c
    # 22.09.2005
    # 23.09.2005
    # 30.09.2005
    # 03.10.2005
    def list_extra_node_types( self, et, ft ):
        gd = self.gel.data['v']
        nodes = self.nodes['v']
        max_ao = nm.amax( nm.sum( nodes.vals, 1 ) )

##         print nodes
##         print 'sdsd', max_ao
##         pause()

        for ii, nt in enumerate( nodes.nts ):
            if (nt[0] == 1): # Edge node.
                edge = gd.edges[nt[1]]
                tpar = float( nodes.vals[ii,edge[1]] )
                key = int( round( tpar / max_ao, 5 ) * 1e5 )
                if not et.has_key( key ): et[key] = len( et )
##                 print ii, nt
##                 print tpar
##                 print nodes.vals[ii], edge

            elif (nt[0] == 2): # Face node.
                face = gd.faces[nt[1]]
                upar = float( nodes.vals[ii,face[1]] )
                vpar = float( nodes.vals[ii,face[-1]] )
                key = [int( round( ii / max_ao, 5 ) * 1e5 )
                       for ii in [upar, vpar]]
                key = tuple( key )
                if not ft.has_key( key ): ft[key] = len( ft )
##                 print upar, vpar
##                 print nodes.vals[ii], face

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

        nts = self.nodes['v'].nts

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
        for key, nodes in self.nodes.iteritems():
            nn[key] = nodes.vals.shape[0]
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
        nodes = self.interp.nodes['v']

        iex = (nodes.nts[:,0] > 0).nonzero()[0]
        qp_coors = nodes.bar_coors[iex,:]

        qp = Struct( vals = qp_coors )
        bf = eval_bf( qp, ginterp.base_funs['v'], ginterp.nodes['v'].vals, 0 )
        bf = bf[:,0,:].copy()

        ##
        # Evaulate extra coordinates with 'bf'.
#        print self.econn.shape, self.sub.n_el

        econn = nm.zeros_like( self.econn[:,iex] ).copy()
        for ii in range( 1, 4 ):
            ix = nm.where( nodes.nts[:,0] == ii )[0]
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
            return (self.region.shape[self.ig].n_cell,)\
                   + self.bf[(iname, 'v', 1)].shape
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
        bf_sg = self.bf[(iname, sd.face_type, 1)]
        n_qp, dim, n_fp = bf_sg.shape
        assert_( n_fp == sd.n_fp )
        
        return sd.n_fa, n_qp, dim + 1, n_fp

    ##
    # c: 05.09.2006, r: 09.05.2008
    def setup_surface_data( self, region ):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        face_indices = region.fis[self.ig]
        
        faces = self.efaces[face_indices[:,1]]
        if faces.size == 0:
            raise RuntimeError, 'region with group with no faces! (%s)'\
                  % region.name
#        print self.efaces
        try:
            ee = self.econn[face_indices[:,0]]
        except:
            pdb.set_trace()

        econn = nm.empty( faces.shape, dtype = nm.int32 )
        for ir, face in enumerate( faces ):
#            print ir, face, ee[ir,face]
            econn[ir] = ee[ir,face]

        ef = econn.flat
        nodes = nm.unique1d( ef )

        aux = -nm.ones( (nm.max( ef ) + 1,), dtype = nm.int32 )
        aux[nodes] = nm.arange( len( nodes ), dtype = nm.int32 )
        leconn = aux[econn].copy()
        assert_( nm.alltrue( nodes[leconn] == econn ) )
        
        n_fa, n_fp = face_indices.shape[0], faces.shape[1]
        face_type = 's%d' % n_fp

        sd = Struct( name = 'surface_data_%s' % region.name,
                     econn = econn, fis = face_indices, n_fa = n_fa, n_fp = n_fp,
                     nodes = nodes, leconn = leconn, face_type = face_type )
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
                if dim == 2:
                    n_fp = len( interp.gel.s_faces[interp.i_key_map[key]] )
                elif dim == 1:
                    n_fp = len( interp.gel.s_edges[interp.i_key_map[key]][0] )
                else:
                    raise NotImplementedError
                geometry = '%d_%d' % (dim, n_fp)
            else:
                geometry = interp.geometry
            vals, weights = integral.get_qp( geometry )()
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
            gkey = self.interp.i_key_map[key]
            base_fun = self.interp.gel.interp.base_funs[gkey]
            nodes = self.interp.gel.interp.nodes[gkey].vals
            bf_key = (iname, 'g' + key, derivative)
        else:
            base_fun = self.interp.base_funs[key]
            nodes = self.interp.nodes[key].vals
            bf_key = (iname, key, derivative)

        if not self.bf.has_key( bf_key ):
            self.bf[bf_key] = eval_bf( qp, base_fun, nodes, derivative )            

        if base_only:
            return self.bf[bf_key]
        else:
            return self.bf[bf_key], qp.weights

    ##
    # 05.09.2006, c
    # 03.05.2007
    def describe_geometry( self, field, geom_request, integral, coors ):
        gtype = geom_request.gtype
        if gtype == 'Volume':
            bf_vg, weights = self.get_base( 'v', 1, integral = integral,
                                          base_only = False )

##             print bf_vg, weights
##             pause()
            vg = gm.VolumeGeometry( *self.get_v_data_shape( integral.name ) )
            vg.mode = gm.GM_Material
            try:
                vg.describe( coors, self.econn, bf_vg, weights )
            except:
#                pdb.set_trace()
                gm.errclear()
                raise
            return vg

        elif (gtype == 'Surface') or (gtype == 'SurfaceExtra'):
##             print self.name
##             print geom_request
            region = geom_request.region
            fa = region.domain.get_neighbour_lists( True )[2]
            region.setup_face_indices( fa )
            region.select_cells_of_surface()
            self.setup_surface_data( region )
            
            sd = self.surface_data[region.name]
##             print sd
##             print integral
            bf_sg, weights = self.get_base( sd.face_type, 1, integral = integral,
                                          base_only = False )
##             print bf_sg, weights 

            sg = gm.SurfaceGeometry( *self.get_s_data_shape( integral.name,
                                                          region.name ) )
            sg.mode = gm.GM_Material

##             print sg
            
            try:
                sg.describe( coors, sd.econn, bf_sg, weights )
            except:
                gm.errclear()
                raise

            if gtype == 'SurfaceExtra':
                sg.alloc_extra_data( self.get_v_data_shape()[2] )

                bf_s = self.get_base( sd.face_type, 0, integral = integral,
                                    from_geometry = True )
##                 print bf_s
                bkey = self.create_bqp( sd.face_type, bf_s, weights,
                                        integral.name )
                # Store bkey in SurfaceData, so that base function can be
                # queried later.
                sd.bkey = bkey

                bf_bg = self.get_base( bkey, 1, integral = integral )
##                 print bf_bg
                sg.evaluate_bfbgm( bf_bg, coors, sd.fis, self.econn )

##                 sg.str( sys.stdout, 0 )
##                 pause()
            
            return sg

        elif gtype == 'Point':
            region = geom_request.region
##             print gtype, region.name, self.ig
##             pause()
            if self.ig == region.igs[0]:
                # Point data only in the first group to avoid multiple
                # assembling of nodes on group boundaries.
                self.setup_point_data( field, region )

        else:
            print 'unknown geometry type:', gtype
            raise ValueError

    ##
    #
    def create_bqp( self, skey, bf_s, weights, iname ):
        interp = self.interp
        gd = interp.gel.data['v']
        bkey = 'b%s' % skey[1:]
        bqpkey = (iname, bkey)
        coors, faces = gd.coors, gd.faces

        vals = _interp_to_faces( coors, bf_s, faces )
        self.qp_coors[bqpkey] = Struct( name = 'BQP_%s' % bkey,
                                     vals = vals, weights = weights )
##         print self.qp_coors[bqpkey]
        interp.nodes[bkey] = interp.nodes['v']
        interp.base_funs[bkey] = interp.base_funs['v']
        return bkey

##
# 14.07.2006, c
class Approximations( Container ):
    """- Region can span over several groups -> different Aproximation
    instances

    - interps and hence node_descs are per region (must have single
    geometry!)
    - no two interps can be in a same group -> no two aps (with different
    regions) can be in a same group -> aps can be uniquely indexed with ig"""

    ##
    # c: 23.11.2007, r: 15.01.2008
    def __init__( self, bases, interps, domain ):

        self.igs = []
        self.interps = {}
        self.aps_per_group = {}
        self.region_names_per_group = {}
        objs = OneTypeList( Approximation )
        for region_name, base_name in bases.iteritems():
##             print region_name, base_name

            try:
                region = domain.regions[region_name]
            except:
                output( 'region %s does not exist' % region_name )
                raise

            interp = interps[base_name]
            if self.interps.has_key( region.name ):
                if self.interps[region.name] is not interp:
                    output( 'interpolation mismatch!' )
                    output( self.interps[region.name].name, interp.name )
                    raise ValueError
            else:
                interp.setup( domain.geom_els[interp.geometry] )
                self.interps[region.name] = interp
            
            for ig in region.igs:
                if ig in self.igs:
                    output( 'base regions must not share groups! (%s, %d)'\
                            % (region.name, ig) )
                    raise ValueError

                self.igs.append( ig )
                
                ap = Approximation( base_name + '_%s_ig%d' % (region.name, ig),
                                    interp, region, ig )
                self.aps_per_group[ig] = ap
                self.region_names_per_group[ig] = region.name
                objs.append( ap )

        self.update( objs )

        self.clear_geometries()

    ##
    # created:       21.12.2007
    # last revision: 21.12.2007
    def clear_geometries( self ):
        self.geometries = {}

    ##
    # c: 23.11.2007, r: 15.01.2008
    def iter_aps( self, igs = None ):
        for ig, ap in self.aps_per_group.iteritems():
            if igs is not None:
                if not ig in igs: continue
            yield self.region_names_per_group[ig], ig, ap

    ##
    # c: 19.07.2006, r: 15.01.2008
    def describe_nodes( self ):
        self.ent = ent = {}
        self.fnt = fnt = {}

        self.node_descs = {}
        for region_name, interp in self.interps.iteritems():
            interp.list_extra_node_types( ent, fnt )
            node_desc = interp.describe_nodes()
            self.node_descs[region_name] = node_desc

##             print ent, fnt
##         pause()

        for region_name, ig, ap in self.iter_aps():
            ap.describe_nodes( self.node_descs[region_name] )
            

    ##
    # c: 23.09.2005, r: 15.01.2008
    def setup_nodes( self ):

        self.edge_oris = {}
        self.face_oris = {}
        for key1, key2, ap in self.iter_aps():
##             print ap
##             print key1, key2
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
        
    ##
    # c: 30.09.2005, r: 04.02.2008
    def setup_global_base( self, domain ):
        """
        efaces: indices of faces into econn.
        """

        node_offset_table = nm.zeros( (4, len( self ) + 1), dtype = nm.int32 )
##         print node_offset_table.shape
        
        i_vertex, i_edge, i_face, i_bubble = 0, 1, 2, 3

        # Global node number.
        iseq = 0

        ##
        # Vertex nodes.
        n_v = domain.shape.n_nod
        cnt_vn = nm.empty( (n_v,), dtype = nm.int32 )
        cnt_vn.fill( -1 )
        
        node_offset_table[i_vertex,0] = iseq
        ia = 0
        for region_name, ig, ap in self.iter_aps():
            region = ap.region
            node_desc = self.node_descs[region_name]
            n_ep = ap.n_ep['v']
            group = region.domain.groups[ig]
            
            ap.econn = nm.zeros( (region.shape[ig].n_cell, n_ep), nm.int32 )
##             print ap.econn.shape
#            pause()
            if node_desc.vertex.size:
                offset = group.shape.n_ep
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
        for region_name, ig, ap in self.iter_aps():
            region = ap.region
            node_desc = self.node_descs[region_name]
            group = region.domain.groups[ig]
            if node_desc.vertex.size:
                offset = group.shape.n_ep
                cells = region.get_cells( ig )
                ap.econn[:,:offset] = remap[group.conn[cells]]
##                 print group.conn, nm.amax( group.conn )
##                 print ap.econn, nm.amax( ap.econn )
##                 pause()

        ed, ned, fa, nfa = domain.get_neighbour_lists()
        entt = self.ent_table
        cnt_en = nm.zeros( (entt.shape[1], ed.n_unique), nm.int32 ) - 1

        ##
        # Edge nodes.
        node_offset_table[i_edge,0] = iseq
        ia = 0
        for region_name, ig, ap in self.iter_aps():
            region = ap.region
            node_desc = self.node_descs[region_name]
            group = region.domain.groups[ig]
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
        for region_name, ig, ap in self.iter_aps():
            node_offset_table[i_face,ia+1] = iseq
            ia += 1

        #    cnt_bn = nm.zeros( (fntt.shape[1], fa.n_unique), nm.int32 ) - 1
        ##
        # Bubble nodes.
        node_offset_table[i_bubble,0] = iseq
        ia = 0
        for region_name, ig, ap in self.iter_aps():
            region = ap.region
            node_desc = self.node_descs[region_name]
            group = region.domain.groups[ig]
            if (node_desc.bubble):
                n_bubble = node_desc.bubble[0].shape[0]
                n_el = region.shape[ig].n_cell
                aux = nm.arange( iseq, iseq + n_bubble * n_el )
                aux.shape = (n_el, n_bubble)
                offset = node_desc.bubble[0][0,0]
                ap.econn[:,offset:] = aux[:,:]
                iseq += n_bubble * n_el

            node_offset_table[i_bubble,ia+1] = iseq
            ia += 1
            
##         print node_offset_table
        if node_offset_table[-1,-1] != iseq:
            raise RuntimeError
        
        self.node_offset_table = node_offset_table

        ia = 0
        for region_name, ig, ap in self.iter_aps():
            ap.node_offsets = self.node_offset_table[:,ia:ia+2]
            ia += 1
##             print ia
##             print ap.econn
##             print ap.node_offsets
#            pause()
            
        for region_name, ig, ap in self.iter_aps():
            node_desc = self.node_descs[region_name]

            gd = ap.interp.gel.data['v']
            ap.efaces = gd.faces.copy()

            if ap.has_extra_edge_nodes:
                nd = node_desc.edge
                efs = []
                for eof in gd.edges_of_faces:
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
                
##             print ap.efaces
##             print ap.interp.base_funs['v'].fun
##             try:
##                 print ap.interp.base_funs['s2'].fun
##                 print ap.interp.base_funs['s3'].fun
##             except:
##                 pass
##             pause()

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
        for region_name, ig, ap in self.iter_aps():
            ap.eval_extra_coor( self.coors, mesh )

##         print self.coors
##         print self.coors.shape
##         pause()

    ##
    # c: 12.10.2005, r: 15.01.2008
    def describe_geometry( self, field, geometries, geom_request, integral,
                          over_write = False ):

        for region_name, ig, ap in self.iter_aps( igs = geom_request.region.igs ):
##             print region_name, ig, ap.name

            # Store integral for possible future base function request.
            ap.integrals[integral.name] = integral
            ap.dim = field.dim

            ##
            # Prepare common bases.
            gtype = geom_request.gtype
            if gtype == 'Volume':
                ap.get_base( 'v', 0, integral = integral )
                ap.get_base( 'v', 1, integral = integral )
            elif (gtype == 'Surface') or (gtype == 'SurfaceExtra'):
                pass

            geom_key = (integral.name, geom_request.gtype,
                        geom_request.region.name, ap.name)
##            print field.name, geom_key

            if geometries.has_key( geom_key ):
                self.geometries[geom_key] = geometries[geom_key]
            else:
                if self.geometries.has_key( geom_key ):
                    geometries[geom_key] = self.geometries[geom_key]
                else:
##                print 'new geometry: %s of %s' % (geom_key, ap.name)
                    geom = ap.describe_geometry( field, geom_request, integral,
                                                 self.coors )
                    self.geometries[geom_key] = geometries[geom_key] = geom
                    # Make an alias Surface -> SurfaceExtra.
                    if geom_request.gtype == 'SurfaceExtra':
                        key2 = list(geom_key)
                        key2[1] = 'Surface'
                        key2 = tuple(key2)
                        self.geometries[key2] = geometries[key2] = geom

        if over_write:
            self.geometries = geometries

    ##
    # created:       21.12.2007
    # last revision: 21.12.2007
    def update_geometry( self, field, regions, geometries ):
        for geom_key, geom in self.geometries.iteritems():
            iname, gtype, tregion_name, ap_name = geom_key
            ap = self[ap_name]
            integral = ap.integrals[iname]
            geom_request = Struct( gtype = gtype, region = regions[tregion_name] )
            geom = ap.describe_geometry( field, geom_request, integral,
                                         self.coors )
            self.geometries[geom_key] = geometries[geom_key] = geom
            
    ##
    # 03.05.2007, c
    # 17.07.2007
    def purge_surface_data( self ):
        for ap in self:
            ap.surface_data = {}
