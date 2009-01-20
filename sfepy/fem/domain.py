from sfepy.base.base import *
from sfepy.base.reader import Reader
from geomElement import GeomElement
from region import Region, sort_by_dependency
import fea
import extmods.meshutils as mu

##
# 14.01.2005
# 22.02.2005
# 04.08.2005
def dm_create_list( groups, sentinel, mode, is_sort ):

    n_gr = len( groups )

    if mode == 0:
        n_obj = [group.shape.n_edge_total for group in groups.itervalues()]
        nn = 2
    else:
        n_obj = [group.shape.n_face_total for group in groups.itervalues()]
        nn = 4
        
    objs = nm.zeros( (sum( n_obj ) + 2, 3 + nn), nm.int32 ) - 1
    gptr = nm.zeros( (n_gr + 1,), nm.int32 )
    eptr = nm.zeros( (n_gr,), nm.int32 )

    ii = 0
    for ig, group in groups.iteritems():
        conn, gel = group.conn, group.gel
        gd = gel.data['v']
        if (mode == 0):
            n_item, items = gd.n_edge, gd.edges
        else:
            n_item, items = gd.n_face, gd.faces

#        print '::', ii
        ii = mu.create_list( ii, objs, ig, conn, items, is_sort )[1]
#        print '::', ii

        eptr[ig] = n_item
        gptr[ig+1] = ii

    aux = nm.repeat( nm.array( [sentinel], nm.int32 ), nn )
    objs[-2,3:] = aux
    objs[-1,3:] = aux + 1

##         print gptr, eptr

    if (ii != sum( n_obj )):
        msg = 'neighbour_list size mismatch! (%d == %d = sum( %s ))'\
              % (ii, sum( n_obj ), n_obj)
        raise ValueError( msg )

    return( (objs, gptr, eptr) )

##
# c: 05.01.2005, r: 04.02.2008
def dm_neighbour_list( obj_in, groups, ic, perm, mode ):

    n_gr = len( groups )
    n_el = [group.shape.n_el for group in groups.itervalues()]
    if mode == 0:
        n_obj = [group.shape.n_edge_total for group in groups.itervalues()]
    else:
        n_obj = [group.shape.n_face_total for group in groups.itervalues()]

    data_s = obj_in.data_s

    pg = nm.zeros( (n_gr + 1,), dtype = nm.int32 )
    pg[1:] = n_el
    pg = nm.cumsum( pg, dtype = nm.int32 )

    pel = nm.zeros( (sum( n_el ) + 1,), nm.int32 )
    for ig in range( n_gr ):
        gd = groups[ig].gel.data['v']
        if (mode == 0):
            n_item = gd.n_edge
        else:
            n_item = gd.n_face
        for ie in range( n_el[ig] ):
            pel[pg[ig]+ie+1] = n_item
    pel = nm.cumsum( pel, dtype = nm.int32 )

    pobj = nm.zeros( (sum( n_obj ) + 1,), dtype = nm.int32 ) + 1
    mu.neighbour_list_ptr( pobj, pg, pel, data_s, ic, mode )
    cnt_pobj = pobj.copy()
    pobj = nm.cumsum( pobj, dtype = nm.int32 )

    objs = nm.zeros( (pobj[-1],), dtype = nm.int32 )
    uid = nm.zeros( (sum( n_obj ),), dtype = nm.int32 )
    cnt = nm.zeros( (sum( n_obj ),), dtype = nm.int32 )
    iu = mu.neighbour_list( objs, uid, cnt, \
                           pg, pel, pobj, data_s, obj_in.uid, ic, perm, mode )[1]


    if (nm.sometrue( cnt_pobj[1:] - cnt )):
        msg = '%s\n%s\n%s\n' % (cnt_pobj[1:], cnt, cnt_pobj[1:] - cnt)
        raise ValueError( msg )

    if (iu != obj_in.n_unique):
        msg = ' '.join( [iu, "==", obj_in.n_unique] )
        raise ValueError( msg )

    return( (pg, pel, pobj, objs, cnt, uid) )

##
# 05.01.2005
# 04.08.2005
# 31.10.2005
def dm_print_neighbour_list( objs, obj_list, pauses = False ):
    pg = obj_list.pg
    pel = obj_list.pel
    pd = obj_list.pd
    data = obj_list.data
    cnt = obj_list.cnt
    for ig in range( len( pg ) - 1 ):
        print ig
        for iel in range( pg[ig+1] - pg[ig] ):
            print "  ", iel
            n_edge = pel[pg[ig]+iel+1] - pel[pg[ig]+iel]
            for ii in range( n_edge ):
                cptr = pel[pg[ig]+iel] + ii
                ptr = pd[cptr]
                print "    ", ii, ":", obj_list.uid[cptr]
                for jj in range( cnt[cptr] ):
                    print "      ", jj, data[ptr+jj], \
                          objs.data[data[ptr+jj]]

            if pauses:
                spause()

##
# 21.12.2005, c
def dm_mark_surface_faces( fa, nfa ):
    """ flag: 0 .. inner, 1 .. triangle, 2 .. quadrangle"""
    flag = nm.zeros( (fa.data.shape[0],), nm.int32 )

    pd = nfa.pd
    data = nfa.data
    cnt = nfa.cnt
    for ii in xrange( len( cnt ) ):
        if (cnt[ii] == 1):
            ptr = pd[ii]
            if (fa.data[data[ptr],-1] == -1):
                flag[data[ptr]] = 1
            else:
                flag[data[ptr]] = 2
    return flag

##
# 17.07.2006, c
class Domain( Struct ):
    """
    Domain is divided into groups, whose purpose is to have homogeneous
    data shapes."""

    ##
    # 17.07.2006, c
    def from_mesh( mesh, component_dir ):

        read = Reader( component_dir )

        geom_els = OneTypeList( GeomElement )
        for desc in mesh.descs:
#            print desc
            
            if not geom_els.find( desc ):
                geom_els.append( read( GeomElement, desc ) )

        interps = {}
        for gel in geom_els:
            name = gel.interpolation

            if interps.has_key( name ):
                gel.interp = interps[name]
            else:
                gel.interp = read( fea.Interpolant, name )
                interps[name] = gel.interp

        obj = Domain( name = mesh.name,
                      mesh = mesh,
                      geom_els = geom_els,
                      geom_interps = interps )
        obj.mat_ids_to_i_gs = {}
        for ig, mat_id in enumerate( mesh.mat_ids ):
            obj.mat_ids_to_i_gs[mat_id[0]] = ig
#        obj.n_nod = obj.mesh.nod0.shape[0]
        
        return obj
    from_mesh = staticmethod( from_mesh )

    def setup_groups( self ):

        for gel in self.geom_els:
            gel.setup()
        
        # Call after gel.setup()
        for ginterp in self.geom_interps.itervalues():
            gel = self.geom_els[ginterp.geometry]
            ginterp.setup( gel )

        n_gr = len( self.mesh.conns )
        n_nod, dim = self.mesh.nod0.shape
        self.shape = Struct( n_gr = len( self.mesh.conns ), n_el = 0,
                             n_nod = n_nod, dim = dim )
        self.groups = {}
        for ii in range( self.shape.n_gr ):
            gel = self.geom_els[self.mesh.descs[ii]] # Shortcut.
            conn = self.mesh.conns[ii]
            vertices = nm.unique1d( conn )

            n_vertex = vertices.shape[0]
            n_el, n_ep = conn.shape
            n_edge = gel.data['v'].n_edge           
            n_edge_total = n_edge * n_el

            if gel.dim == 3:
                n_face = gel.data['v'].n_face
                n_face_total = n_face * n_el
            else:
                n_face = n_face_total = 0

            shape = Struct( n_vertex = n_vertex, n_el = n_el, n_ep = n_ep,
                            n_edge = n_edge, n_edge_total = n_edge_total,
                            n_face = n_face, n_face_total = n_face_total )
            self.groups[ii] = Struct( ig = ii,
                                      vertices = vertices,
                                      conn = conn,
                                      gel = gel,
                                      shape = shape )
            self.shape.n_el += n_el

    ##
    # c: 22.11.2007, r: 25.03.2008
    def iter_groups( self, igs = None ):
        if igs is None:
            for ig in xrange( self.shape.n_gr ): # sorted by ig.
                yield self.groups[ig]
        else:
            for ig in igs:
                yield ig, self.groups[ig]

    ##
    # c: 25.03.2008, r: 25.03.2008
    def get_cell_offsets( self ):
        offs = {}
        off = 0
        for group in self.iter_groups():
            ig = group.ig
            offs[ig] = off
            off += group.shape.n_el
        return offs

    ##
    # 22.08.2006, c
    def get_mesh_coors( self ):
        return self.mesh.nod0[:,:-1]

    ##
    # 30.08.2007, c
    def get_conns( self ):
        return self.mesh.conns

    
    ##
    # 08.06.2006, c
    # 17.07.2006
    def fix_element_orientation( self ):

        coors = self.mesh.nod0
        for ii, group in self.groups.iteritems():

            ori, conn = group.gel.orientation, group.conn

            itry = 0
            while itry < 2:
                flag = -nm.ones( conn.shape[0], dtype = nm.int32 )

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                mu.orient_elements( flag, conn, coors,
                                   ori.v_roots, ori.v_vecs,
                                   ori.swap_from, ori.swap_to )
    #            print flag
                if nm.alltrue( flag == 0 ):
                    if itry > 0: output( '...corrected' )
                    itry = -1
                    break

                output( 'warning: bad element orienation, trying to correct...' )
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError, "elements cannot be oriented! (%d, %s)"\
                      % (ii, self.mesh.descs[ii] )
            elif flag[0] == -1:
                output( 'warning: element orienation not checked' )

    ##
    # 19.07.2006, c
    def get_orientation( self, ii, mode = 'edges' ):
        group = self.groups[ii]
        gd = group.gel.data['v']

        if mode == 'edges':
            ori = nm.zeros( (group.shape.n_el, gd.n_edge), nm.int32 )
            mu.orient_edges( ori, group.conn, gd.edges );
        elif mode == 'faces':
            output( 'orient faces' )
            raise NotImplementedError
        else:
            output( mode )
            raise ValueError
        
        return ori

    def has_faces( self ):
        return sum( [group.shape.n_face
                     for group in self.iter_groups()] ) > 0
        

    ##
    # c: 04.08.2005, r: 20.02.2008
    def setup_neighbour_lists( self, create_edge_list = 1, create_face_list = 1 ):
        mode_names = ['edges', 'faces']
        
        is_face = self.has_faces()
        flags = [create_edge_list, create_face_list and is_face]
        sort_cols = [[3,4], [3,4,5,6]]

        for mode in range( 2 ):
            if flags[mode]:
                output( 'setting up domain %s...' % mode_names[mode] )

                tt = time.clock()
                obj = Struct()
                obj.data, obj.gptr, obj.eptr \
                          = dm_create_list( self.groups,
                                           self.mesh.nod0.shape[0], mode, 1 )

#                print "t = ", time.clock() - tt
                ii = nm.arange( obj.data.shape[0], dtype = nm.int32 );
                ii.shape = (ii.shape[0], 1)
                aux = nm.concatenate( (obj.data, ii), 1 ).copy()
    ##             print aux.flags['contiguous']
                mu.sort_rows( aux, nm.array( sort_cols[mode], nm.int32 ) )
    ##             print "->", aux

#                print "t = ", time.clock() - tt
                obj.perm = perm = aux[:,-1].copy()
                obj.data_s = aux[:,:-1].copy()
                aux = nm.arange( perm.shape[0], dtype = nm.int32 )
                obj.perm_i = nm.zeros_like( obj.perm )
                obj.perm_i[perm] = aux

##                 print perm
##                 print obj.perm_i

                ic = nm.where( nm.sum( nm.absolute( \
                    obj.data_s[1:,3:] - obj.data_s[:-1,3:] ), 1 ), 0, 1 )
		ic = nm.asarray( ic, dtype = nm.int32 )
##                 print ic, len( ic )
                obj.n_data = len( ic ) - 1
                obj.n_unique = len( ic ) - nm.sum( ic ) - 1
                obj.unique_list = nm.asarray( nm.where( ic[:-1] == 0 )[0],
					     dtype = nm.int32 )
#               print "t = ", time.clock() - tt

                assert_( len( obj.unique_list ) == obj.n_unique )
#                print obj.n_unique, obj.unique_list, obj.unique_list.shape
                ii = nm.cumsum( ic[:-1] == 0, dtype = nm.int32 )
                obj.uid = ii.copy()
                obj.uid[0], obj.uid[1:] = 0, ii[:-1]
                obj.uid_i = obj.uid[obj.perm_i[:-2]]
##                 print obj
##                 debug()

                nobj = Struct()
                nobj.pg, nobj.pel, nobj.pd, nobj.data, nobj.cnt, nobj.uid \
                         = dm_neighbour_list( obj, self.groups, ic, perm, mode )
#                print "t = ", time.clock() - tt


                if mode == 0:
#                    dm_print_neighbour_list( obj, nobj, pauses = True )
                    self.ed, self.ned = obj, nobj
                else:
                    self.fa, self.nfa = obj, nobj

                output( '...done in %.2f s' % (time.clock() - tt) )
        if not is_face:
            self.fa, self.nfa = None, None

    ##
    # 19.07.2006
    # 24.08.2006
    def get_neighbour_lists( self, force_faces = False ):
        if force_faces and not self.fa:
            return self.ed, self.ned, self.ed, self.ned
        else:
            return self.ed, self.ned, self.fa, self.nfa

    ##
    # c: 31.10.2005, r: 20.02.2008
    def create_regions( self, region_defs, funmod = None ):
        from sfepy.fem.parseReg import create_bnf, visit_stack, print_stack,\
             ParseException

        output( 'creating regions...' )
        tt = time.clock()
        regions = OneTypeList( Region )

        ##
        # 14.06.2006, c
        # 15.06.2006
        # 19.02.2007
        # 02.03.2007
        # 02.05.2007
        # 30.05.2007
        # 05.06.2007
        def region_leaf( domain, rdef, funmod ):
            def _region_leaf( level, op ):

                token, details = op['token'], op['orig']
                if token != 'KW_Region':
                    parse_def = token + '<' + ' '.join( details ) + '>'
##                     conns = [group.conn for group in domain.groups.itervalues()]
##                     vertex_groups = [group.vertices
##                                     for group in domain.groups.itervalues()]
                    region = Region( 'leaf', rdef, domain, parse_def )

                if token == 'KW_Region':
                    details = details[1][2:]
                    aux = regions.find( details )
                    if not aux:
                        raise ValueError, 'region %s does not exist' % details
                    else:
                        if rdef[:4] == 'copy':
                            region = aux.copy()
                        else:
                            region = aux

                elif token == 'KW_All':
                    region.set_vertices( nm.arange( domain.mesh.nod0.shape[0],
                                                   dtype = nm.int32 ) )
                elif token == 'E_NIR':
                    where = details[2]
                    
                    if where[0] == '[':
                        out = nm.array( eval( where ), dtype = nm.int32 )
                        assert_( nm.amin( out ) >= 0 )
                        assert_( nm.amax( out ) < domain.mesh.nod0.shape[0] )
                    else:
                        x = domain.mesh.nod0[:,0]
                        y = domain.mesh.nod0[:,1]
                        if domain.mesh.dim == 3:
                            z = domain.mesh.nod0[:,2]
                        else:
                            z = None
                        coor_dict = {'x' : x, 'y' : y, 'z': z}
                        
                        out = nm.where( eval( where, {}, coor_dict ) )[0]
                    region.set_vertices( out )
                    
                elif token == 'E_NOS':

                    if domain.fa: # 3D.
                        fa, nfa = domain.fa, domain.nfa
                    else:
                        fa, nfa = domain.ed, domain.ned
                        
                    flag = dm_mark_surface_faces( fa, nfa )
                    ii = nm.where( flag > 0 )[0]
                    aux = nm.unique1d( fa.data[ii,3:].ravel() )
                    if aux[0] == -1: # Triangular faces have -1 as 4. point.
                        aux = aux[1:]
                    region.can_cells = False
                    region.set_vertices( aux )

                elif token == 'E_NBF':
                    where = details[2]
                    
                    x = domain.mesh.nod0[:,0]
                    if domain.shape.dim > 1:
                        y = domain.mesh.nod0[:,1]
                        if domain.shape.dim > 2:
                            z = domain.mesh.nod0[:,2]
                        else:
                            z = None
                    else:
                        y = None
                    aux = {'x' : x, 'y' : y, 'z': z}
                        
                    fun = 'funmod.' + where
#                    print fun
                    out = nm.where( eval( fun, {'funmod' : funmod}, aux ) )[0]

                    region.set_vertices( out )

                elif token == 'E_EBF':
                    where = details[2]
                    
                    aux = {'domain' : domain}
                        
                    fun = 'funmod.' + where
#                    print fun
                    out = eval( fun, {'funmod' : funmod}, aux )
                    print out
                    region.set_cells( out )

                elif token == 'E_EOG':

                    group = int( details[3] )

                    ig = domain.mat_ids_to_i_gs[group]
                    group = domain.groups[ig]
                    region.set_from_group( ig, group.vertices, group.shape.n_el )

                elif token == 'E_ONIR':
                    aux = regions[details[3][2:]]
                    region.set_vertices( aux.all_vertices[0:1] )

                elif token == 'E_NI':
                    region.set_vertices( nm.array( [int( details[1] )],
                                                  dtype = nm.int32 ) )

                else:
                    output( 'token "%s" unkown - check regions!' % token )
                    raise NotImplementedError
                return region
            
            return _region_leaf

        ##
        # 14.06.2006, c
        # 15.06.2006
        def region_op( domain, rdef, funmod ):
            def _region_op( level, op, item1, item2 ):

                token = op['token']
                if token == 'OA_SubN':
                    return item1.sub_n( item2 )
                elif token == 'OA_SubE':
                    return item1.sub_e( item2 )
                elif token == 'OA_AddN':
                    return item1.add_n( item2 )
                elif token == 'OA_AddE':
                    return item1.add_e( item2 )
                elif token == 'OA_IntersectN':
                    return item1.intersect_n( item2 )
                elif token == 'OA_IntersectE':
                    return item1.intersect_e( item2 )
                else:
                    raise NotImplementedError, token
            return _region_op

        stack = []
        bnf = create_bnf( stack )

        ##
        # Sort region definitions by dependencies.
        depends = re.compile( 'r\.([a-zA-Z_0-9]+)' ).search
        graph = {}
        name_to_sort_name = {}
        for sort_name, rdef in region_defs.iteritems():
            name, sel = rdef.name, rdef.select
#            print sort_name, name, sel
            if name_to_sort_name.has_key( name ):
                raise 'region %s/%s already defined!' % (sort_name, name)
            name_to_sort_name[name] = sort_name

            if not graph.has_key( name ):
                graph[name] = [0]

            while len( sel ):
                aux = depends( sel )
                if aux:
                    graph[name].append( aux.group()[2:] )
                    sel = sel[aux.end():]
                else:
                    sel = ''
#        print graph

        sorted_regions = sort_by_dependency( graph )
#        print sorted_regions
        
        ##
        # Define regions.
        for name in sorted_regions:
            sort_name = name_to_sort_name[name]
            rdef = region_defs[sort_name]

            stack[:] = []
            try:
                out = bnf.parseString( rdef.select )
            except ParseException:
                print 'parsing failed:', rdef
                raise

#            print_stack( copy( stack ) )

            region = visit_stack( stack, region_op( self, rdef.select, funmod ),
                                 region_leaf( self, rdef.select, funmod ) )
            if hasattr( rdef, 'forbid' ):
                fb = re.compile( '^group +\d+(\s+\d+)*$' ).match( rdef.forbid )
                if fb:
                    groups = rdef.forbid[5:].strip().split()
                    forbid = [int( ii ) for ii in groups]
                else:
                    raise SyntaxError, 'bad forbid: %s' % rdef.forbid
                forbidden_igs = [self.mat_ids_to_i_gs[mat_id] for mat_id in forbid]
##                 print forbidden_igs
##                 pause()
                region.delete_groups( forbidden_igs )
            if hasattr( rdef, 'can_cells' ):
                region.switch_cells( rdef.can_cells )
            region.complete_description( self.ed, self.fa )

            region.type_name = region.name
            region.name = rdef.name
            region.sort_name = sort_name
            
            output( ' ', region.type_name, region.name, region.sort_name )
#            print region.definition
#            print region.parse_def
            regions.append( region )

        # Sort by definition name.
        regions.sort( cmp = lambda i1, i2: cmp( i1.sort_name, i2.sort_name ) )
        self.regions = regions
        output( '...done in %.2f s' % (time.clock() - tt) )

        return regions

    ##
    # 26.07.2007, c
    def get_diameter( self ):
        bbox = self.mesh.get_bounding_box()
        return (bbox[1,:] - bbox[0,:]).max()

    ##
    # 31.07.2007, c
    def get_element_diameters( self, ig, cells, vg, mode, square = True ):
        group = self.groups[ig]
        diameters = nm.empty( (len( cells ), 1, 1, 1), dtype = nm.float64 )
        if vg is None:
            diameters.fill( 1.0 )
        else:
            vg.get_element_diameters( diameters, group.gel.data['v'].edges,
                                    self.get_mesh_coors().copy(), group.conn,
                                    cells, mode )
        if square:
            out = diameters.squeeze()
        else:
            out = nm.sqrt( diameters.squeeze() )

        return out
            
    ##
    # 29.08.2007, re-c from 00.01.18, r: 26.03.2008
    def surface_faces( self ):

        if not self.fa:
            print "no faces defined!"
            raise ValueError

        fa = Struct()
        fa.data, fa.gptr, fa.eptr \
                 = dm_create_list( self.groups, self.mesh.nod0.shape[0], 1, 0 )

        flag = dm_mark_surface_faces( fa, self.nfa )

        surf_faces = []
        itri = nm.where( flag == 1 )[0]
        if itri.size:
            surf_faces.append( fa.data[itri,3:6] )
        itet = nm.where( flag == 2 )[0]
        if itet.size:
            surf_faces.append( fa.data[itet,3:7] )

        isurf = nm.where( flag >= 1 )[0]
        if isurf.size:
            lst = fa.data[isurf,0:3]

        return lst, surf_faces
