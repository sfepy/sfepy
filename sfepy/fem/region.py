from sfepy.base.base import *

_depends = re.compile( 'r\.([a-zA-Z_0-9]+)' ).findall

def get_parents(selector):
    """Given a region selector, return names of regions it is based on."""
    parents = _depends(selector)

    return parents

def get_dependency_graph(region_defs):
    """Return a dependency graph and a name-sort name mapping for given
    region definitions."""
    graph = {}
    name_to_sort_name = {}
    for sort_name, rdef in region_defs.iteritems():
        name, sel = rdef.name, rdef.select
##         print sort_name, name, sel
        if name_to_sort_name.has_key( name ):
            msg = 'region %s/%s already defined!' % (sort_name, name)
            raise ValueError(msg)
        name_to_sort_name[name] = sort_name

        if not graph.has_key( name ):
            graph[name] = [0]

        for parent in get_parents(sel):
            graph[name].append(parent)
##     print graph
    return graph, name_to_sort_name

##
# 15.06.2006, c
# 17.07.2006
# 04.09.2006
def sort_by_dependency( graph ):

    out = []

    n_nod = len( graph )
    idone = 0
    idone0 = -1
    while idone < n_nod:

        dep_removed = 0
        for node, deps in graph.iteritems():
#            print '--', node, deps
            if (len( deps ) == 1) and not deps[0]:
                out.append( node )
                deps[0] = 1
                idone += 1
            elif not deps[0]:
#                print '--->', deps
                for ii, dep in enumerate( deps[1:] ):
                    if graph[dep][0]:
                        ir = deps.index( dep )
                        deps.pop( ir )
                        dep_removed += 1
#                print '---<', deps

##         print graph
##         print out
##         print idone, idone0, n_nod, dep_removed
##         pause()

        if (idone <= idone0) and not dep_removed:
            raise ValueError, 'circular dependency'
        idone0 = idone

    return out

##
# 15.06.2006, c
def _join( def1, op, def2 ):
    return '(' + def1 + ' ' + op + ' ' + def2 + ')'

##
# 31.10.2005, c
class Region( Struct ):

    ##
    # 14.06.2006, c
    # 15.06.2006
    # 23.02.2007
    def __init__( self, name, definition, domain, parse_def ):
        """conns, vertex_groups are links to domain data"""
        Struct.__init__(self,
                        name = name, definition = definition,
                        n_v_max = domain.shape.n_nod, domain = domain,
                        parse_def = parse_def, all_vertices = None,
                        igs = [], vertices = {}, edges = {}, faces = {},
                        cells = {}, fis = {},
                        volume = {}, surface = {}, length = {},
                        can_cells = True, must_update = True,
                        is_complete = False,
                        mirror_region = None, ig_map = None,
                        ig_map_i = None)

    ##
    # 15.06.2006, c
    def light_copy( self, name, parse_def ):
        return Region( name, self.definition, self.domain, parse_def )

    ##
    # c: 15.06.2006, r: 04.02.2008
    def update_groups( self, force = False ):
        """Vertices common to several groups are listed only in all of them -
        fa, ed.unique_indx contain no edge/face duplicates already."""
        if self.must_update or force:

            self.igs = []
            self.vertices = {}
            self.cells = {}

            for group in self.domain.iter_groups():
                ig = group.ig
                vv = nm.intersect1d( group.vertices, self.all_vertices )
                if len( vv ) == 0: continue

                self.igs.append( ig )
                self.vertices[ig] = vv

                if self.can_cells:
                    mask = nm.zeros( self.n_v_max, nm.int32 )
                    mask[vv] = 1

                    conn = group.conn
                    aux = nm.sum( mask[conn], 1, dtype = nm.int32 )
                    rcells = nm.where( aux == conn.shape[1] )[0]
                    self.cells[ig] = nm.asarray( rcells, dtype = nm.int32 )
        self.must_update = False

    ##
    # 15.06.2006, c
    def update_vertices( self ):
        self.all_vertices = nm.zeros( (0,), nm.int32 )
        self.vertices = {}
        for ig, group in self.domain.iter_groups( self.igs ):
            rcells = self.cells[ig]
            conn = group.conn
            nods = conn[rcells,:].ravel()
            aux = nm.unique1d( nods )
            self.vertices[ig] = aux
            self.all_vertices = nm.unique1d( nm.r_[self.all_vertices, aux] )
        
    ##
    # 15.06.2006, c
    def set_vertices( self, vertices ):

        self.all_vertices = nm.array(vertices, dtype=nm.int32)
        self.update_groups( force = True )
        self.is_complete = False

    ##
    # c: 15.06.2006, r: 14.07.2008
    def set_cells( self, cells ):

        self.igs = []
        self.cells = {}
        for ig, rcells in cells.iteritems():
            self.cells[ig] = nm.array( rcells, dtype = nm.int32, ndmin = 1 )
            self.igs.append( ig )
        self.update_vertices()
        self.is_complete = False
        self.must_update = False

    ##
    # 15.06.2006, c
    def set_from_group( self, ig, vertices, n_cell ):

        self.igs = [ig]
        self.cells = {ig : nm.arange( n_cell, dtype = nm.int32 )}
        self.vertices = {ig: vertices.copy()}
        self.all_vertices = vertices.copy()
        self.must_update = False

    ##
    # c: 23.02.2007, r: 22.01.2008
    def delete_groups( self, digs ):
        """self.complete_description must be called after!"""
        for ig in digs:
            try:
                del self.vertices[ig]
                del self.cells[ig]
                self.igs.remove( ig )
            except KeyError:
                pass

    ##
    # 17.07.2007, c
    def switch_cells( self, can_cells ):
        if self.can_cells:
            self.can_cells = can_cells
            if not can_cells:
                self.cells = {}
        else:
            self.can_cells = can_cells
            if can_cells:
                self.update_groups( force = True )
        
    ##
    # 15.06.2006, c
    # 02.08.2006
    # 21.02.2007
    # 23.02.2007
    # 30.03.2007
    def complete_description( self, ed, fa ):
        """self.edges, self.faces list edge/face indices per group (pointers to
        ed.data, fa.data) - repetitions among groups are possible."""
        ##
        # Get edges, faces, etc. par subdomain.
        edges = ed.data

        if fa is not None:
            faces = fa.data.copy()
            # Treat both 3, 4 node faces.
            ii = nm.where( faces[:,-1] == -1 )[0]
            faces[ii,-1] = faces[ii,3]

        self.edges = {}
        self.faces = {}
        self.shape = {}
        for ig, group in self.domain.iter_groups( self.igs ):
            vv = self.vertices[ig]
            if self.cells.has_key( ig ):
                n_cell = self.cells[ig].shape[0]
            else:
                n_cell = 0
            self.shape[ig] = Struct( n_vertex = vv.shape[0],
                                     n_cell = n_cell )
            if len( vv ) == 0: continue

            mask = nm.zeros( self.n_v_max, nm.int32 )
            mask[vv] = 1

            ied = nm.arange( ed.gptr[ig], ed.gptr[ig+1], dtype = nm.int32 )
            aux = nm.sum( mask[edges[ied,3:5]], 1 )
            # Points to ed.data.
            redges = ied[nm.where( aux == 2 )[0]]
            self.edges[ig] = redges
            if fa is None: continue
            
            ifa = nm.arange( fa.gptr[ig], fa.gptr[ig+1], dtype = nm.int32 )
            aux = nm.sum( mask[faces[ifa,3:7]], 1 )
            # Points to fa.data.
            rfaces = ifa[nm.where( aux == 4 )[0]]
            self.faces[ig] = rfaces

            self.shape[ig].n_edge = redges.shape[0]
            self.shape[ig].n_face = rfaces.shape[0]
            
        self.is_complete = True

    def setup_face_indices(self, reset=True):
        """
        Initialize an array (per group) of (iel, ifa) for each face.
        """
        if reset or not self.fis:
            fa = self.domain.get_neighbour_lists(force_faces=True)[2]

            if self.faces:
                faces = self.faces
            else:
                faces = self.edges

            self.fis = {}
            for ig in self.igs:
                rfaces = faces[ig]
                aux = fa.data[rfaces]
                assert_( nm.all( aux[:,0] == ig ) )
                fi = aux[:,1:3].copy()
                self.fis[ig] = fi

    ##
    # 05.09.2006, c
    # 22.02.2007
    # 17.07.2007
    def select_cells( self, n_verts ):
        """Select cells containing at least n_verts[ii] vertices per group ii."""
        if not self.can_cells:
            print 'region %s cannot have cells!' % self.name
            raise ValueError
        
        self.cells = {}
        for ig, group in self.domain.iter_groups( self.igs ):
            vv = self.vertices[ig]
            if len( vv ) == 0: continue
            
            mask = nm.zeros( self.n_v_max, nm.int32 )
            mask[vv] = 1

            aux = nm.sum( mask[group.conn], 1 )
            rcells = nm.where( aux >= n_verts[ig] )[0]
#            print rcells.shape
            self.cells[ig] = rcells

    def select_cells_of_surface(self, reset=True):
        """
        Select cells corresponding to faces (or edges in 2D).
        """
        if not self.can_cells:
            raise ValueError('region %s cannot have cells!' % self.name)

        self.setup_face_indices(reset=reset)

        if self.faces:
            faces = self.faces
        else:
            faces = self.edges

        self.cells = {}
        for ig in self.igs:
            rcells = self.fis[ig][:,0]
            self.cells[ig]= rcells

    ##
    # 02.03.2007, c
    def copy( self ):
        """Vertices-based copy."""
        tmp = self.light_copy( 'copy', self.parse_def )
        tmp.set_vertices( copy( self.all_vertices ) )
        return tmp
        
    ##
    # 15.06.2006, c
    def sub_n( self, other ):
        tmp = self.light_copy( 'op',
                              _join( self.parse_def, '-n', other.parse_def ) )
        tmp.set_vertices( nm.setdiff1d( self.all_vertices,
                                       other.all_vertices ) )
        
        return tmp

    ##
    # 15.06.2006, c
    def add_n( self, other ):
        tmp = self.light_copy( 'op',
                              _join( self.parse_def, '+n', other.parse_def ) )
        tmp.set_vertices( nm.union1d( self.all_vertices,
                                     other.all_vertices ) )
        
        return tmp

    ##
    # 15.06.2006, c
    def intersect_n( self, other ):
        tmp = self.light_copy( 'op',
                              _join( self.parse_def, '*n', other.parse_def ) )
        tmp.set_vertices( nm.intersect1d( self.all_vertices,
                                         other.all_vertices ) )
        
        return tmp

    ##
    # c: 15.06.2006, r: 15.04.2008
    def sub_e( self, other ):
        tmp = self.light_copy( 'op',
                              _join( self.parse_def, '-e', other.parse_def ) )
        for ig in self.igs:
            if ig not in other.igs:
                tmp.igs.append( ig )
                tmp.cells[ig] = self.cells[ig].copy()
                continue
            
            aux = nm.setdiff1d( self.cells[ig], other.cells[ig] )
            if not len( aux ): continue
            tmp.cells[ig] = aux
            tmp.igs.append( ig )

        tmp.update_vertices()
        return tmp

    ##
    # 15.06.2006, c
    def add_e( self, other ):
        tmp = self.light_copy( 'op',
                              _join( self.parse_def, '+e', other.parse_def ) )
        for ig in self.igs:
            tmp.igs.append( ig )
            if ig not in other.igs:
                tmp.cells[ig] = self.cells[ig].copy()
                continue

            tmp.cells[ig] = nm.union1d( self.cells[ig],
                                        other.cells[ig] )

        for ig in other.igs:
            if ig in tmp.igs: continue
            tmp.igs.append( ig )
            tmp.cells[ig] = other.cells[ig].copy()

        tmp.update_vertices()
        return tmp

    ##
    # 15.06.2006, c
    # 20.02.2007
    def intersect_e( self, other ):
        tmp = self.light_copy( 'op',
                              _join( self.parse_def, '*e', other.parse_def ) )
        for ig in self.igs:
            if ig not in other.igs: continue
            aux = nm.intersect1d( self.cells[ig], other.cells[ig] )
            if not len( aux ): continue
            tmp.igs.append( ig )
            tmp.cells[ig] = aux

        tmp.update_vertices()
        return tmp

    def setup_mirror_region(self):
        """
        Find the corresponding mirror region, set up element mapping.
        """
        for reg in self.domain.regions:
            if (reg is not self) and \
                   (len(reg.igs) == len(self.igs)) and \
                   nm.all(self.all_vertices == reg.all_vertices):
                mirror_region = reg
                break
        else:
            raise ValueError('cannot find mirror region! (%s)' % self.name)

        ig_map = {}
        ig_map_i = {}
        for igr in self.igs:
            for igc in mirror_region.igs:
                if nm.all(self.vertices[igr] ==
                          mirror_region.vertices[igc]):
                    ig_map[igc] = igr
                    ig_map_i[igr] = igc
                    break
            else:
                raise ValueError('cannot find mirror region group! (%d)' \
                                 % igr)

        self.mirror_region = mirror_region
        self.ig_map = ig_map
        self.ig_map_i = ig_map_i

    def get_mirror_region(self):
        return self.mirror_region, self.ig_map, self.ig_map_i

    def get_field_nodes(self, field, merge=False, clean=False,
                        warn=False, igs=None):
        """
        Get nodes of the field contained in the region.

        Notes
        -----
        For one edge node type only! (should index row of cnt_en...)
        """
        if igs is None:
            igs = self.igs
        cnt_en = field.cnt_en

        nods = []
        node_descs = field.get_node_descs( self )
        for ig, node_desc in node_descs.iteritems():
            if not ig in igs:
                nods.append( None )
                continue
            
            nnew = nm.empty( (0,), dtype = nm.int32 )
            if node_desc.vertex.size:
                nnew = nm.concatenate( (nnew, field.remap[self.vertices[ig]]) )

            if node_desc.edge:
                ed = field.domain.ed
                # ed.uid_i[self.edges[ii]] == ed.uid[ed.perm_i[self.edges[ii]]]
                enods = cnt_en[:cnt_en.shape[0],ed.uid_i[self.edges[ig]]].ravel()
                enods = nm.compress( (enods >= 0), enods )
                nnew = nm.concatenate( (nnew, enods) )

            if node_desc.face:
                print self.name, field.name
                raise NotImplementedError

            if node_desc.bubble and self.can_cells:
                noft = field.aps.node_offset_table
                ia = field.aps.igs.index( ig )
                enods = self.cells[ig] + noft[3,ia]
                nnew = nm.concatenate( (nnew, enods) )

            nods.append( nnew )

        if merge:
            nods = [nn for nn in nods if nn is not None]
            nods = nm.unique1d( nm.hstack( nods ) )

        elif clean:
            for nn in nods[:]:
                if nn is None:
                    nods.remove(nn)
                    if warn is not None:
                        output(warn + ('%s' % region.name))
            
        return nods

    ##
    # 22.02.2007, c
    def get_vertices( self, ig ):
        return self.vertices[ig]

    ##
    # 05.06.2007, c
    def get_edges( self, ig ):
        return self.edges[ig]
        
    ##
    # 05.06.2007, c
    def get_faces( self, ig ):
        return self.faces[ig]
        
    ##
    # 05.06.2007, c
    def get_cells( self, ig ):
        return self.cells[ig]

    def iter_cells(self):
        ii = 0
        for ig, cells in self.cells.iteritems():
            for iel in cells:
                yield ig, ii, iel
                ii += 1
        
    ##
    # created:       28.05.2007
    # last revision: 11.12.2007
    def has_cells( self ):

        if self.can_cells:
            for cells in self.cells.itervalues():
                if cells.size:
                    return True
            return False
        else:
            return False

    def has_cells_if_can( self ):
        if self.can_cells:
            for cells in self.cells.itervalues():
                if cells.size:
                    return True
            return False
        else:
            return True

    def update_geometry_info( self, field, key, mode = 'volume' ):
        """
        key: iname, aregion_name, ig
        TODO: surfaces, lengths
        ?call for all regions & fields in describe_geometry()?"""
        if self.has_cells():
            val = self.volume.setdefault( field.name, {} )
        else:
            self.volume[field.name] = 0.0
            return

        iname, ig = key
        aps = field.aps
        geometries = aps.geometries
        ap = aps.aps_per_group[ig]
        g_key = (iname, mode, self.name, ap.name)

        vg = geometries[g_key]

        volume = vg.variable( 2 )
        volume.shape = (nm.prod( volume.shape ),)
        if mode == 'volume':
            val[key] = nm.sum( volume[self.cells[ig]] )
        else:
            val[key] = nm.sum( volume[:] )
        self.volume[field.name] = val

    def get_volume(self, field, key=None, update=False, mode='volume'):
        if update:
            self.update_geometry_info( field, key, mode )

        if key is None:
            return self.volume[field.name]
        else:
            return self.volume[field.name][key]

    def contains( self, other ):
        """Tests only igs for now!!!"""
        return set( other.igs ).issubset( set( self.igs ) )

    ##
    # c: 25.03.2008, r: 25.03.2008
    def get_cell_offsets( self ):
        offs = {}
        off = 0
        for ig in self.igs:
            offs[ig] = off
            off += self.shape[ig].n_cell
        return offs

    def get_charfun( self, by_cell = False, val_by_id = False ):
        """
        Return the characteristic function of the region as a vector of values
        defined either in the mesh nodes (by_cell == False) or cells. The
        values are either 1 (val_by_id == False) or sequential id + 1.
        """
        if by_cell:
            chf = nm.zeros( (self.domain.shape.n_el,), dtype = nm.float64 )
            offs = self.get_cell_offsets()
            for ig, cells in self.cells.iteritems():
                iel = offs[ig] + cells
                if val_by_id:
                    chf[iel] = iel + 1
                else:
                    chf[iel] = 1.0
        else:
            chf = nm.zeros( (self.domain.shape.n_nod,), dtype = nm.float64 )
            if val_by_id:
                chf[self.all_vertices] = self.all_vertices + 1
            else:
                chf[self.all_vertices] = 1.0

        return chf
