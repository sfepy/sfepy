import re
from copy import copy

import numpy as nm

from sfepy.base.base import assert_, Struct

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

def _try_delete(obj, ig):
    try:
        del obj[ig]
    except KeyError:
        pass

##
# 31.10.2005, c
class Region( Struct ):

    @staticmethod
    def from_vertices(vertices, domain, name='region',
                      igs=None, can_cells=False, surface_integral=False):
        """
        Create a new region containing given vertices.

        Parameters
        ----------
        vertices : array
            The array of vertices.
        domain : Domain instance
            The domain containing the vertices.
        name : str, optional
            The name of the region.
        igs : list, optional
            The allowed element groups. Other groups will be ignored,
            even though the region might have vertices in them - the
            same effect the 'forbid' flag has.
        can_cells : bool, optional
            If True, the region can have cells.
        surface_integral : bool, optional
            If True, then each region surface facet (edge in 2D, face in
            3D) can be listed only in one group.

        Returns
        -------
        obj : Region instance
            The new region.
        """
        obj = Region(name, 'given vertices', domain, '')

        obj.set_vertices(vertices)

        if igs is not None:
            forbidden = nm.setdiff1d(obj.igs, igs)
            obj.delete_groups(forbidden)

        obj.switch_cells(can_cells) 
        obj.complete_description(domain.ed, domain.fa,
                                 surface_integral=surface_integral)

        return obj

    @staticmethod
    def from_faces(faces, domain, name='region',
                   igs=None, can_cells=False):
        """
        Create a new region containing given faces.

        Parameters
        ----------
        faces : array
            The array with indices to `domain.fa`.
        domain : Domain instance
            The domain containing the faces.
        name : str, optional
            The name of the region.
        igs : list, optional
            The allowed element groups. Other groups will be ignored,
            even though the region might have faces in them - the
            same effect the 'forbid' flag has.
        can_cells : bool, optional
            If True, the region can have cells.

        Returns
        -------
        obj : Region instance
            The new region.
        """
        obj = Region(name, 'given faces', domain, '')
        obj.set_faces(faces, igs=igs, can_cells=can_cells)

        return obj

    def __init__(self, name, definition, domain, parse_def):
        """
        Create region instance.

        Parameters
        ----------
        name : str
            The region name, either given, or automatic for intermediate
            regions.
        definition : str
            The region selector definition.
        domain : Domain instance
            The domain of the region.
        parse_def : str
            The parsed definition of the region.

        Notes
        -----
        conns, vertex_groups are links to domain data.
        """
        Struct.__init__(self,
                        name=name, definition=definition,
                        n_v_max=domain.shape.n_nod, domain=domain,
                        parse_def=parse_def, all_vertices=None,
                        igs=[], vertices={}, edges={}, faces={},
                        cells={}, fis={},
                        can_cells=True, true_cells={}, must_update=True,
                        is_complete=False,
                        mirror_region=None, ig_map=None,
                        ig_map_i=None)

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
                    self.true_cells[ig] = True

                else:
                    self.true_cells[ig] = False

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
            aux = nm.unique( nods )
            self.vertices[ig] = aux
            self.all_vertices = nm.unique( nm.r_[self.all_vertices, aux] )
        
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

    def set_from_group(self, ig, vertices, n_cell):
        """
        Set region to contain the given element group.
        """
        self.igs = [ig]
        self.cells = {ig : nm.arange( n_cell, dtype = nm.int32 )}
        self.true_cells[ig] = True
        self.vertices = {ig: vertices.copy()}
        self.all_vertices = vertices.copy()
        self.must_update = False

    def set_faces(self, faces, igs=None, can_cells=False):
        """
        Set region data using given faces. The region description is
        complete afterwards.

        Parameters
        ----------
        faces : array
            The array with indices to `domain.fa`.
        igs : list, optional
            The allowed element groups. Other groups will be ignored,
            even though the region might have faces in them.
        can_cells : bool, optional
            If True, the region can have cells.
        """
        faces = nm.asarray(faces)

        fa = self.domain.fa
        ed = self.domain.ed

        indices = fa.indices[faces]
        facets = fa.facets[faces]

        faces_igs = indices[:, 0]

        self.igs = nm.unique(faces_igs)
        if igs is not None:
            self.igs = nm.intersect1d(self.igs, igs)

        all_vertices = []
        self.vertices = {}
        self.edges = {}
        self.faces = {}
        self.cells = {}

        mask = nm.zeros(self.n_v_max, dtype=nm.bool)

        for ig, group in self.domain.iter_groups(self.igs):
            n_fp = fa.n_fps[ig]

            ii = faces_igs == ig

            vv = nm.unique(facets[ii, :n_fp])
            if len(vv) == 0: continue

            self.vertices[ig] = vv

            all_vertices.append(vv)

            self.faces[ig] = faces[ii]
            self.edges[ig] = ed.get_complete_facets(vv, ig, mask)

            if can_cells:
                mask.fill(False)
                mask[vv] = True

                conn = group.conn
                aux = nm.sum(mask[conn], 1, dtype=nm.int32)
                rcells = nm.where(aux == conn.shape[1])[0]
                self.cells[ig] = nm.asarray(rcells, dtype=nm.int32)
                self.true_cells[ig] = True

            else:
                self.true_cells[ig] = False

        self.all_vertices = nm.unique(nm.hstack(all_vertices))

        self.update_shape()
        self.is_complete = True
        self.must_update = False

    def delete_groups(self, digs):
        """
        Delete given element groups from the region.
        """
        for ig in digs:
            _try_delete(self.vertices, ig)
            _try_delete(self.cells, ig)
            _try_delete(self.faces, ig)
            _try_delete(self.edges, ig)
            try:
                self.igs.remove(ig)
            except ValueError:
                pass

        self.update_shape()

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
        
    def complete_description(self, ed, fa, surface_integral=False):
        """
        Complete the region description by listing edges and faces for
        each element group.

        Parameters
        ----------
        ed : Facets instance
            The edge facets.
        fa : Facets instance
            The face facets.
        surface_integral : bool
            If True, the each region surface facet (edge in 2D, face in
            3D) can be listed only in one group. Sub-entities are
            updated accordingly (vertices in 2D, vertices and edges in
            3D).

        Notes
        ------
        If `surface_integral` is False, `self.edges`, `self.faces` simply
        list edge/face indices per group (pointers to `ed.facets`,
        `fa.facets`) - repetitions among groups are possible.
        """
        ##
        # Get edges, faces, etc. par subdomain.
        mask = nm.zeros(self.n_v_max, dtype=nm.bool)

        self.edges = {}
        self.faces = {}

        if surface_integral:
            if self.domain.shape.dim == 2:
                allowed = nm.ones(ed.n_unique, dtype=nm.bool)
                facets = ed
                surf = self.edges

            else:
                allowed = nm.ones(fa.n_unique, dtype=nm.bool)
                facets = fa
                surf = self.faces

            # Get unique surface facets.
            empty_igs = []
            for ig, group in self.domain.iter_groups(self.igs):
                vv = self.vertices[ig]
                if len(vv) == 0: continue

                mask.fill(False)
                mask[vv] = True

                ifacets = facets.get_complete_facets(vv, ig, mask)
                ii = facets.uid_i[ifacets]
                surf[ig] = ifacets[allowed[ii]]
                allowed[ii] = False

                if not len(surf[ig]):
                    empty_igs.append(ig)

            self.delete_groups(empty_igs)

            # Update vertices, cells, and, in 3D, edges.
            for ig, group in self.domain.iter_groups(self.igs):
                n_fp = facets.n_fps[ig]
                vv = nm.unique(facets.facets[surf[ig], :n_fp])

                self.vertices[ig] = vv

                mask.fill(False)
                mask[vv] = True

                conn = group.conn
                aux = nm.sum(mask[conn], 1, dtype=nm.int32)
                rcells = nm.where(aux == conn.shape[1])[0]
                self.cells[ig] = nm.asarray(rcells, dtype=nm.int32)
                self.true_cells[ig] = True

                if self.domain.shape.dim == 3:
                    self.edges[ig] = ed.get_complete_facets(vv, ig, mask)

        else:
            for ig, group in self.domain.iter_groups(self.igs):
                vv = self.vertices[ig]
                if len(vv) == 0: continue

                mask.fill(False)
                mask[vv] = True

                # Points to ed.facets.
                self.edges[ig] = ed.get_complete_facets(vv, ig, mask)

                if fa is None: continue

                # Points to fa.facets.
                self.faces[ig] = fa.get_complete_facets(vv, ig, mask)

        self.update_shape()

        self.is_complete = True

    def update_shape(self):
        """
        Update shape of each group according to region vertices, edges,
        faces and cells.
        """
        aux = nm.array([])

        self.shape = {}
        for ig in self.igs:
            n_vertex = self.vertices.get(ig, aux).shape[0]
            n_edge = self.edges.get(ig, aux).shape[0]
            n_face = self.faces.get(ig, aux).shape[0]
            n_cell = self.cells.get(ig, aux).shape[0]

            self.shape[ig] = Struct(n_vertex=n_vertex,
                                    n_edge=n_edge,
                                    n_face=n_face,
                                    n_cell=n_cell)

    def setup_face_indices(self, reset=True):
        """
        Initialize an array (per group) of (iel, ifa) for each face.
        """
        if reset or not self.fis:
            fa = self.domain.get_facets(force_faces=True)[1]

            if self.faces:
                faces = self.faces
            else:
                faces = self.edges

            self.fis = {}
            for ig in self.igs:
                rfaces = faces[ig]
                fi = fa.indices[rfaces]
                assert_(nm.all(fi[:,0] == ig))
                self.fis[ig] = fi[:,1:].copy()

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
            self.true_cells[ig] = False

    def select_cells_of_surface(self, reset=True):
        """
        Select cells corresponding to faces (or edges in 2D).
        """
        if not self.can_cells:
            raise ValueError('region %s cannot have cells!' % self.name)

        self.setup_face_indices(reset=reset)

        self.cells = {}
        for ig in self.igs:
            rcells = self.fis[ig][:,0]
            self.cells[ig] = nm.ascontiguousarray(rcells)
            self.true_cells[ig] = False

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

    def create_mapping(self, kind, ig):
        """
        Create mapping from reference elements to physical elements,
        given the integration kind ('v' or 's').

        This mapping can be used to compute the physical quadrature
        points.

        Returns
        -------
        mapping : VolumeMapping or SurfaceMapping instance
            The requested mapping.
        """
        from sfepy.fem.mappings import VolumeMapping, SurfaceMapping
        from sfepy.fem.fe_surface import FESurface

        coors = self.domain.get_mesh_coors()
        if kind == 's':
            coors = coors[self.all_vertices]

        gel = self.domain.groups[ig].gel
        conn = self.domain.groups[ig].conn

        if kind == 'v':
            cells = self.cells[ig]

            mapping = VolumeMapping(coors, conn[cells], gel=gel)

        elif kind == 's':
            aux = FESurface('aux', self, gel.get_surface_entities(),
                            conn , ig)
            mapping = SurfaceMapping(coors, aux.leconn, gel=gel.surface_facet)

        return mapping

    def get_n_cells(self, ig=None, is_surface=False):
        """
        Get number of region cells.

        Parameters
        ----------
        ig : int, optional
            The group index. If None, counts from all groups are added
            together.
        is_surface : bool
            If True, number of edges or faces according to domain
            dimension is returned instead.

        Returns
        -------
        n_cells : int
            The number of cells.
        """
        if ig is not None:
            if is_surface:
                if self.domain.shape.dim == 2:
                    return self.shape[ig].n_edge

                else:
                    return self.shape[ig].n_face

            else:
                return self.shape[ig].n_cell

        else:
            return sum(self.get_n_cells(ig, is_surface=is_surface)
                       for ig in self.igs)

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

    def get_surface_entities(self, ig):
        """
        Return either region edges (in 2D) or faces (in 3D) .
        """
        if self.domain.shape.dim == 2:
            return self.edges[ig]

        else:
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
