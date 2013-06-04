import re
from copy import copy

import numpy as nm

from sfepy.base.base import assert_, Struct

_depends = re.compile('r\.([a-zA-Z_0-9.]+)').findall

def get_parents(selector):
    """
    Given a region selector, return names of regions it is based on.
    """
    parents = _depends(selector)

    return parents

def get_dependency_graph(region_defs):
    """
    Return a dependency graph and a name-sort name mapping for given
    region definitions.
    """
    graph = {}
    name_to_sort_name = {}
    for sort_name, rdef in region_defs.iteritems():
        name, sel = rdef.name, rdef.select
        if name_to_sort_name.has_key(name):
            msg = 'region %s/%s already defined!' % (sort_name, name)
            raise ValueError(msg)
        name_to_sort_name[name] = sort_name

        if not graph.has_key(name):
            graph[name] = [0]

        for parent in get_parents(sel):
            graph[name].append(parent)

    return graph, name_to_sort_name

def sort_by_dependency(graph):
    out = []

    n_nod = len(graph)
    idone = 0
    idone0 = -1
    while idone < n_nod:

        dep_removed = 0
        for node, deps in graph.iteritems():

            if (len(deps) == 1) and not deps[0]:
                out.append(node)
                deps[0] = 1
                idone += 1

            elif not deps[0]:

                for ii, dep in enumerate(deps[1:]):
                    if graph[dep][0]:
                        ir = deps.index(dep)
                        deps.pop(ir)
                        dep_removed += 1

        if (idone <= idone0) and not dep_removed:
            raise ValueError, 'circular dependency'
        idone0 = idone

    return out

def _join(def1, op, def2):
    return '(' + def1 + ' ' + op + ' ' + def2 + ')'

def _try_delete(obj, ig):
    try:
        del obj[ig]
    except KeyError:
        pass

class Region(Struct):
    """
    Region defines a subset of a FE domain.

    Created: 31.10.2005
    """
    __can = {
        'cell'        : (1, 1, 1, 1),
        'face'        : (1, 1, 1, 0),
        'edge'        : (1, 1, 0, 0),
        'vertex'      : (1, 0, 0, 0),
        'cell_only'   : (0, 0, 0, 1),
        'face_only'   : (0, 0, 1, 0),
        'edge_only'   : (0, 1, 0, 0),
        'vertex_only' : (1, 0, 0, 0),
    }

    __facet_kinds = {
        2 : {'facet' : 'edge', 'facet_only' : 'edge_only'},
        3 : {'facet' : 'face', 'facet_only' : 'face_only'},
    }

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

    def __init__(self, name, definition, domain, parse_def, kind='cell'):
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
                        domain=domain, parse_def=parse_def,
                        n_v_max=domain.shape.n_nod, dim=domain.shape.dim,
                        igs=[], entities=[None] * 4,
                        fis={},
                        must_update=True,
                        is_complete=False,
                        mirror_region=None, ig_map=None,
                        ig_map_i=None)
        self.set_kind(kind)

    def set_kind(self, kind):
        self.kind = kind
        if 'facet' in kind:
            self.true_kind = self.__facet_kinds[self.dim][kind]

        else:
            self.true_kind = kind

        can = [bool(ii) for ii in self.__can[self.true_kind]]
        can[2] = can[2] and (self.dim == 3)
        self.can = can

        self.can_vertices = can[0]
        self.can_edges = can[1]
        self.can_faces = can[2]
        self.can_cells = can[3]

        for ii, ican in enumerate(self.can):
            if not ican:
                self.entities[ii] = nm.empty(0, dtype=nm.uint32)

    def light_copy(self, name, parse_def):
        return Region(name, self.definition, self.domain, parse_def,
                      kind=self.kind)

    def setup_from_highest(self, dim):
        """
        Setup entities of topological dimension `dim` using the available
        entities of the highest topological dimension.
        """
        if not self.can[dim]: return

        for idim in range(self.dim + 1):
            if self.entities[idim] is not None:
                if self.entities[idim].shape[0] > 0:
                    break

        if idim <= dim:
            msg = 'setup_from_highest() can be used only with dim < %d'
            raise ValueError(msg % idim)

        cmesh = self.domain.cmesh
        cmesh.setup_connectivity(idim, dim)

        incident = cmesh.get_incident(dim, self.entities[idim], idim)
        self.entities[dim] = nm.unique(incident)

    def setup_from_vertices(self, dim):
        """
        Setup entities of topological dimension `dim` using the region
        vertices.
        """
        if not self.can[dim]: return

        cmesh = self.domain.cmesh
        cmesh.setup_connectivity(dim, 0)
        vv = self.vertices
        self.entities[dim] = cmesh.get_complete(dim, vv, 0)

    @property
    def vertices(self):
        if self.entities[0] is None:
            self.setup_from_highest(0)
        return self.entities[0]

    @vertices.setter
    def vertices(self, vals):
        if self.can_vertices:
            self.entities[0] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have vertices!' % self.name)

    @property
    def edges(self):
        if self.entities[1] is None:
            self.setup_from_vertices(1)
        return self.entities[1]

    @edges.setter
    def edges(self, vals):
        if self.can_edges:
            self.entities[1] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have edges!' % self.name)

    @property
    def faces(self):
        if self.dim == 2:
            raise AttributeError('2D region has no faces!')

        if self.entities[2] is None:
            self.setup_from_vertices(2)
        return self.entities[2]

    @faces.setter
    def faces(self, vals):
        if self.can_faces:
            self.entities[2] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have faces!' % self.name)

    @property
    def facets(self):
        if self.dim == 3:
            return self.faces

        else:
            return self.edges

    @facets.setter
    def facets(self, vals):
        if self.dim == 3:
            self.faces = vals

        else:
            self.edges = vals

    @property
    def cells(self):
        if self.entities[self.dim] is None:
            self.setup_from_vertices(self.dim)
        return self.entities[self.dim]

    @cells.setter
    def cells(self, vals):
        if self.can_cells:
            self.entities[self.dim] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have cells!' % self.name)

    def delete_zero_faces(self, eps=1e-14):
        """
        Delete faces with zero area.
        """
        from sfepy.linalg import get_face_areas

        fa = self.domain.fa

        for ig, faces in self.faces.iteritems():
            fav = fa.facets[faces]

            areas = get_face_areas(fav, self.domain.mesh.coors)
            self.faces[ig] = faces[areas > eps]

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

    def select_cells(self, n_verts):
        """
        Select cells containing at least n_verts[ii] vertices per group ii.
        """
        if not self.can_cells:
            raise ValueError('region %s cannot have cells!' % self.name)

        self.cells = {}
        for ig, group in self.domain.iter_groups(self.igs):
            vv = self.vertices[ig]
            if len(vv) == 0: continue

            mask = nm.zeros(self.n_v_max, nm.int32)
            mask[vv] = 1

            aux = nm.sum(mask[group.conn], 1)
            rcells = nm.where(aux >= n_verts[ig])[0]
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

    def copy(self):
        """
        Vertices-based copy.
        """
        tmp = self.light_copy('copy', self.parse_def)
        tmp.set_vertices(copy(self.all_vertices))
        return tmp

    def sub_n(self, other):
        tmp = self.light_copy('op',
                              _join(self.parse_def, '-n', other.parse_def))
        tmp.set_vertices(nm.setdiff1d(self.all_vertices, other.all_vertices))

        return tmp

    def add_n(self, other):
        tmp = self.light_copy('op',
                              _join(self.parse_def, '+n', other.parse_def))
        tmp.set_vertices(nm.union1d(self.all_vertices, other.all_vertices))

        return tmp

    def intersect_n(self, other):
        tmp = self.light_copy('op',
                              _join(self.parse_def, '*n', other.parse_def))
        tmp.set_vertices(nm.intersect1d(self.all_vertices, other.all_vertices))

        return tmp

    def sub_e(self, other):
        tmp = self.light_copy('op',
                              _join(self.parse_def, '-e', other.parse_def))
        for ig in self.igs:
            if ig not in other.igs:
                tmp.igs.append(ig)
                tmp.cells[ig] = self.cells[ig].copy()
                continue

            aux = nm.setdiff1d(self.cells[ig], other.cells[ig])
            if not len(aux): continue
            tmp.cells[ig] = aux
            tmp.igs.append(ig)

        tmp.update_vertices()
        return tmp

    def add_e(self, other):
        tmp = self.light_copy('op',
                              _join(self.parse_def, '+e', other.parse_def))
        for ig in self.igs:
            tmp.igs.append(ig)
            if ig not in other.igs:
                tmp.cells[ig] = self.cells[ig].copy()
                continue

            tmp.cells[ig] = nm.union1d(self.cells[ig], other.cells[ig])

        for ig in other.igs:
            if ig in tmp.igs: continue
            tmp.igs.append(ig)
            tmp.cells[ig] = other.cells[ig].copy()

        tmp.update_vertices()
        return tmp

    def intersect_e(self, other):
        tmp = self.light_copy('op',
                              _join(self.parse_def, '*e', other.parse_def))
        for ig in self.igs:
            if ig not in other.igs: continue
            aux = nm.intersect1d(self.cells[ig], other.cells[ig])
            if not len(aux): continue
            tmp.igs.append(ig)
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
                raise ValueError('cannot find mirror region group! (%d)' % igr)

        self.mirror_region = mirror_region
        self.ig_map = ig_map
        self.ig_map_i = ig_map_i

        if self.domain.shape.dim == 2:
            self.domain.ed.setup_group_interfaces()

        elif self.domain.shape.dim == 3:
            self.domain.fa.setup_group_interfaces()

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
                if self.domain.groups[ig].shape.dim == 2:
                    return self.shape[ig].n_edge

                else:
                    return self.shape[ig].n_face

            else:
                return self.shape[ig].n_cell

        else:
            return sum(self.get_n_cells(ig, is_surface=is_surface)
                       for ig in self.igs)

    def get_vertices_of_cells(self, return_per_group=False):
        """
        Return all vertices, that are in some cell of the region.

        Parameters
        ----------
        return_per_group : bool
            It True, return also a dict of vertices per element group.
        """
        all_vertices = nm.zeros((0,), dtype=nm.int32)
        vertices = {}
        for ig, group in self.domain.iter_groups(self.igs):
            rcells = self.cells[ig]
            conn = group.conn
            nods = conn[rcells,:].ravel()
            vertices[ig] = nm.unique(nods)
            all_vertices = nm.unique(nm.r_[all_vertices, vertices[ig]])

        if return_per_group:
            out = all_vertices, vertices

        else:
            out = all_vertices

        return out

    def get_vertices(self, ig):
        return self.vertices[ig]

    def get_edges(self, ig):
        return self.edges[ig]

    def get_faces(self, ig):
        return self.faces[ig]

    def get_surface_entities(self, ig):
        """
        Return either region edges (in 2D) or faces (in 3D) .
        """
        if self.domain.shape.dim == 2:
            return self.edges[ig]

        else:
            return self.faces[ig]

    def get_cells(self, ig, true_cells_only=True):
        """
        Get cells of the region.

        Raises ValueError if true_cells_only is True and the cells are not true
        cells (e.g. surface integration region).
        """
        if true_cells_only and not self.true_cells[ig]:
            msg = 'region %s has not true cells! (surface integration?)' \
                  % self.name
            raise ValueError(msg)
        return self.cells[ig]

    def iter_cells(self):
        ii = 0
        for ig, cells in self.cells.iteritems():
            for iel in cells:
                yield ig, ii, iel
                ii += 1

    def has_cells(self):

        if self.can_cells:
            for cells in self.cells.itervalues():
                if cells.size:
                    return True
            return False
        else:
            return False

    def has_cells_if_can(self):
        if self.can_cells:
            for cells in self.cells.itervalues():
                if cells.size:
                    return True
            return False
        else:
            return True

    def contains(self, other):
        """
        Tests only igs for now!!!
        """
        return set(other.igs).issubset(set(self.igs))

    def get_cell_offsets(self):
        offs = {}
        off = 0
        for ig in self.igs:
            offs[ig] = off
            off += self.shape[ig].n_cell
        return offs

    def get_charfun(self, by_cell=False, val_by_id=False):
        """
        Return the characteristic function of the region as a vector of values
        defined either in the mesh nodes (by_cell == False) or cells. The
        values are either 1 (val_by_id == False) or sequential id + 1.
        """
        if by_cell:
            chf = nm.zeros((self.domain.shape.n_el,), dtype=nm.float64)
            offs = self.get_cell_offsets()
            for ig, cells in self.cells.iteritems():
                iel = offs[ig] + cells
                if val_by_id:
                    chf[iel] = iel + 1
                else:
                    chf[iel] = 1.0
        else:
            chf = nm.zeros((self.domain.shape.n_nod,), dtype=nm.float64)
            if val_by_id:
                chf[self.all_vertices] = self.all_vertices + 1
            else:
                chf[self.all_vertices] = 1.0

        return chf

    def get_edge_graph(self):
        """
        Return the graph of region edges as a sparse matrix having uid(k) + 1
        at (i, j) if vertex[i] is connected with vertex[j] by the edge k.

        Degenerate edges are ignored.
        """
        from scipy.sparse import csr_matrix

        ed = self.domain.ed

        rows, cols, vals = [], [], []
        for ig, edges in self.edges.iteritems():
            e_verts = ed.facets[edges]
            ii = nm.where(e_verts[:, 0] != e_verts[:, 1])[0]
            edges = edges[ii]
            e_verts = e_verts[ii]

            vals.append(ed.uid_i[edges] + 1)
            rows.append(e_verts[:, 0])
            cols.append(e_verts[:, 1])

        vals, indx = nm.unique(nm.concatenate(vals), return_index=True)
        rows = nm.concatenate(rows)[indx]
        cols = nm.concatenate(cols)[indx]

        num = self.all_vertices.max() + 1
        graph = csr_matrix((vals, (rows, cols)), shape=(num, num))

        nnz = graph.nnz
        # Symmetrize.
        graph = graph + graph.T
        assert_(graph.nnz == 2 * nnz)

        return graph
