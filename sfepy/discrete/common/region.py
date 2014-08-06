import re
from copy import copy

import numpy as nm

from sfepy.base.base import output, assert_, Struct

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

        if rdef.get('parent', None) is not None:
            graph[name].append(rdef.parent)

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
                    if not dep in graph:
                        msg = 'dependency %s of region %s does not exist!'
                        raise ValueError(msg % (dep, node))

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

    Region kinds:

    - cell_only, facet_only, face_only, edge_only, vertex_only - only the
      specified entities are included, others are empty sets (so that the
      operators are still defined)
    - cell, facet, face, edge, vertex - entities of higher dimension are not
      included

    The 'cell' kind is the most general and it is the default.

    Region set-like operators: + (union), - (difference), * (intersection),
    followed by one of ('v', 'e', 'f', 'c', and 's') for vertices, edges,
    faces, cells, and facets.

    Notes
    -----
    Functions depending on `ig` are adapters for current code that should be
    removed after new assembling is done.

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

    __op_to_fun = {
        '+' : nm.union1d,
        '-' : nm.setdiff1d,
        '*' : nm.intersect1d,
    }

    @staticmethod
    def from_vertices(vertices, domain, name='region', kind='cell'):
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
        kind : str, optional
            The kind of the region.

        Returns
        -------
        obj : Region instance
            The new region.
        """
        obj = Region(name, 'given vertices', domain, '', kind=kind)
        obj.vertices = vertices

        return obj

    @staticmethod
    def from_facets(facets, domain, name='region', kind='facet', parent=None):
        """
        Create a new region containing given facets.

        Parameters
        ----------
        facets : array
            The array with indices to unique facets.
        domain : Domain instance
            The domain containing the facets.
        name : str, optional
            The name of the region.
        kind : str, optional
            The kind of the region.
        parent : str, optional
            The name of the parent region.

        Returns
        -------
        obj : Region instance
            The new region.
        """
        obj = Region(name, 'given faces', domain, '', kind=kind, parent=parent)
        obj.facets = facets

        return obj

    def __init__(self, name, definition, domain, parse_def, kind='cell',
                 parent=None):
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
        kind : str
            The region kind - one of 'cell', 'facet', 'face', 'edge', 'vertex',
            'cell_only', ..., 'vertex_only'.
        parent : str, optional
            The name of the parent region.
        """
        tdim = domain.shape.tdim
        Struct.__init__(self,
                        name=name, definition=definition,
                        domain=domain, parse_def=parse_def,
                        n_v_max=domain.shape.n_nod, dim=domain.shape.dim,
                        tdim=tdim, kind_tdim=None,
                        entities=[None] * (tdim + 1),
                        kind=None, parent=parent, shape=None,
                        mirror_region=None, ig_map=None,
                        ig_map_i=None)
        self.set_kind(kind)

    def set_kind(self, kind):
        if kind == self.kind: return

        self.kind = kind
        if 'facet' in kind:
            self.true_kind = self.__facet_kinds[self.tdim][kind]

        else:
            self.true_kind = kind

        can = [bool(ii) for ii in self.__can[self.true_kind]]

        self.can_vertices = can[0]
        self.can_edges = can[1]

        if self.tdim == 2:
            self.can = (can[0], can[1], can[3])
            self.can_cells = can[2]

        else:
            self.can = can
            self.can_faces = can[2]
            self.can_cells = can[3]

        for ii, ican in enumerate(self.can):
            if not ican:
                self.entities[ii] = nm.empty(0, dtype=nm.uint32)

        self._igs = None

        self.set_kind_tdim()

    def set_kind_tdim(self):
        if 'vertex' in self.true_kind:
            self.kind_tdim = 0

        elif 'edge' in self.true_kind:
            self.kind_tdim = 1

        elif 'face' in self.true_kind:
            self.kind_tdim = 2

        elif 'cell' in self.true_kind:
            self.kind_tdim = self.tdim

    @property
    def vertices(self):
        if self.entities[0] is None:
            self._access(1)
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
            if 'edge' in self.true_kind:
                self.setup_from_vertices(1)

            else:
                self._access(2)
                self.setup_from_highest(1)

        return self.entities[1]

    @edges.setter
    def edges(self, vals):
        if self.can_edges:
            self.entities[1] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have edges!' % self.name)

    @property
    def faces(self):
        if self.tdim == 2:
            raise AttributeError('2D region has no faces!')

        if self.entities[2] is None:
            if 'face' in self.true_kind:
                self.setup_from_vertices(2)

            else:
                self._access(3)
                self.setup_from_highest(2)

        return self.entities[2]

    @faces.setter
    def faces(self, vals):
        if self.can_faces:
            self.entities[2] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have faces!' % self.name)

    @property
    def facets(self):
        if self.tdim == 3:
            return self.faces

        else:
            return self.edges

    @facets.setter
    def facets(self, vals):
        if self.tdim == 3:
            self.faces = vals

        else:
            self.edges = vals

    @property
    def cells(self):
        if self.entities[self.tdim] is None:
            self.setup_from_vertices(self.tdim)
        return self.entities[self.tdim]

    @cells.setter
    def cells(self, vals):
        if self.can_cells:
            self.entities[self.tdim] = nm.asarray(vals, dtype=nm.uint32)

        else:
            raise ValueError('region "%s" cannot have cells!' % self.name)

    def _access(self, dim):
        """
        Helper to access region entities of dimension `dim`.
        """
        if dim == 1:
            self.edges

        elif dim == 2:
            if self.tdim == 3:
                self.faces

            else:
                self.cells

        else:
            self.cells

    def setup_from_highest(self, dim, allow_lower=True):
        """
        Setup entities of topological dimension `dim` using the available
        entities of the highest topological dimension.
        """
        if not self.can[dim]: return

        for idim in range(self.tdim, -1, -1):
            if self.entities[idim] is not None:
                if self.entities[idim].shape[0] > 0:
                    break

        else:
            msg = 'region "%s" has no entities!'
            raise ValueError(msg % self.name)

        cmesh = self.domain.cmesh
        if idim <= dim:
            if not allow_lower:
                msg = 'setup_from_highest() can be used only with dim < %d'
                raise ValueError(msg % idim)

            cmesh.setup_connectivity(dim, idim)
            ents = self.get_entities(idim)
            self.entities[dim] = cmesh.get_complete(dim, ents, idim)

        else:
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

    def finalize(self):
        """
        Initialize the entities corresponding to the region kind and regenerate
        all already existing (accessed) entities of lower topological dimension
        from the kind entities.
        """
        self._access(self.kind_tdim)
        for idim in range(self.kind_tdim - 1, -1, -1):
            if self.can[idim] and self.entities[idim] is not None:
                try:
                    self.setup_from_highest(idim, allow_lower=False)

                except ValueError, exc:
                    msg = '\n'.join((exc.message,
                                     'fix region kind? (region: %s, kind: %s)'
                                     % (self.name, self.kind)))
                    raise ValueError(msg)

    def eval_op_vertices(self, other, op):
        parse_def = _join(self.parse_def, '%sv' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def)
        tmp.vertices = self.__op_to_fun[op](self.vertices, other.vertices)

        return tmp

    def eval_op_edges(self, other, op):
        parse_def = _join(self.parse_def, '%se' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def)
        tmp.edges = self.__op_to_fun[op](self.edges, other.edges)

        return tmp

    def eval_op_faces(self, other, op):
        parse_def = _join(self.parse_def, '%sf' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def)
        tmp.faces = self.__op_to_fun[op](self.faces, other.faces)

        return tmp

    def eval_op_facets(self, other, op):
        parse_def = _join(self.parse_def, '%ss' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def)
        tmp.facets = self.__op_to_fun[op](self.facets, other.facets)

        return tmp

    def eval_op_cells(self, other, op):
        parse_def = _join(self.parse_def, '%sc' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def)
        tmp.cells = self.__op_to_fun[op](self.cells, other.cells)

        return tmp

    def light_copy(self, name, parse_def):
        return Region(name, self.definition, self.domain, parse_def,
                      kind=self.kind)

    def copy(self):
        """
        Vertices-based copy.
        """
        tmp = self.light_copy('copy', self.parse_def)
        tmp.vertices = copy(self.vertices)

        return tmp

    def delete_zero_faces(self, eps=1e-14):
        raise NotImplementedError

    @property
    def igs(self):
        """
        Cell group indices according to region kind.
        """
        if self.parent is not None:
            self._igs = self.domain.regions[self.parent].igs

        elif self._igs is None:
            if 'vertex' in self.true_kind:
                self._igs = self.domain.cmesh.get_igs(self.vertices, 0)

            elif 'edge' in self.true_kind:
                self._igs = self.domain.cmesh.get_igs(self.edges, 1)

            elif 'face' in self.true_kind:
                self._igs = self.domain.cmesh.get_igs(self.faces, 2)

            elif 'cell' in self.true_kind:
                self._igs = self.domain.cmesh.get_igs(self.cells, self.tdim)

            if not len(self._igs):
                output('warning: region %s of %s kind has empty group indices!'
                       % (self.name, self.kind))

        return self._igs

    def update_shape(self):
        """
        Update shape of each group according to region vertices, edges,
        faces and cells.
        """
        get = self.domain.cmesh.get_from_cell_group

        self.shape = {}
        for ig in self.igs:
            n_vertex = get(ig, 0, self.vertices).shape[0]
            n_edge = get(ig, 1, self.edges).shape[0]
            n_cell = get(ig, self.tdim, self.cells).shape[0]
            if self.tdim == 3:
                n_face = get(ig, 2, self.faces).shape[0]

            else:
                n_face = 0

            self.shape[ig] = Struct(n_vertex=n_vertex,
                                    n_edge=n_edge,
                                    n_face=n_face,
                                    n_cell=n_cell)

    def get_entities(self, dim, ig=None):
        """
        Return mesh entities of dimension `dim`, and optionally with the cell
        group `ig`.
        """
        if ig is not None:
            out = self.domain.cmesh.get_from_cell_group(ig, dim,
                                                        self.entities[dim])

        else:
            if dim <= self.tdim:
                self._access(dim)
                out = self.entities[dim]

            else:
                out = nm.empty(0, dtype=nm.uint32)

        return out

    def get_vertices_of_cells(self):
        """
        Return all vertices, that are in some cell of the region.
        """
        vertices = self.domain.cmesh.get_incident(0, self.cells, self.tdim)

        return nm.unique(vertices)

    def get_vertices(self, ig):
        out = self.domain.cmesh.get_from_cell_group(ig, 0, self.vertices)
        return out

    def get_edges(self, ig):
        out = self.domain.cmesh.get_from_cell_group(ig, 1, self.edges)
        return out

    def get_faces(self, ig):
        out = self.domain.cmesh.get_from_cell_group(ig, 2, self.faces)
        return out

    def get_facets(self, ig):
        """
        Return either region edges (in 2D) or faces (in 3D) .
        """
        if self.tdim == 2:
            return self.get_edges(ig)

        else:
            return self.get_faces(ig)

    def get_cells(self, ig, true_cells_only=True, offset=True):
        """
        Get cells of the region.

        Raises ValueError if `true_cells_only` is True and the region kind does
        not allow cells (e.g. surface integration region). For
        `true_cells_only` equal to False, cells incident to facets are returned
        if the region itself contains no cells.

        If `offset` is True, the cell group offset is subtracted from the cell
        ids.
        """
        cmesh = self.domain.cmesh

        if self.cells.shape[0] == 0:
            if true_cells_only:
                msg = 'region %s has not true cells! (has kind: %s)' \
                      % (self.name, self.kind)
                raise ValueError(msg)

            else:
                # Has to be consistent with get_facet_indices()!
                cmesh.setup_connectivity(self.tdim - 1, self.tdim)
                out = cmesh.get_incident(self.tdim, self.facets, self.tdim - 1)

                igs = cmesh.cell_groups[out]
                ic = nm.where(igs == ig)
                out = out[ic]

        else:
            out = cmesh.get_from_cell_group(ig, self.tdim, self.cells)

        if offset:
            out -= self.domain.cell_offsets[ig]

        return out

    def get_facet_indices(self, ig, offset=True, force_ig=True):
        """
        Return an array (per group) of (iel, ifa) for each facet. A facet can
        be in 1 (surface) or 2 (inner) cells.

        If `offset` is True, the cell group offset is subtracted from the cell
        ids.

        If `force_ig` is True, only the cells with the given `ig` are used.
        """
        cmesh = self.domain.cmesh
        facets = self.get_facets(ig)
        cells, offs = cmesh.get_incident(self.tdim, facets, self.tdim - 1,
                                         ret_offsets=True)
        ii = cmesh.get_local_ids(facets, self.tdim - 1, cells, offs, self.tdim)
        fis = nm.c_[cells, ii]

        if force_ig:
            igs = cmesh.cell_groups[cells]
            ic = nm.where(igs == ig)
            fis = fis[ic]

        if offset:
            fis[:, 0] -= self.domain.cell_offsets[ig]

        return fis

    def setup_mirror_region(self):
        """
        Find the corresponding mirror region, set up element mapping.
        """
        regions = self.domain.regions

        parent = regions[self.parent]
        for reg in regions:
            mirror_parent = regions.find(reg.parent)
            if mirror_parent is None: continue
            if ((reg is not self)
                and (len(reg.igs) == len(self.igs))
                and (not len(nm.intersect1d(parent.igs, mirror_parent.igs)))
                and nm.all(self.vertices == reg.vertices)):
                mirror_region = reg
                break
        else:
            raise ValueError('cannot find mirror region! (%s)' % self.name)

        ig_map = {}
        ig_map_i = {}
        for igr in parent.igs:
            v1 = self.get_vertices(igr)
            for igc in mirror_parent.igs:
                v2 = self.get_vertices(igc)
                if nm.all(v1 == v2):
                    ig_map[igc] = igr
                    ig_map_i[igr] = igc
                    break
            else:
                raise ValueError('cannot find mirror region group! (%d)' % igr)

        self._igs = parent._igs
        mirror_region._igs = mirror_parent._igs

        self.mirror_region = mirror_region
        self.ig_map = ig_map
        self.ig_map_i = ig_map_i

    def get_mirror_region(self):
        return self.mirror_region, self.ig_map, self.ig_map_i

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

    def iter_cells(self):
        ii = 0
        off = 0
        for ig in self.igs:
            n_cell = self.shape[ig].n_cell
            for iel in self.cells[off : off + n_cell]:
                yield ig, ii, iel - off
                ii += 1

            off += n_cell

    def has_cells(self):
        return self.cells.size > 0

    def contains(self, other):
        """
        Return True in the region contains the `other` region.

        The check is performed using entities corresponding to the other region
        kind.
        """
        tdim = other.kind_tdim

        se = self.get_entities(tdim)
        oe = other.entities[tdim]

        return len(nm.intersect1d(se, oe))

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
        defined either in the mesh vertices (by_cell == False) or cells. The
        values are either 1 (val_by_id == False) or sequential id + 1.
        """
        if by_cell:
            chf = nm.zeros((self.domain.shape.n_el,), dtype=nm.float64)
            if val_by_id:
                chf[self.cells] = self.cells + 1
            else:
                chf[self.cells] = 1.0

        else:
            chf = nm.zeros((self.domain.shape.n_nod,), dtype=nm.float64)
            if val_by_id:
                chf[self.vertices] = self.vertices + 1
            else:
                chf[self.vertices] = 1.0

        return chf

    def get_edge_graph(self):
        """
        Return the graph of region edges as a sparse matrix having uid(k) + 1
        at (i, j) if vertex[i] is connected with vertex[j] by the edge k.

        Degenerate edges are ignored.
        """
        from scipy.sparse import csr_matrix

        cmesh = self.domain.cmesh

        e_verts = cmesh.get_incident(0, self.edges, 1)
        e_verts.shape = (e_verts.shape[0] / 2, 2)

        ii = nm.where(e_verts[:, 0] != e_verts[:, 1])[0]
        edges = self.edges[ii]
        e_verts = e_verts[ii]

        vals = edges + 1
        rows = e_verts[:, 0]
        cols = e_verts[:, 1]

        num = self.vertices.max() + 1
        graph = csr_matrix((vals, (rows, cols)), shape=(num, num))

        nnz = graph.nnz
        # Symmetrize.
        graph = graph + graph.T
        assert_(graph.nnz == 2 * nnz)

        return graph
