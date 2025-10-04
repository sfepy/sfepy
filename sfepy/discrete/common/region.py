import re

import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.base.compat import in1d

_depends = re.compile(r'r\.([a-zA-Z_\-0-9.]+)').findall

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
    for sort_name, rdef in region_defs.items():
        name, sel = rdef.name, rdef.select
        if name in name_to_sort_name:
            msg = 'region %s/%s already defined!' % (sort_name, name)
            raise ValueError(msg)
        name_to_sort_name[name] = sort_name

        if name not in graph:
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
        for node, deps in graph.items():

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
            raise ValueError('circular dependency')
        idone0 = idone

    return out

def are_disjoint(r1, r2):
    """
    Check if the regions `r1` and `r2` are disjoint.

    Uses vertices for the check - `*_only` regions not allowed.
    """
    return len(nm.intersect1d(r1.vertices, r2.vertices,
                              assume_unique=True)) == 0

def _join(def1, op, def2):
    return '(' + def1 + ' ' + op + ' ' + def2 + ')'

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
        1 : {'facet' : 'vertex', 'facet_only' : 'vertex_only'},
        2 : {'facet' : 'edge',   'facet_only' : 'edge_only'},
        3 : {'facet' : 'face',   'facet_only' : 'face_only'},
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

    @staticmethod
    def from_cells(cells, domain, name='region', kind='cell', parent=None):
        """
        Create a new region containing given cells.

        Parameters
        ----------
        cells : array
            The array of cells.
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
        obj = Region(name, 'given cells', domain, '', kind=kind, parent=parent)
        obj.cells = cells

        return obj

    def __init__(self, name, definition, domain, parse_def, kind='cell',
                 parent=None, tdim=None):
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
        tdim : int
            The topological dimension of the cells.
        """
        if tdim is None:
            tdim = domain.shape.tdim
        Struct.__init__(self,
                        name=name, definition=definition,
                        domain=domain, parse_def=parse_def,
                        n_v_max=domain.shape.n_nod, dim=domain.shape.dim,
                        tdim=tdim, kind_tdim=None,
                        entities=[None] * (tdim + 1),
                        kind=None, parent=parent, shape=None,
                        mirror_regions={}, mirror_maps={}, is_empty=False,
                        cmesh=domain.cmesh_tdim[tdim])
        self.set_kind(kind)

    def set_kind(self, kind):
        if kind == self.kind: return

        if self.__facet_kinds[self.tdim]['facet'] == kind:
            kind = 'facet'

        self.kind = kind
        if 'facet' in kind:
            self.true_kind = self.__facet_kinds[self.tdim][kind]

        else:
            self.true_kind = kind

        can = [bool(ii) for ii in self.__can[self.true_kind]]

        self.can_vertices = can[0]
        self.can_cells = can[3]

        if self.tdim == 1:
            self.can = (can[0], can[3])
            self.can_edges = False
            self.can_faces = False

        elif self.tdim == 2:
            self.can = (can[0], can[1], can[3])
            self.can_edges = can[1]
            self.can_faces = False

        else:
            self.can = can
            self.can_edges = can[1]
            self.can_faces = can[2]

        for ii, ican in enumerate(self.can):
            if not ican:
                self.entities[ii] = nm.empty(0, dtype=nm.uint32)

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
        if self.tdim <= 1:
            raise AttributeError('1D region has no edges!')

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
        if self.tdim <= 2:
            raise AttributeError('1D or 2D region has no faces!')

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

        elif self.tdim == 2:
            return self.edges

        else:
            return self.vertices

    @facets.setter
    def facets(self, vals):
        if self.tdim == 3:
            self.faces = vals

        elif self.tdim == 2:
            self.edges = vals

        else:
            self.vertices = vals

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
        if dim == 0:
            self.vertices

        elif dim == 1:
            if self.tdim == 1:
                self.cells

            else:
                self.edges

        elif dim == 2:
            if self.tdim == 3:
                self.faces

            else:
                self.cells

        else:
            self.cells

    def setup_from_highest(self, dim, allow_lower=True, allow_empty=False):
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
            if not (allow_empty or self.is_empty):
                msg = 'region "%s" has no entities!'
                raise ValueError(msg % self.name)

            if self.entities[dim] is None:
                self.entities[dim] = nm.empty(0, dtype=nm.uint32)

            self.is_empty = True

            return

        cmesh = self.cmesh
        if idim <= dim:
            if not (allow_lower or allow_empty):
                msg = 'setup_from_highest() can be used only with dim < %d'
                raise ValueError(msg % idim)

            if allow_lower:
                cmesh.setup_connectivity(dim, idim)
                ents = self.get_entities(idim)
                self.entities[dim] = cmesh.get_complete(dim, ents, idim)

            else:
                for idim in range(self.kind_tdim - 1, -1, -1):
                    self.entities[idim] = nm.empty(0, dtype=nm.uint32)

                self.is_empty = True

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

        cmesh = self.cmesh
        cmesh.setup_connectivity(dim, 0)
        vv = self.vertices
        self.entities[dim] = cmesh.get_complete(dim, vv, 0)

    def finalize(self, allow_empty=False):
        """
        Initialize the entities corresponding to the region kind and regenerate
        all already existing (accessed) entities of lower topological dimension
        from the kind entities.
        """
        self._access(self.kind_tdim)
        is_empty = self.entities[self.kind_tdim].shape[0] == 0
        if allow_empty:
            self.is_empty = is_empty

        elif is_empty:
            raise ValueError('region "%s" has no entities!' % self.name)

        for idim in range(self.kind_tdim - 1, -1, -1):
            if self.can[idim] and self.entities[idim] is not None:
                try:
                    self.setup_from_highest(idim, allow_lower=False,
                                            allow_empty=allow_empty)

                except ValueError as exc:
                    msg = '\n'.join((str(exc),
                                     'fix region kind? (region: %s, kind: %s)'
                                     % (self.name, self.kind)))
                    raise ValueError(msg)

    def eval_op_vertices(self, other, op):
        parse_def = _join(self.parse_def, '%sv' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def, tdim=self.tdim)
        tmp.vertices = self.__op_to_fun[op](self.vertices, other.vertices)

        return tmp

    def eval_op_edges(self, other, op):
        parse_def = _join(self.parse_def, '%se' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def, tdim=self.tdim)
        tmp.edges = self.__op_to_fun[op](self.edges, other.edges)

        return tmp

    def eval_op_faces(self, other, op):
        parse_def = _join(self.parse_def, '%sf' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def, tdim=self.tdim)
        tmp.faces = self.__op_to_fun[op](self.faces, other.faces)

        return tmp

    def eval_op_facets(self, other, op):
        parse_def = _join(self.parse_def, '%ss' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def, tdim=self.tdim)
        tmp.facets = self.__op_to_fun[op](self.facets, other.facets)

        return tmp

    def eval_op_cells(self, other, op):
        parse_def = _join(self.parse_def, '%sc' % op, other.parse_def)
        tmp = self.light_copy('op', parse_def, tdim=self.tdim)
        tmp.cells = self.__op_to_fun[op](self.cells, other.cells)

        return tmp

    def light_copy(self, name, parse_def, tdim=None):
        return Region(name, self.definition, self.domain, parse_def,
                      kind=self.kind, tdim=tdim)

    def copy(self):
        """
        Make a copy based on the region kind.
        """
        tmp = self.light_copy('copy', self.parse_def, tdim=self.tdim)
        tmp.entities[self.kind_tdim] = self.get_entities(self.kind_tdim).copy()

        return tmp

    def delete_zero_faces(self, eps=1e-14):
        raise NotImplementedError

    def update_shape(self):
        """
        Update shape of each group according to region vertices, edges,
        faces and cells.
        """
        n_vertex = self.vertices.shape[0]
        n_cell = self.cells.shape[0]
        n_edge = self.edges.shape[0] if self.tdim > 1 else 0
        n_face = self.faces.shape[0] if self.tdim == 3 else 0
        n_facet = self.facets.shape[0]

        self.shape = Struct(n_vertex=n_vertex,
                            n_edge=n_edge,
                            n_face=n_face,
                            n_facet=n_facet,
                            n_cell=n_cell)

    def get_entities(self, dim):
        """
        Return mesh entities of dimension `dim`.
        """
        if dim <= self.tdim:
            self._access(dim)
            out = self.entities[dim]

        else:
            out = nm.empty(0, dtype=nm.uint32)

        return out

    def get_cells(self, true_cells_only=True):
        """
        Get cells of the region.

        Raises ValueError if `true_cells_only` is True and the region kind does
        not allow cells. For `true_cells_only` equal to False, cells incident
        to facets are returned if the region itself contains no cells. Obeys
        parent region, if given.
        """
        if self.kind != 'cell':
            if true_cells_only:
                msg = 'region %s has not true cells! (has kind: %s)' \
                      % (self.name, self.kind)
                raise ValueError(msg)

            else:
                # Has to be consistent with get_facet_indices()!
                cmesh = self.cmesh
                cmesh.setup_connectivity(self.tdim - 1, self.tdim)
                out = cmesh.get_incident(self.tdim, self.facets, self.tdim - 1)

                if self.parent is not None:
                    pcells = self.domain.regions[self.parent].cells
                    ip = in1d(out, pcells, assume_unique=False)
                    out = out[ip]

        else:
            out = self.cells

        return out

    def get_cell_indices(self, cells, true_cells_only=True):
        """
        Return indices of `cells` in the region cells.

        Raises ValueError if `true_cells_only` is True and the region kind does
        not allow cells. For `true_cells_only` equal to False, cells incident
        to facets are returned if the region itself contains no cells.

        Raises ValueError if all `cells` are not in the region cells.
        """
        fcells = self.get_cells(true_cells_only=true_cells_only)

        iin = nm.isin(cells, fcells)
        if not iin.all():
            raise ValueError(f'some cells are not in region {self.name} cells!')

        iis = nm.searchsorted(fcells, cells)
        assert_((fcells[iis[iin]] == cells[iin]).all())

        ii = nm.where(iin, iis, -1)
        return ii

    def get_facet_indices(self):
        """
        Return an array (per group) of (iel, ifa) for each facet. A facet can
        be in 1 (surface) or 2 (inner) cells.
        """
        cmesh = self.cmesh
        cmesh.setup_connectivity(self.tdim - 1, self.tdim)

        facets = self.facets
        cells, offs = cmesh.get_incident(self.tdim, facets, self.tdim - 1,
                                         ret_offsets=True)

        if self.parent is not None:
            pcells = self.domain.regions[self.parent].cells
            ip = in1d(cells, pcells, assume_unique=False)
            cells = cells[ip]

            counts = nm.diff(offs).astype(nm.int32)
            pos = nm.repeat(nm.arange(facets.shape[0], dtype=nm.int32), counts)
            new_counts = nm.bincount(pos, weights=ip).astype(nm.uint32)
            offs = nm.cumsum(nm.r_[0, new_counts], dtype=nm.uint32)

        ii = cmesh.get_local_ids(facets, self.tdim - 1, cells, offs, self.tdim)
        fis = nm.c_[cells, ii]

        return fis

    def setup_mirror_region(self, mirror_name=None, ret_name=False):
        """
        Find the corresponding mirror region, set up element mapping.
        """
        from sfepy.discrete.fem.mesh import find_map

        regions = self.domain.regions
        eopts = self.extra_options

        if (mirror_name is None) and (eopts is not None)\
            and ('mirror_region' in eopts):
            mirror_name = eopts['mirror_region']

        if mirror_name is not None:
            if mirror_name in self.mirror_regions:
                return mirror_name if ret_name else None

            mreg = regions[mirror_name]
            if self.vertices.shape[0] != mreg.vertices.shape[0]:
                raise ValueError('%s: incompatible mirror region! (%s)'
                    % (self.name, mreg.name))
            coors = self.cmesh.coors
            coors1 = coors[self.vertices, :]
            coors2 = coors[mreg.vertices, :]
            shift = ((nm.sum(coors2, axis=0) - nm.sum(coors1, axis=0))
                / coors1.shape[0])
            vmap1, vmap2 = find_map(coors1, coors2 - shift, join=False)
            if vmap1.shape[0] != coors1.shape[0]:
                print(coors1[vmap1])
                print(coors2[vmap2])
                raise ValueError('cannot match vertices!')

            _ = self.cmesh.get_conn_as_graph(self.dim - 1, 0)
            cc1 = self.cmesh.get_centroids(self.dim - 1)
            sfacets = self.cells if self.tdim < self.dim else self.facets
            cc2 = mreg.cmesh.get_centroids(mreg.dim - 1)
            mfacets = mreg.cells if mreg.tdim < mreg.dim else mreg.facets
            fmap1, fmap2 = find_map(cc1[sfacets], cc2[mfacets] - shift,
                                    join=False)
            if fmap1.shape[0] != sfacets.shape[0]:
                print(cc1[sfacets][fmap1])
                print(cc2[mfacets][fmap2])
                raise ValueError('cannot match facets!')

            mirror_map_i = nm.zeros_like(fmap1)
            mirror_map_i[fmap1] = fmap2
            mirror_map = nm.zeros_like(fmap1)
            mirror_map[fmap2] = fmap1
        else:
            for reg in regions:
                mirror_parent = regions.find(reg.parent)
                if mirror_parent is None: continue
                if ((reg is not self)
		    and (len(self.vertices)) == len(reg.vertices)
                    and nm.all(self.vertices == reg.vertices)):
                    mreg = reg
                    mirror_map = mirror_map_i = None
                    break
            else:
                raise ValueError('cannot find mirror region! (%s)' % self.name)

        self.mirror_regions[mreg.name] = mreg
        self.mirror_maps[mreg.name] = mirror_map
        mreg.mirror_regions[self.name] = self
        mreg.mirror_maps[self.name] = mirror_map_i

        if ret_name:
            return mreg.name

    def get_mirror_region(self, name):
        return self.mirror_regions[name]

    def get_n_cells(self, is_surface=False):
        """
        Get number of region cells.

        Parameters
        ----------
        is_surface : bool
            If True, number of edges or faces according to domain
            dimension is returned instead.

        Returns
        -------
        n_cells : int
            The number of cells.
        """
        if is_surface:
            return self.shape.n_facet

        else:
            return self.shape.n_cell

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

        cmesh = self.cmesh

        e_verts = cmesh.get_incident(0, self.edges, 1)
        e_verts.shape = (e_verts.shape[0] // 2, 2)

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
