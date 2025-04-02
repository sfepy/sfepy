# -*- Mode: Python -*-
# cython: language_level=3
"""
C Mesh data structures and functions.
"""
cimport cython

import numpy as np

np.import_array()

cdef class CConnectivity:
    """
    Notes
    -----

    The memory is allocated/freed in C - this class just wraps NumPy arrays
    around that data without copying.
    """
    def __cinit__(self, num, n_incident):
        self.num = num
        self.n_incident = n_incident

        self.offset = 0

    cdef _set_conn(self, MeshConnectivity *conn):
        # This cannot be in __cinit__, as a C pointer needs to be passed
        # around.
        cdef np.npy_intp shape[1]

        self.conn = conn

        shape[0] = <np.npy_intp> self.num + 1
        self.offsets = np.PyArray_SimpleNewFromData(1, shape,
                                                    np.NPY_UINT32,
                                                    <void *> conn.offsets)

        shape[0] = <np.npy_intp> self.n_incident
        self.indices = np.PyArray_SimpleNewFromData(1, shape,
                                                    np.NPY_UINT32,
                                                    <void *> conn.indices)

    def __str__(self):
        return 'CConnectivity: num: %d, n_incident %d' \
               % (self.num, self.n_incident)

    def cprint(self):
        conn_print(self.conn, stdout)

cdef _create_cconn(MeshConnectivity *pconn, num, n_incident, what):
    if conn_alloc(pconn, num, n_incident):
        raise MemoryError('cannot allocate %s connectivity!' % what)
    cconn = CConnectivity(num, n_incident)
    cconn._set_conn(pconn)
    return cconn

cdef class CMesh:

    @classmethod
    def from_data(cls, coors, vertex_groups, conns, mat_ids, descs,
                  copy_coors=True):
        """
        Fill CMesh data using Python data.
        """
        cdef uint32 tdim
        cdef float64[:, ::1] _coors
        cdef uint32[::1] _cell_types
        cdef MeshConnectivity *pconn

        self = CMesh()

        n_e_ps = np.array([conn.shape[1] for conn in conns])
        n_els = np.array([conn.shape[0] for conn in conns])
        n_el = np.sum(n_els)

        # Max. topological dimension of cells.
        tdim = 0
        for desc in descs:
            tdim = max(tdim, int(desc[0]))

        # Geometry coordinates.
        self.n_coor, self.dim = coors.shape
        if (self.dim < 1) or (self.dim > 3):
            raise ValueError('CMesh geometry dimension must be 1, 2 or 3! (%d)'
                             % self.dim)
        _coors = self.coors = coors.copy() if copy_coors else coors

        mesh_set_coors(self.mesh, &_coors[0, 0], self.n_coor, self.dim, tdim)

        self.vertex_groups = vertex_groups

        # Cell-vertex (D -> 0) connectivity.
        self.n_el = n_el
        self.tdim = tdim
        self.mesh.topology.num[self.tdim] = self.n_el

        _cell_types = self.cell_types = np.empty(self.n_el, dtype=np.uint32)
        self.mesh.topology.cell_types = &_cell_types[0]

        # Length of connectivity.
        n_incident = (n_e_ps * n_els).sum()

        ii = self._get_conn_indx(self.tdim, 0)
        cconn = _create_cconn(self.mesh.topology.conn[ii],
                              self.n_el, n_incident, 'D -> 0')

        self.cell_groups = mat_ids

        indices = []
        offsets = []
        ict = 0
        for ig, conn in enumerate(conns):
            n_el, n_ep = conn.shape

            off = np.empty(n_el, dtype=np.uint32)
            off.fill(n_ep)
            offsets.append(off)
            indices.append(conn.ravel())

            if descs[ig] in self.key_to_index:
                self.cell_types[ict:ict+n_el] = self.key_to_index[descs[ig]]

            else:
                self.cell_types[ict:ict+n_el] = 6 # Higher order mesh.

            ict += n_el

        indices = np.concatenate(indices)
        offsets = np.concatenate(offsets)

        cconn.indices[:] = indices
        cconn.offsets[0] = 0
        cconn.offsets[1:] = np.cumsum(offsets)

        self.conns = [None] * (self.mesh.topology.max_dim + 1)**2
        self.conns[ii] = cconn

        self._update_num()

        return self

    def __cinit__(self):
        mesh_init(self.mesh)
        self.num = np.zeros(4, dtype=np.uint32)

        self.key_to_index = {
            '1_2' : 0,
            '2_3' : 1,
            '2_4' : 2,
            '3_4' : 3,
            '3_8' : 4,
            '3_6' : 5,
        }

    def __dealloc__(self):
        mesh_free(self.mesh)

    def create_new(self,
                   uint32[::1] entities=None,
                   int32 dent=0, localize=False):
        """
        Create a new CMesh instance, with cells corresponding to the given
        `entities` of dimension `dent`.

        Parameters
        ----------
        entities : array, optional
            The selected topological entities of the mesh to be in the new
            mesh. If not given, a copy of the mesh based on the cell-vertex
            connectivity is returned.
        dent : int, optional
            The topological dimension of the entities.
        localize : bool
            If True, strip the vertices not used in the the resulting sub-mesh
            cells and renumber the connectivity.

        Returns
        -------
        cmesh : CMesh
            The new mesh with the cell-vertex connectivity. Other
            connectivities have to be created and local entities need to be set
            manually.
        """
        cdef float64[:, ::1] _coors
        cdef uint32[::1] _ct
        cdef uint32 ii, tdim, n_v, num
        cdef (uint32 *) pct, poffsets
        cdef uint32[::1] indices
        cdef uint32[::1] offsets
        cdef Indices _entities[1]
        cdef MeshConnectivity _incident[1]
        cdef CMesh cmesh

        if entities is None: # Make a copy based cell-vertex connectivity.
            cmesh = CMesh()

            cmesh.n_coor = self.n_coor
            cmesh.dim = self.dim
            cmesh.n_el = self.n_el
            cmesh.tdim = self.tdim

            _coors = cmesh.coors = self.coors.copy()
            mesh_set_coors(cmesh.mesh, &_coors[0, 0],
                           cmesh.n_coor, cmesh.dim, cmesh.tdim)

            cmesh.vertex_groups = self.vertex_groups.copy()
            cmesh.cell_groups = self.cell_groups.copy()

            cmesh.mesh.topology.num[cmesh.tdim] = cmesh.n_el

            _ct = cmesh.cell_types = self.cell_types.copy()
            cmesh.mesh.topology.cell_types = &_ct[0]

            ii = cmesh._get_conn_indx(cmesh.tdim, 0)
            scconn = self.conns[ii]
            cconn = _create_cconn(cmesh.mesh.topology.conn[ii],
                                  cmesh.n_el, scconn.n_incident, 'D -> 0')
            cconn.indices[:] = scconn.indices
            cconn.offsets[:] = scconn.offsets

            cmesh.conns = [None] * (cmesh.mesh.topology.max_dim + 1)**2
            cmesh.conns[ii] = cconn

        else: # Create a submesh based on entities as cells.
            if dent < 1:
                raise ValueError('dimension of entities must be >= 1! (%d)'
                                 % dent)

            cmesh = CMesh()

            cmesh.dim = self.dim
            cmesh.n_el = entities.shape[0]
            cmesh.tdim = dent

            _entities.num = entities.shape[0]
            _entities.indices = &entities[0]

            cmesh.mesh.topology.num[cmesh.tdim] = cmesh.n_el

            num = mesh_count_incident(self.mesh, 0, _entities, dent)

            indices = np.empty(num, dtype=np.uint32)
            offsets = np.empty(_entities.num + 1, dtype=np.uint32)
            _incident.num = _entities.num
            _incident.n_incident = num
            _incident.indices = &indices[0]
            _incident.offsets = &offsets[0]
            mesh_get_incident(self.mesh, _incident, 0, _entities, dent)

            ii0 = (cmesh.tdim + 1) * cmesh.tdim + 0
            cconn = _create_cconn(cmesh.mesh.topology.conn[ii0],
                                  cmesh.n_el, num,'%d -> 0' % dent)
            cconn.offsets[:] = offsets

            if cmesh.tdim == self.tdim: # Cell-based submesh.
                cmesh.cell_groups = self.cell_groups[entities].copy()

                _ct = cmesh.cell_types = self.cell_types[entities].copy()
                cmesh.mesh.topology.cell_types = &_ct[0]

            else: # Face- or edge- based submesh.
                # Cell groups are not unique -> not preserved at all!
                cmesh.cell_groups = np.zeros(cmesh.n_el, dtype=np.int32)

                _ct = cmesh.cell_types = np.empty(cmesh.n_el, dtype=np.uint32)
                pct = cmesh.mesh.topology.cell_types = &_ct[0]
                poffsets = &_incident.offsets[0]
                for ii in range(<uint32>cmesh.n_el):
                    n_v = poffsets[ii + 1] - poffsets[ii]
                    if n_v == 2:
                        pct[ii] = 0
                    elif n_v == 3:
                        pct[ii] = 1
                    elif n_v == 4:
                        pct[ii] = 2

            if localize:
                vertices = np.unique(indices)
                n_new_vertex = vertices.shape[0]

                cmesh.vertex_groups = self.vertex_groups[vertices].copy()

                remap = np.empty(vertices[-1] + 1, dtype=np.uint32)
                remap.fill(n_new_vertex)
                remap[vertices] = np.arange(n_new_vertex, dtype=np.uint32)

                cconn.indices[:] = remap[indices]

                cmesh.n_coor = n_new_vertex

                _coors = cmesh.coors = self.coors[vertices].copy()
                mesh_set_coors(cmesh.mesh, &_coors[0, 0],
                               cmesh.n_coor, cmesh.dim, cmesh.tdim)

            else:
                cmesh.vertex_groups = self.vertex_groups.copy()

                cconn.indices[:] = indices

                cmesh.n_coor = self.n_coor

                _coors = cmesh.coors = self.coors.copy()
                mesh_set_coors(cmesh.mesh, &_coors[0, 0],
                               cmesh.n_coor, cmesh.dim, cmesh.tdim)

            ii = cmesh._get_conn_indx(cmesh.tdim, 0)
            if ii != ii0:
                raise AssertionError('wrong connectivity index!(%d == %d)'
                                     % (ii, ii0))
            cmesh.conns = [None] * (cmesh.mesh.topology.max_dim + 1)**2
            cmesh.conns[ii] = cconn

        cmesh._update_num()

        return cmesh

    def set_local_entities(self, gels):
        cdef (MeshConnectivity *) pedges, pfaces

        self.mesh.entities.num = len(self.key_to_index)

        self.entities = {}

        for key, gel in gels.iteritems():
            try:
                ii = self.key_to_index[key]

            except KeyError:
                continue

            # Reference element edges.
            if gel.n_edge > 0:
                n_incident = gel.n_edge * gel.edges.shape[1]
                cedges = _create_cconn(self.mesh.entities.edges[ii],
                                       gel.n_edge, n_incident, 'local edge')

                cedges.indices[:] = gel.edges.ravel()
                cedges.offsets[0] = 0
                nums = np.empty(gel.n_edge, dtype=np.uint32)
                nums.fill(gel.edges.shape[1])
                cedges.offsets[1:] = np.cumsum(nums)

            else:
                cedges = None

            # Reference element faces.
            if gel.n_face > 0:
                n_incident = gel.n_face * gel.faces.shape[1]
                cfaces = _create_cconn(self.mesh.entities.faces[ii],
                                       gel.n_face, n_incident, 'local face')

                cfaces.indices[:] = gel.faces.ravel()
                cfaces.offsets[0] = 0
                nums = np.empty(gel.n_face, dtype=np.uint32)
                nums.fill(gel.faces.shape[1])
                cfaces.offsets[1:] = np.cumsum(nums)

            else:
                cfaces = None

            self.entities[key] = (cedges, cfaces)

    def get_local_entities(self, key):
        return self.entities[key]

    def setup_entities(self):
        """
        Set up mesh edge (2D and 3D) and face connectivities (3D only) as well
        as their orientations.
        """
        cdef np.npy_intp shape[1]

        if not self.entities:
            msg = 'CMesh.setup_entities() must be called after'\
                  ' CMesh.set_local_entities()!'
            raise ValueError(msg)

        if self.tdim == 1:
            pass

        else:
            self.setup_connectivity(1, 0)

            ii = self._get_conn_indx(self.mesh.topology.max_dim, 1)
            shape[0] = <np.npy_intp> self.conns[ii].n_incident
            ptr = self.mesh.topology.edge_oris
            self.edge_oris = np.PyArray_SimpleNewFromData(1, shape,
                                                          np.NPY_UINT32,
                                                          <void *> ptr)

            if self.tdim == 3:
                self.setup_connectivity(2, 0)

                ii = self._get_conn_indx(self.mesh.topology.max_dim, 2)
                shape[0] = <np.npy_intp> self.conns[ii].n_incident
                ptr = self.mesh.topology.face_oris
                self.face_oris = np.PyArray_SimpleNewFromData(1, shape,
                                                              np.NPY_UINT32,
                                                              <void *> ptr)

                self.facet_oris = self.face_oris

            else:
                self.facet_oris = self.edge_oris

    def setup_connectivity(self, d1, d2):
        cdef MeshConnectivity *pconn

        ii = self._get_conn_indx(d1, d2)
        self.conns[ii] = None

        # This call can create several (intermediate) connectivities.
        mesh_setup_connectivity(self.mesh, d1, d2)

        self._update_pointers()
        self._update_num()

    def _update_pointers(self):
        cdef MeshConnectivity *pconn
        cdef uint32 ii

        for ii in range((self.mesh.topology.max_dim + 1)**2):
            pconn = self.mesh.topology.conn[ii]
            if (pconn.num > 0) and ((self.conns[ii] is None)
                                    or (pconn.num != self.conns[ii].num)):
                cconn = CConnectivity(pconn.num, pconn.n_incident)
                cconn._set_conn(pconn)
                self.conns[ii] = cconn

    def _update_num(self):
        self.num[0] = self.n_coor
        self.num[self.mesh.topology.max_dim] = self.n_el

        for idim in range(1, self.mesh.topology.max_dim):
            ii = self._get_conn_indx(idim, 0)
            conn = self.conns[ii]
            if conn is None:
                self.num[idim] = 0
                self.mesh.topology.num[idim] = 0

            else:
                self.num[idim] = conn.num
                self.mesh.topology.num[idim] = conn.num

    def free_connectivity(self, d1, d2):
        ii = self._get_conn_indx(d1, d2)
        if self.conns[ii] is None:
            return

        self.conns[ii] = None

        mesh_free_connectivity(self.mesh, d1, d2)
        self._update_num()

        max_dim = self.mesh.topology.max_dim
        if (max_dim > d1 > 0) and (d2 == 0):
            self.free_connectivity(max_dim, d1)

        if (d1 == max_dim) and (max_dim > d2 > 0):
            self.free_connectivity(d2, 0)

    def _get_conn_indx(self, d1, d2):
        return (self.mesh.topology.max_dim + 1) * d1 + d2

    def get_conn(self, d1, d2):
        ii = self._get_conn_indx(d1, d2)
        return self.conns[ii]

    def get_cell_conn(self):
        return self.get_conn(self.tdim, 0)

    def get_conn_as_graph(self, d1, d2):
        """
        Get d1 -> d2 connectivity as a sparse matrix graph (values = ones).

        For safety, creates a copy of the connectivity arrays. The connectivity
        is created if necessary.
        """
        import scipy.sparse as sps

        self.setup_connectivity(d1, d2)
        conn = self.get_conn(d1, d2)

        graph = sps.csr_matrix((np.ones(conn.indices.shape[0], dtype=bool),
                                np.array(conn.indices, copy=True,
                                         dtype=np.int32),
                                np.array(conn.offsets, copy=True,
                                         dtype=np.int32)))

        return graph

    def __str__(self):
        return 'CMesh: n_coor: %d, dim %d, tdim: %d, n_el %d' \
               % (self.n_coor, self.dim, self.tdim, self.n_el)

    def cprint(self, int32 header_only=1):
        mesh_print(self.mesh, stdout, header_only)

    def get_surface_facets(self):
        """
        Get facets (edges in 2D, faces in 3D) on the mesh surface.
        """
        self.setup_connectivity(self.tdim - 1, self.tdim)
        conn = self.get_conn(self.tdim - 1, self.tdim)

        ii = np.where(np.diff(conn.offsets) == 1)[0]

        return ii.astype(np.uint32)

    def get_incident(self, int32 dim,
                     uint32[::1] entities not None,
                     int32 dent, ret_offsets=False):
        """
        Get non-unique entities `indices` of dimension `dim` that are contained
        in entities of dimension `dent` listed in `entities`. As each of
        entities can be in several entities of dimension `dent`, `offsets`
        array is returned optionally.
        """
        cdef Indices _entities[1]
        cdef MeshConnectivity _incident[1]
        cdef np.ndarray[uint32, mode='c', ndim=1] indices
        cdef np.ndarray[uint32, mode='c', ndim=1] offsets
        cdef uint32 num

        if not entities.shape[0] > 0:
            empty = np.empty(0, dtype=np.uint32)
            if ret_offsets:
                return empty, empty

            else:
                return empty

        _entities.num = entities.shape[0]
        _entities.indices = &entities[0]

        num = mesh_count_incident(self.mesh, dim, _entities, dent)

        indices = np.empty(num, dtype=np.uint32)
        offsets = np.empty(_entities.num + 1, dtype=np.uint32)

        if num > 0:
            _incident.num = _entities.num
            _incident.n_incident = num
            _incident.indices = &indices[0]
            _incident.offsets = &offsets[0]
            mesh_get_incident(self.mesh, _incident, dim, _entities, dent)

        else:
            offsets[:] = 0

        if ret_offsets:
            return indices, offsets

        else:
            return indices

    def get_local_ids(self,
                      uint32[::1] entities not None,
                      int32 dent,
                      uint32[::1] incident not None,
                      uint32[::1] offsets not None,
                      int32 dim):
        """
        Get local ids of `entities` of dimension `dent` in non-unique entities
        `incident` of dimension `dim` (with given `offsets` per `entities`)
        incident to `entities`, see `mesh_get_incident()`.

        The function searches `entities` in `incident` -> `entities`
        connectivity for each non-unique entity in `incident`.
        """
        cdef Indices[1] _entities, _local_ids
        cdef MeshConnectivity _incident[1]
        cdef np.ndarray[uint32, mode='c', ndim=1] out

        if not entities.shape[0] > 0:
            return np.empty(0, dtype=np.uint32)

        _entities.num = entities.shape[0]
        _entities.indices = &entities[0]

        _incident.num = _entities.num
        _incident.n_incident = incident.shape[0]
        _incident.indices = &incident[0]
        _incident.offsets = &offsets[0]

        out = np.empty(_incident.n_incident, dtype=np.uint32)
        _local_ids.num = _incident.n_incident
        _local_ids.indices = &out[0]
        mesh_get_local_ids(self.mesh, _local_ids, _entities, dent, _incident, dim)

        return out

    def get_complete(self, int32 dim,
                     uint32[::1] entities not None,
                     int32 dent):
        """
        Get entities of dimension `dim` that are completely given by entities
        of dimension `dent` listed in `entities`.
        """
        cdef Mask mask[1]
        cdef Indices _entities[1]
        cdef np.ndarray[uint32, mode='c', ndim=1] out
        cdef uint32 *_out
        cdef uint32 ii, ic

        if not entities.shape[0] > 0:
            return np.empty(0, dtype=np.uint32)

        _entities.num = entities.shape[0]
        _entities.indices = &entities[0]

        mesh_select_complete(self.mesh, mask, dim, _entities, dent)

        out = np.empty(mask.n_true, dtype=np.uint32)

        if mask.n_true > 0:
            _out = &out[0]

            ic = 0
            for ii in range(mask.num):
                if mask.mask[ii]:
                    _out[ic] = ii
                    ic += 1

        pyfree(mask.mask)

        return out

    def get_orientations(self, int32 dim, codim=None):
        """
        Get orientations of entities of dimension `dim`. Alternatively,
        co-dimension can be specified using `codim` argument.
        """
        if codim is not None:
            dim = self.tdim - codim

        if dim == 1:
            return self.edge_oris

        elif dim == 2:
            return self.face_oris

        else:
            raise ValueError('only edges or faces have orientations! (%d)'
                             % dim)

    def get_centroids(self, dim):
        """
        Return the coordinates of centroids of mesh entities with dimension
        `dim`.
        """
        cdef np.ndarray[float64, mode='c', ndim=2] out

        if dim == 0:
            return self.coors

        else:
            out = np.empty((self.mesh.topology.num[dim], self.dim),
                           dtype=np.float64)
            mesh_get_centroids(self.mesh, &out[0, 0], dim)

        return out

    def get_volumes(self, dim):
        """
        Return the volumes of mesh entities with dimension `dim` > 0.
        """
        cdef np.ndarray[float64, mode='c', ndim=1] out

        if dim == 0:
            raise ValueError('vertices have no volume!')

        else:
            out = np.empty((self.mesh.topology.num[dim],),
                           dtype=np.float64)
            mesh_get_volumes(self.mesh, &out[0], dim)

        return out

    def get_facet_normals(self, int32 which=-1):
        """
        Return the normals of facets for each mesh cell. The normals can be
        accessed using the cell-facet connectivity.

        If `which` is -1, two normals of each quadrilateral face are
        averaged. If it is 0 or 1, the corresponding normal is used.
        """
        cdef np.ndarray[float64, mode='c', ndim=2] out

        td = self.tdim
        self.setup_connectivity(td, td - 1)
        ccf = self.get_conn(td, td - 1)

        out = np.empty((ccf.n_incident, self.dim), dtype=np.float64)
        if self.dim == 1:
            # fix for 1D normals, relies on nice ordering of 1D mesh
            # i.e. cell facets are uniformly indexed from left to right
            # globally and locally
            out = np.tile(np.array([-1, 1], dtype=np.float64),
                          ccf.n_incident//2)[:, None]
        else:
            mesh_get_facet_normals(self.mesh, &out[0, 0], which)

        return out

def cmem_statistics():
    mem_statistics(0, '', '', '')

def get_cmem_usage():
    cur_usage = mem_get_cur_usage()
    max_usage = mem_get_max_usage()
    n_frags =  mem_get_n_frags()

    return cur_usage, max_usage, n_frags

## Utils. ##
cdef extern from 'meshutils.h':
    int32 c_orient_elements \
          'orient_elements'(int32 *flag, int32 flag_n_row,
                            Mesh *mesh, Indices *cells, int32 dcells,
                            int32 *v_roots, int32 v_roots_n_row,
                            int32 *v_vecs,
                            int32 v_vecs_n_row, int32 v_vecs_n_col,
                            int32 *swap_from,
                            int32 swap_from_n_row, int32 swap_from_n_col,
                            int32 *swap_to,
                            int32 swap_to_n_row, int32 swap_to_n_col)

@cython.boundscheck(False)
def orient_elements(int32[::1] flag not None,
                    CMesh cmesh not None,
                    uint32[::1] cells not None,
                    int32 dcells,
                    int32[::1] v_roots not None,
                    int32[:, ::1] v_vecs not None,
                    int32[:, ::1] swap_from not None,
                    int32[:, ::1] swap_to not None):
    """
    Swap element nodes so that its volume is positive.
    """
    cdef Indices _cells[1]

    _cells.num = cells.shape[0]
    _cells.indices = &cells[0]

    return c_orient_elements(&flag[0], flag.shape[0],
                             cmesh.mesh, _cells, dcells,
                             &v_roots[0], v_roots.shape[0],
                             &v_vecs[0, 0], v_vecs.shape[0], v_vecs.shape[1],
                             &swap_from[0, 0],
                             swap_from.shape[0], swap_from.shape[1],
                             &swap_to[0, 0],
                             swap_to.shape[0], swap_to.shape[1])
