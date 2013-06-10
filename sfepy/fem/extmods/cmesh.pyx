# -*- Mode: Python -*-
"""
C Mesh data structures and functions.
"""
cimport cython

import numpy as np
cimport numpy as np

from libc.stdio cimport FILE, stdout

from types cimport uint32, int32, float64, complex128

np.import_array()

cdef extern from 'string.h':
    void *memcpy(void *dest, void *src, size_t n)

cdef extern from 'common.h':
    void *pyalloc(size_t size)
    void pyfree(void *pp)
    void mem_statistics(int lineNo, char *funName,
                        char *fileName, char *dirName)
    size_t mem_get_cur_usage()
    size_t mem_get_max_usage()
    size_t mem_get_n_frags()

cdef extern from 'mesh.h':
    ctypedef struct Indices:
        uint32 *indices
        uint32 num

    ctypedef struct Mask:
        char *mask
        uint32 num
        uint32 n_true

    ctypedef struct MeshGeometry:
        uint32 num
        uint32 dim
        float64 *coors

    ctypedef struct MeshTopology:
        uint32 max_dim
        uint32 num[4]
        uint32 *cell_types
        uint32 *face_oris
        uint32 *edge_oris
        MeshConnectivity *conn[16]

    ctypedef struct MeshConnectivity:
        uint32 num
        uint32 n_incident
        uint32 *indices
        uint32 *offsets
        uint32 offset

    ctypedef struct LocalEntities:
        uint32 num
        MeshConnectivity **edges
        MeshConnectivity **faces

    ctypedef struct Mesh:
        MeshGeometry geometry[1]
        MeshTopology topology[1]
        LocalEntities entities[1]

    cdef int32 mesh_init(Mesh *mesh)
    cdef int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only)

    cdef int32 conn_alloc(MeshConnectivity *conn,
                          uint32 num, uint32 n_incident)
    cdef int32 conn_free(MeshConnectivity *conn)
    cdef int32 conn_print(MeshConnectivity *conn, FILE *file)

    cdef int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim)

    cdef int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2)
    cdef int32 mesh_free_connectivity(Mesh *mesh, int32 d1, int32 d2)

    cdef uint32 mesh_count_incident(Mesh *mesh, int32 dim,
                                    Indices *entities, int32 dent)
    cdef int32 mesh_get_incident(Mesh *mesh, Indices *incident, int32 dim,
                                 Indices *entities, int32 dent)
    cdef int32 mesh_select_complete(Mesh *mesh, Mask *mask, int32 dim,
                                    Indices *entities, int32 dent)

cdef class CConnectivity:
    """
    Notes
    -----

    The memory is allocated/freed in C - this class just wraps NumPy arrays
    around that data without copying.
    """
    cdef MeshConnectivity *conn

    cdef public np.ndarray indices
    cdef public np.ndarray offsets
    cdef public int num, n_incident, offset

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
    cdef Mesh mesh[1]

    cdef readonly np.ndarray coors
    cdef readonly np.ndarray cell_types
    cdef readonly np.ndarray cell_groups # ig for each cell.
    cdef readonly list conns
    cdef readonly dict entities
    cdef readonly int n_coor, dim, n_el
    cdef readonly np.ndarray num # Numbers of topological entities.
    cdef readonly np.ndarray face_oris # Allocated in C.
    cdef readonly np.ndarray edge_oris # Allocated in C.
    cdef readonly np.ndarray facet_oris # face_oris in 3D, edge_oris in 2D

    cdef readonly dict key_to_index

    @classmethod
    def from_mesh(cls, mesh):
        """
        Fill data from a Python mesh.
        """
        cdef np.ndarray[float64, mode='c', ndim=2] _coors
        cdef np.ndarray[uint32, mode='c', ndim=1] _cell_types
        cdef MeshConnectivity *pconn

        self = CMesh()

        # Geometry coordinates.
        self.n_coor, self.dim = mesh.coors.shape
        if (self.dim < 2) or (self.dim > 3):
            raise ValueError('CMesh geometry dimension must be 2 or 3! (%d)'
                             % self.dim)
        _coors = self.coors = mesh.coors.copy()
        mesh_set_coors(self.mesh, &_coors[0, 0], self.n_coor, self.dim)

        # Cell-vertex (D -> 0) connectivity.
        self.n_el = mesh.n_el
        self.mesh.topology.num[self.dim] = self.n_el

        _cell_types = self.cell_types = np.empty(self.n_el, dtype=np.uint32)
        self.mesh.topology.cell_types = &_cell_types[0]

        # Length of connectivity.
        n_incident = (mesh.n_e_ps * mesh.n_els).sum()

        ii = self._get_conn_indx(self.dim, 0)
        cconn = _create_cconn(self.mesh.topology.conn[ii],
                              self.n_el, n_incident, 'D -> 0')

        self.cell_groups = np.empty(self.n_el, dtype=np.uint32)

        indices = []
        offsets = []
        ict = 0
        for ig, conn in enumerate(mesh.conns):
            n_el, n_ep = conn.shape

            off = np.empty(n_el, dtype=np.uint32)
            off.fill(n_ep)
            offsets.append(off)
            indices.append(conn.ravel())

            self.cell_types[ict:ict+n_el] = self.key_to_index[mesh.descs[ig]]
            self.cell_groups[ict:ict+n_el] = ig

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
        }

    def set_local_entities(self, gels):
        cdef MeshConnectivity *pedges, *pfaces

        self.mesh.entities.num = len(self.key_to_index)

        self.entities = {}

        for key, gel in gels.iteritems():
            ii = self.key_to_index[key]

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
        Set up mesh edge and face connectivities (3D only) as well as their
        orientations.
        """
        cdef np.npy_intp shape[1]

        if not self.entities:
            msg = 'CMesh.setup_entities() must be called after'\
                  ' CMesh.set_local_entities()!'
            raise ValueError(msg)

        self.setup_connectivity(1, 0)

        shape[0] = <np.npy_intp> self.num[1]
        ptr = self.mesh.topology.edge_oris
        self.edge_oris = np.PyArray_SimpleNewFromData(1, shape,
                                                      np.NPY_UINT32,
                                                      <void *> ptr)

        if self.dim == 3:
            self.setup_connectivity(2, 0)

            shape[0] = <np.npy_intp> self.num[2]
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
        return self.get_conn(self.dim, 0)

    def get_conn_as_graph(self, d1, d2):
        """
        Get d1 -> d2 connectivity as a sparse matrix graph (values = ones).

        For safety, creates a copy of the connectivity arrays. The connectivity
        is created if necessary.
        """
        import scipy.sparse as sps

        self.setup_connectivity(d1, d2)
        conn = self.get_conn(d1, d2)

        graph = sps.csr_matrix((np.ones(conn.indices.shape[0], dtype=np.bool),
                                np.array(conn.indices, copy=True,
                                         dtype=np.int32),
                                np.array(conn.offsets, copy=True,
                                         dtype=np.int32)))

        return graph

    def __str__(self):
        return 'CMesh: n_coor: %d, dim %d, n_el %d' \
               % (self.n_coor, self.dim, self.n_el)

    def cprint(self, int32 header_only=1):
        mesh_print(self.mesh, stdout, header_only)

    def get_surface_facets(self):
        """
        Get facets (edges in 2D, faces in 3D) on the mesh surface.
        """
        self.setup_connectivity(self.dim - 1, self.dim)
        conn = self.get_conn(self.dim - 1, self.dim)

        ii = np.where(np.diff(conn.offsets) == 1)[0]

        return ii

    def get_incident(self, int32 dim,
                     np.ndarray[uint32, mode='c', ndim=1] entities not None,
                     int32 dent):
        """
        Get non-unique entities of dimension `dim` that are contained in
        entities of dimension `dent` listed in `entities`.
        """
        cdef Indices _entities[1], _incident[1]
        cdef np.ndarray[uint32, mode='c', ndim=1] out
        cdef uint32 *_out
        cdef uint32 num

        _entities.num = entities.shape[0]
        _entities.indices = &entities[0]

        num = mesh_count_incident(self.mesh, dim, _entities, dent)

        out = np.empty(num, dtype=np.uint32)
        _incident.num = num
        _incident.indices = &out[0]
        mesh_get_incident(self.mesh, _incident, dim, _entities, dent)

        return out

    def get_complete(self, int32 dim,
                     np.ndarray[uint32, mode='c', ndim=1] entities not None,
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

    def get_from_cell_group(self, int32 ig, int32 dim,
                            np.ndarray[uint32, mode='c', ndim=1] entities=None):
        """
        Get entities of dimension `dim` that are contained in cells of group
        `ig` and are listed in `entities`. If `entities` is None, all entities
        are used.

        Adapter function to be removed after new assembling is done.
        """
        cdef np.ndarray[uint32, mode='c', ndim=1] cells
        cdef np.ndarray[uint32, mode='c', ndim=1] candidates
        cdef np.ndarray[uint32, mode='c', ndim=1] out

        cells = np.where(self.cell_groups == ig)[0].astype(np.uint32)
        if cells.shape[0] == 0:
            raise ValueError('group %d does not exist!' % ig)

        if dim == self.dim:
            candidates = cells

        else:
            self.setup_connectivity(self.dim, dim)
            candidates = self.get_incident(dim, cells, self.dim)

        if entities is not None:
            out = np.intersect1d(candidates, entities)

        else:
            out = candidates

        return out

    def get_igs(self,
                np.ndarray[uint32, mode='c', ndim=1] entities not None,
                int32 dim):
        """
        Get cell groups of incident to entities of dimension `dim`.

        Adapter function to be removed after new assembling is done.
        """
        cdef np.ndarray[uint32, mode='c', ndim=1] cells

        if dim == self.dim:
            cells = entities

        else:
            self.setup_connectivity(dim, self.dim)
            cells = self.get_incident(self.dim, entities, dim)

        igs = np.unique(self.cell_groups[cells])

        return igs

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
                            int32 *conn, int32 conn_n_row, int32 conn_n_col,
                            float64 *coors,
                            int32 coors_n_row, int32 coors_n_col,
                            int32 *v_roots, int32 v_roots_n_row,
                            int32 *v_vecs,
                            int32 v_vecs_n_row, int32 v_vecs_n_col,
                            int32 *swap_from,
                            int32 swap_from_n_row, int32 swap_from_n_col,
                            int32 *swap_to,
                            int32 swap_to_n_row, int32 swap_to_n_col)

    int32 mesh_graph(int32 *p_nnz, int32 **p_prow, int32 **p_icol,
                     int32 n_row, int32 n_col, int32 n_gr, int32 *n_el,
                     int32 *n_epr, int32 **conn_r,
                     int32 *n_epc, int32 **conn_c)

    int32 c_graph_components \
          'graph_components'(int32 *p_n_comp,
                             int32 *flag, int32 flag_len,
                             int32 *row, int32 row_len,
                             int32 *col, int32 col_len,
                             int32 *pos, int32 pos_len)

@cython.boundscheck(False)
def orient_elements(np.ndarray[int32, mode='c', ndim=1] flag not None,
                    np.ndarray[int32, mode='c', ndim=2] conn not None,
                    np.ndarray[float64, mode='c', ndim=2] coors not None,
                    np.ndarray[int32, mode='c', ndim=1] v_roots not None,
                    np.ndarray[int32, mode='c', ndim=2] v_vecs not None,
                    np.ndarray[int32, mode='c', ndim=2] swap_from not None,
                    np.ndarray[int32, mode='c', ndim=2] swap_to not None):
    """
    Swap element nodes so that its volume is positive.
    """
    return c_orient_elements(&flag[0], flag.shape[0],
                             &conn[0, 0], conn.shape[0], conn.shape[1],
                             &coors[0, 0], coors.shape[0], coors.shape[1],
                             &v_roots[0], v_roots.shape[0],
                             &v_vecs[0, 0], v_vecs.shape[0], v_vecs.shape[1],
                             &swap_from[0, 0],
                             swap_from.shape[0], swap_from.shape[1],
                             &swap_to[0, 0],
                             swap_to.shape[0], swap_to.shape[1])

@cython.boundscheck(False)
def create_mesh_graph(int n_row, int n_col, int n_gr, rconns, cconns):
    """
    Create sparse (CSR) graph corresponding to given row and column
    connectivities.

    Parameters
    ----------
    n_row : int
        The number of row connectivity nodes.
    n_col : int
        The number of column connectivity nodes.
    n_gr : int
        The number of element groups.
    rconns : list of arrays
        The list of length `n_gr` of row connectivities.
    cconns : list of arrays
        The list of length `n_gr` of column connectivities.

    Returns
    -------
    nnz : int
        The number of graph nonzeros.
    prow : array
        The array of CSR row pointers.
    icol : array
        The array of CSR column indices.
    """
    cdef int ii
    cdef int32 nnz
    cdef int32 *pprow, *picol
    cdef np.ndarray[int32, mode='c', ndim=1] prow
    cdef np.ndarray[int32, mode='c', ndim=1] icol
    cdef np.ndarray[int32, mode='c', ndim=2] rconn
    cdef np.ndarray[int32, mode='c', ndim=2] cconn
    cdef int32 *n_el = <int32 *> pyalloc(n_gr * sizeof(int32))
    cdef int32 *n_epr = <int32 *> pyalloc(n_gr * sizeof(int32))
    cdef int32 *n_epc = <int32 *> pyalloc(n_gr * sizeof(int32))
    cdef int32 **conn_r = <int32 **> pyalloc(n_gr * sizeof(int32 *))
    cdef int32 **conn_c = <int32 **> pyalloc(n_gr * sizeof(int32 *))

    try:
        for ii in range(0, n_gr):
            rconn = rconns[ii]
            cconn = cconns[ii]

            if rconn.shape[0] != cconn.shape[0]:
                msg = 'connectivities must have the same number of elements!' \
                      ' (%d == %d)' % (rconn.shape[0], cconn.shape[0])
                raise ValueError(msg)

            n_el[ii] = rconn.shape[0]
            n_epr[ii] = rconn.shape[1]
            n_epc[ii] = cconn.shape[1]
            conn_r[ii] = &rconn[0, 0]
            conn_c[ii] = &cconn[0, 0]

        mesh_graph(&nnz, &pprow, &picol, n_row, n_col, n_gr, n_el,
                   n_epr, conn_r, n_epc, conn_c)

        prow = np.empty((n_row + 1,), dtype=np.int32)
        icol = np.empty((nnz,), dtype=np.int32)

        memcpy(&prow[0], pprow, (n_row + 1) * sizeof(int32))
        memcpy(&icol[0], picol, nnz * sizeof(int32))

        pyfree(pprow)
        pyfree(picol)

    finally:
        pyfree(n_el)
        pyfree(n_epr)
        pyfree(n_epc)
        pyfree(conn_r)
        pyfree(conn_c)

    return nnz, prow, icol

@cython.boundscheck(False)
def graph_components(int n_nod,
                     np.ndarray[int32, mode='c', ndim=1] row not None,
                     np.ndarray[int32, mode='c', ndim=1] col not None):
    """
    Determine connected compoments of a compressed sparse graph.

    Returns
    -------
    n_comp : int
        The number of components.
    flag : array
        The flag marking for each node its component.
    """
    cdef int32 n_comp
    cdef np.ndarray[int32, mode='c', ndim=1] flag = np.empty((n_nod,),
                                                             dtype=np.int32)
    cdef np.ndarray[int32, mode='c', ndim=1] pos = np.empty((n_nod,),
                                                            dtype=np.int32)

    c_graph_components(&n_comp, &flag[0], flag.shape[0],
                       &row[0], row.shape[0], &col[0], col.shape[0],
                       &pos[0], pos.shape[0])

    return n_comp, flag
