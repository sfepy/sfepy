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

cdef extern from 'mesh.h':
    ctypedef struct MeshGeometry:
        uint32 num
        uint32 dim
        float64 *coors

    ctypedef struct MeshTopology:
        uint32 max_dim
        uint32 num[16]
        MeshConnectivity *conn[16]

    ctypedef struct MeshConnectivity:
        uint32 num
        uint32 n_incident
        uint32 *indices
        uint32 *offsets
        uint32 offset

    ctypedef struct Mesh:
        MeshGeometry geometry[1]
        MeshTopology topology[1]

    cdef int32 mesh_init(Mesh *mesh)
    cdef int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only)

    cdef int32 conn_alloc(MeshConnectivity *conn,
                          uint32 num, uint32 n_incident)
    cdef int32 conn_free(MeshConnectivity *conn)
    cdef int32 conn_print(MeshConnectivity *conn, FILE *file)

    cdef int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim)

    cdef int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2)
    cdef int32 mesh_free_connectivity(Mesh *mesh, int32 d1, int32 d2)

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

cdef class CMesh:
    cdef Mesh mesh[1]

    cdef readonly np.ndarray coors
    cdef readonly list conns
    cdef readonly int n_coor, dim, n_el

    @classmethod
    def from_mesh(cls, mesh):
        """
        Fill data from a Python mesh.
        """
        cdef np.ndarray[float64, mode='c', ndim=2] _coors
        cdef MeshConnectivity *pconn

        self = CMesh()

        # Geometry coordinates.
        self.n_coor, self.dim = mesh.coors.shape
        _coors = self.coors = mesh.coors.copy()
        mesh_set_coors(self.mesh, &_coors[0, 0], self.n_coor, self.dim)

        # Cell-vertex (D -> 0) connectivity.
        self.n_el = mesh.n_el
        self.mesh.topology.num[self.dim] = self.n_el

        # Length of connectivity.
        n_incident = (mesh.n_e_ps * mesh.n_els).sum()

        ii = self._get_conn_indx(self.dim, 0)
        pconn = self.mesh.topology.conn[ii]
        if conn_alloc(pconn, self.n_el, n_incident):
            raise MemoryError('cannot allocate D -> 0 connectivity!')
        cconn = CConnectivity(self.n_el, n_incident)
        cconn._set_conn(pconn)

        indices = []
        offsets = []
        for ig, conn in enumerate(mesh.conns):
            n_el, n_ep = conn.shape

            off = np.empty(n_el, dtype=np.uint32)
            off.fill(n_ep)
            offsets.append(off)
            indices.append(conn.ravel())

        indices = np.concatenate(indices)
        offsets = np.concatenate(offsets)

        cconn.indices[:] = indices
        cconn.offsets[0] = 0
        cconn.offsets[1:] = np.cumsum(offsets)

        self.conns = [None] * (self.mesh.topology.max_dim + 1)**2
        self.conns[ii] = cconn

        return self

    def __cinit__(self):
        mesh_init(self.mesh)

    def setup_connectivity(self, d1, d2):
        cdef MeshConnectivity *pconn

        ii = self._get_conn_indx(d1, d2)
        self.conns[ii] = None

        mesh_setup_connectivity(self.mesh, d1, d2)

        pconn = self.mesh.topology.conn[ii]
        cconn = CConnectivity(pconn.num, pconn.n_incident)
        cconn._set_conn(pconn)

        self.conns[ii] = cconn

    def free_connectivity(self, d1, d2):
        cdef MeshConnectivity *pconn

        ii = self._get_conn_indx(d1, d2)
        self.conns[ii] = None

        mesh_free_connectivity(self.mesh, d1, d2)

    def _get_conn_indx(self, d1, d2):
        return (self.mesh.topology.max_dim + 1) * d1 + d2

    def get_conn(self, d1, d2):
        ii = self._get_conn_indx(d1, d2)
        return self.conns[ii]

    def get_cell_conn(self):
        return self.get_conn(self.dim, 0)

    def __str__(self):
        return 'CMesh: n_coor: %d, dim %d, n_el %d' \
               % (self.n_coor, self.dim, self.n_el)

    def cprint(self, int32 header_only=1):
        mesh_print(self.mesh, stdout, header_only)

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
