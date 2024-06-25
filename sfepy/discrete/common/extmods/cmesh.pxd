# -*- Mode: Python -*-
# cython: language_level=3
"""
C Mesh data structures and functions.
"""
cimport numpy as np

from libc.stdio cimport FILE, stdout

from sfepy.discrete.common.extmods.types cimport (uint32, int32, float64,
                                                  complex128)

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
    cdef int32 mesh_free(Mesh *mesh)
    cdef int32 mesh_print(Mesh *mesh, FILE *file, int32 header_only)

    cdef int32 conn_alloc(MeshConnectivity *conn,
                          uint32 num, uint32 n_incident)
    cdef int32 conn_free(MeshConnectivity *conn)
    cdef int32 conn_print(MeshConnectivity *conn, FILE *file)

    cdef int32 mesh_set_coors(Mesh *mesh, float64 *coors, int32 num, int32 dim,
                              int32 tdim)

    cdef int32 mesh_setup_connectivity(Mesh *mesh, int32 d1, int32 d2)
    cdef int32 mesh_free_connectivity(Mesh *mesh, int32 d1, int32 d2)

    cdef uint32 mesh_count_incident(Mesh *mesh, int32 dim,
                                    Indices *entities, int32 dent)
    cdef int32 mesh_get_incident(Mesh *mesh,
                                 MeshConnectivity *incident, int32 dim,
                                 Indices *entities, int32 dent)
    cdef int32 mesh_get_local_ids(Mesh *mesh, Indices *local_ids,
                                  Indices *entities, int32 dent,
                                  MeshConnectivity *incident, int32 dim)
    cdef int32 mesh_select_complete(Mesh *mesh, Mask *mask, int32 dim,
                                    Indices *entities, int32 dent)
    cdef int32 mesh_get_centroids(Mesh *mesh, float64 *ccoors, int32 dim)
    cdef int32 mesh_get_volumes(Mesh *mesh, float64 *volumes, int32 dim)
    cdef int32 mesh_get_facet_normals(Mesh *mesh, float64 *normals,
                                      int32 which)

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

    cdef _set_conn(self, MeshConnectivity *conn)

cdef class CMesh:
    cdef Mesh mesh[1]

    cdef readonly np.ndarray coors
    cdef readonly np.ndarray vertex_groups
    cdef readonly np.ndarray cell_types
    cdef readonly np.ndarray cell_groups # ig for each cell.
    cdef readonly list conns
    cdef readonly dict entities
    cdef readonly int n_coor, dim, n_el, tdim
    cdef readonly np.ndarray num # Numbers of topological entities.
    cdef readonly np.ndarray face_oris # Allocated in C.
    cdef readonly np.ndarray edge_oris # Allocated in C.
    cdef readonly np.ndarray facet_oris # face_oris in 3D, edge_oris in 2D

    cdef readonly dict key_to_index
