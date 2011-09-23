# -*- Mode: Python -*-
"""
Low level mesh functions employing element connectivity.
"""
cimport cython

import numpy as np
cimport numpy as np

ctypedef np.complex128_t complex128
ctypedef np.float64_t float64
ctypedef np.int32_t int32

cdef extern from 'string.h':
    void *memcpy(void *dest, void *src, size_t n)

cdef extern from 'common.h':
    void *pyalloc(size_t size)
    void pyfree(void *pp)

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
