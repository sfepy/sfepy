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
