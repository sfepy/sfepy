# -*- Mode: Python -*-
cimport cython

@cython.boundscheck(False)
cdef inline void array2fmfield4(FMField *out,
                                np.ndarray[float64, mode='c', ndim=4] arr):
    cdef int32 n_cell, n_lev, n_row, n_col
    cdef int32 ii, stride

    sh = arr.shape
    n_cell, n_lev, n_row, n_col = sh[0], sh[1], sh[2], sh[3]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0, 0, 0, 0])
    out.stride = 1

@cython.boundscheck(False)
cdef inline void array2fmfield3(FMField *out,
                                np.ndarray[float64, mode='c', ndim=3] arr):
    cdef int32 n_cell = 1, n_lev, n_row, n_col
    cdef int32 ii, stride

    sh = arr.shape
    n_lev, n_row, n_col = sh[0], sh[1], sh[2]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0, 0, 0])
    out.stride = 1

@cython.boundscheck(False)
cdef inline void array2fmfield2(FMField *out,
                                np.ndarray[float64, mode='c', ndim=2] arr):
    cdef int32 n_cell = 1, n_lev = 1, n_row, n_col
    cdef int32 ii, stride

    sh = arr.shape
    n_row, n_col = sh[0], sh[1]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0, 0])
    out.stride = 1

@cython.boundscheck(False)
cdef inline void array2fmfield1(FMField *out,
                                np.ndarray[float64, mode='c', ndim=1] arr):
    cdef int32 n_cell = 1, n_lev = 1, n_row = 1, n_col
    cdef int32 ii, stride

    sh = arr.shape
    n_col = sh[0]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0])
    out.stride = 1
