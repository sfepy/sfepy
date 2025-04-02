# -*- Mode: Python -*-
# cython: language_level=3
cimport cython

@cython.boundscheck(False)
cdef inline int array2fmfield4(FMField *out, float64[:, :, :, ::1] arr) \
                               except -1:
    cdef int32 n_cell, n_lev, n_row, n_col
    cdef int32 ii

    sh = arr.shape
    n_cell, n_lev, n_row, n_col = sh[0], sh[1], sh[2], sh[3]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0, 0, 0, 0])

@cython.boundscheck(False)
cdef inline int array2fmfield3(FMField *out, float64[:, :, ::1] arr) except -1:
    cdef int32 n_cell = 1, n_lev, n_row, n_col
    cdef int32 ii

    sh = arr.shape
    n_lev, n_row, n_col = sh[0], sh[1], sh[2]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0, 0, 0])

@cython.boundscheck(False)
cdef inline int array2fmfield2(FMField *out, float64[:, ::1] arr) except -1:
    cdef int32 n_cell = 1, n_lev = 1, n_row, n_col
    cdef int32 ii

    sh = arr.shape
    n_row, n_col = sh[0], sh[1]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0, 0])

@cython.boundscheck(False)
cdef inline int array2fmfield1(FMField *out, float64[::1] arr) except -1:
    cdef int32 n_cell = 1, n_lev = 1, n_row = 1, n_col
    cdef int32 ii

    sh = arr.shape
    n_col = sh[0]

    out.nAlloc = -1
    fmf_pretend(out, n_cell, n_lev, n_row, n_col, &arr[0])

@cython.boundscheck(False)
cdef inline int array2pint2(int32 **out, int32 *n_row, int32 *n_col,
                            int32[:, ::1] arr) except -1:
    out[0] = &arr[0, 0]
    n_row[0] = arr.shape[0]
    n_col[0] = arr.shape[1]

@cython.boundscheck(False)
cdef inline int array2pint1(int32 **out, int32 *n_row,
                            int32[::1] arr) except -1:
    out[0] = &arr[0]
    n_row[0] = arr.shape[0]

@cython.boundscheck(False)
cdef inline int array2puint2(uint32 **out, uint32 *n_row, uint32 *n_col,
                             uint32[:, ::1] arr) except -1:
    out[0] = &arr[0, 0]
    n_row[0] = arr.shape[0]
    n_col[0] = arr.shape[1]

@cython.boundscheck(False)
cdef inline int array2puint1(uint32 **out, uint32 *n_row,
                             uint32[::1] arr) except -1:
    out[0] = &arr[0]
    n_row[0] = arr.shape[0]
