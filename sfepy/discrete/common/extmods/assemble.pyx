# -*- Mode: Python -*-
# cython: language_level=3
"""
Low level finite element assembling functions.
"""
cimport cython

import numpy as np
cimport numpy as np

from sfepy.discrete.common.extmods.types cimport int32, uint32, float64, complex128

@cython.boundscheck(False)
def assemble_vector(float64[::1] vec not None,
                    float64[:, :, :, ::1] vec_in_els not None,
                    int32[::1] iels not None,
                    float64 sign,
                    int32[:, ::1] conn not None):
    cdef int32 ii, iel, ir, irg
    cdef int32 num = iels.shape[0]
    cdef int32 n_ep = conn.shape[1]
    # Allow both row or column vectors.
    cdef int32 cell_size = vec_in_els.shape[2] * vec_in_els.shape[3]
    cdef (int32 *) pconn0, pconn
    cdef int32 *piels = &iels[0]
    cdef float64 *val = &vec[0]
    cdef (float64 *) vec_in_el0, vec_in_el

    assert num == vec_in_els.shape[0]

    pconn0 = &conn[0, 0]
    vec_in_el0 = &vec_in_els[0, 0, 0, 0]

    for ii in range(0, num):
        iel = piels[ii]

        pconn = pconn0 + iel * n_ep
        vec_in_el = vec_in_el0 + ii * cell_size

        for ir in range(0, n_ep):
            irg = pconn[ir]
            if irg < 0: continue

            val[irg] += sign * vec_in_el[ir]

@cython.boundscheck(False)
def assemble_vector_complex(complex128[::1] vec not None,
                            complex128[:, :, :, ::1] vec_in_els not None,
                            int32[::1] iels not None,
                            complex128 sign,
                            int32[:, ::1] conn not None):
    cdef int32 ii, iel, ir, irg
    cdef int32 num = iels.shape[0]
    cdef int32 n_ep = conn.shape[1]
    # Allow both row or column vectors.
    cdef int32 cell_size = vec_in_els.shape[2] * vec_in_els.shape[3]
    cdef (int32 *) pconn0, pconn
    cdef int32 *piels = &iels[0]
    cdef complex128 *val = &vec[0]
    cdef (complex128 *) vec_in_el0, vec_in_el

    assert num == vec_in_els.shape[0]

    pconn0 = &conn[0, 0]
    vec_in_el0 = &vec_in_els[0, 0, 0, 0]

    for ii in range(0, num):
        iel = piels[ii]

        pconn = pconn0 + iel * n_ep
        vec_in_el = vec_in_el0 + ii * cell_size

        for ir in range(0, n_ep):
            irg = pconn[ir]
            if irg < 0: continue

            val[irg] += sign * vec_in_el[ir]

@cython.boundscheck(False)
def assemble_matrix(float64[::1] mtx not None,
                    int32[::1] prows not None,
                    int32[::1] cols not None,
                    float64[:, :, :, ::1] mtx_in_els not None,
                    int32[::1] iels not None,
                    float64 sign,
                    int32[:, ::1] row_conn not None,
                    int32[:, ::1] col_conn not None):
    cdef int32 ii, iel, ir, ic, irg, icg, ik, iloc
    cdef int32 num = iels.shape[0]
    cdef int32 n_epr = row_conn.shape[1]
    cdef int32 n_epc = col_conn.shape[1]
    cdef int32 cell_size = mtx_in_els.shape[2] * mtx_in_els.shape[3]
    cdef (int32 *) prow_conn0, pcol_conn0, prow_conn, pcol_conn
    cdef int32 *piels = &iels[0]
    cdef int32 *_prows = &prows[0]
    cdef int32 *_cols = &cols[0]
    cdef float64 *val = &mtx[0]
    cdef (float64 *) mtx_in_el0, mtx_in_el

    assert num == mtx_in_els.shape[0]

    prow_conn0 = &row_conn[0, 0]
    pcol_conn0 = &col_conn[0, 0]
    mtx_in_el0 = &mtx_in_els[0, 0, 0, 0]

    for ii in range(0, num):
        iel = piels[ii]

        prow_conn = prow_conn0 + iel * n_epr
        pcol_conn = pcol_conn0 + iel * n_epc
        mtx_in_el = mtx_in_el0 + ii * cell_size

        for ir in range(0, n_epr):
            irg = prow_conn[ir]
            if irg < 0: continue

            for ic in range(0, n_epc):
                icg = pcol_conn[ic]
                if icg < 0: continue

                iloc = n_epc * ir + ic

                for ik in range(_prows[irg], _prows[irg + 1]):
                    if _cols[ik] == icg:
                        val[ik] += sign * mtx_in_el[iloc]
                        break

                else:
                    msg = 'matrix item (%d, %d) does not exist!' % (irg, icg)
                    raise IndexError(msg)

@cython.boundscheck(False)
def assemble_matrix_complex(complex128[::1] mtx not None,
                            int32[::1] prows not None,
                            int32[::1] cols not None,
                            complex128[:, :, :, ::1] mtx_in_els not None,
                            int32[::1] iels not None,
                            complex128 sign,
                            int32[:, ::1]
                            row_conn not None,
                            int32[:, ::1]
                            col_conn not None):
    cdef int32 ii, iel, ir, ic, irg, icg, ik, iloc
    cdef int32 num = iels.shape[0]
    cdef int32 n_epr = row_conn.shape[1]
    cdef int32 n_epc = col_conn.shape[1]
    cdef int32 cell_size = mtx_in_els.shape[2] * mtx_in_els.shape[3]
    cdef (int32 *) prow_conn0, pcol_conn0, prow_conn, pcol_conn
    cdef int32 *piels = &iels[0]
    cdef int32 *_prows = &prows[0]
    cdef int32 *_cols = &cols[0]
    cdef complex128 *val = &mtx[0]
    cdef (complex128 *) mtx_in_el0, mtx_in_el

    assert num == mtx_in_els.shape[0]

    prow_conn0 = &row_conn[0, 0]
    pcol_conn0 = &col_conn[0, 0]
    mtx_in_el0 = &mtx_in_els[0, 0, 0, 0]

    for ii in range(0, num):
        iel = piels[ii]

        prow_conn = prow_conn0 + iel * n_epr
        pcol_conn = pcol_conn0 + iel * n_epc
        mtx_in_el = mtx_in_el0 + ii * cell_size

        for ir in range(0, n_epr):
            irg = prow_conn[ir]
            if irg < 0: continue

            for ic in range(0, n_epc):
                icg = pcol_conn[ic]
                if icg < 0: continue

                iloc = n_epc * ir + ic

                for ik in range(_prows[irg], _prows[irg + 1]):
                    if _cols[ik] == icg:
                        val[ik] += sign * mtx_in_el[iloc]
                        break

                else:
                    msg = 'matrix item (%d, %d) does not exist!' % (irg, icg)
                    raise IndexError(msg)

cdef inline int32 bsearch(int32 *cols,
                          int32 i0,
                          int32 i1,
                          int32 ii) nogil:
    cdef int32 im

    while i0 <= i1:
        im = <int32> ((<uint32> i0 + <uint32> i1) >> 1)

        if cols[im] == ii:
            return im

        elif cols[im] < ii:
            i0 = im + 1

        else:
            i1 = im - 1

    return -1

@cython.boundscheck(False)
def assemble_matrix_b(float64[::1] mtx not None,
                      int32[::1] prows not None,
                      int32[::1] cols not None,
                      float64[:, :, :, ::1] mtx_in_els not None,
                      int32[::1] iels not None,
                      float64 sign,
                      int32[:, ::1] row_conn not None,
                      int32[:, ::1] col_conn not None):
    cdef int32 ii, iel, ir, ic, irg, icg, ik, iloc
    cdef int32 num = iels.shape[0]
    cdef int32 n_epr = row_conn.shape[1]
    cdef int32 n_epc = col_conn.shape[1]
    cdef int32 cell_size = mtx_in_els.shape[2] * mtx_in_els.shape[3]
    cdef (int32 *) prow_conn0, pcol_conn0, prow_conn, pcol_conn
    cdef int32 *piels = &iels[0]
    cdef int32 *_prows = &prows[0]
    cdef int32 *_cols = &cols[0]
    cdef float64 *val = &mtx[0]
    cdef (float64 *) mtx_in_el0, mtx_in_el

    assert num == mtx_in_els.shape[0]

    prow_conn0 = &row_conn[0, 0]
    pcol_conn0 = &col_conn[0, 0]
    mtx_in_el0 = &mtx_in_els[0, 0, 0, 0]

    for ii in range(0, num):
        iel = piels[ii]

        prow_conn = prow_conn0 + iel * n_epr
        pcol_conn = pcol_conn0 + iel * n_epc
        mtx_in_el = mtx_in_el0 + ii * cell_size

        for ir in range(0, n_epr):
            irg = prow_conn[ir]
            if irg < 0: continue

            for ic in range(0, n_epc):
                icg = pcol_conn[ic]
                if icg < 0: continue

                iloc = n_epc * ir + ic

                ik = bsearch(_cols, _prows[irg], _prows[irg + 1], icg)
                if ik >= 0:
                    val[ik] += sign * mtx_in_el[iloc]

                else:
                    msg = 'matrix item (%d, %d) does not exist!' % (irg, icg)
                    raise IndexError(msg)

@cython.boundscheck(False)
def assemble_matrix_complex_b(complex128[::1] mtx not None,
                              int32[::1] prows not None,
                              int32[::1] cols not None,
                              complex128[:, :, :, ::1] mtx_in_els not None,
                              int32[::1] iels not None,
                              complex128 sign,
                              int32[:, ::1]
                              row_conn not None,
                              int32[:, ::1]
                              col_conn not None):
    cdef int32 ii, iel, ir, ic, irg, icg, ik, iloc
    cdef int32 num = iels.shape[0]
    cdef int32 n_epr = row_conn.shape[1]
    cdef int32 n_epc = col_conn.shape[1]
    cdef int32 cell_size = mtx_in_els.shape[2] * mtx_in_els.shape[3]
    cdef (int32 *) prow_conn0, pcol_conn0, prow_conn, pcol_conn
    cdef int32 *piels = &iels[0]
    cdef int32 *_prows = &prows[0]
    cdef int32 *_cols = &cols[0]
    cdef complex128 *val = &mtx[0]
    cdef (complex128 *) mtx_in_el0, mtx_in_el

    assert num == mtx_in_els.shape[0]

    prow_conn0 = &row_conn[0, 0]
    pcol_conn0 = &col_conn[0, 0]
    mtx_in_el0 = &mtx_in_els[0, 0, 0, 0]

    for ii in range(0, num):
        iel = piels[ii]

        prow_conn = prow_conn0 + iel * n_epr
        pcol_conn = pcol_conn0 + iel * n_epc
        mtx_in_el = mtx_in_el0 + ii * cell_size

        for ir in range(0, n_epr):
            irg = prow_conn[ir]
            if irg < 0: continue

            for ic in range(0, n_epc):
                icg = pcol_conn[ic]
                if icg < 0: continue

                iloc = n_epc * ir + ic

                ik = bsearch(_cols, _prows[irg], _prows[irg + 1], icg)
                if ik >= 0:
                    val[ik] += sign * mtx_in_el[iloc]

                else:
                    msg = 'matrix item (%d, %d) does not exist!' % (irg, icg)
                    raise IndexError(msg)
