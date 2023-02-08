# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
cimport cython

cimport numpy as np
import numpy as np

from libc.stdio cimport FILE, stdout

from _fmfield cimport (FMField, array2fmfield4, array2fmfield3,
                       array2fmfield2, array2fmfield1)

from types cimport int32, float64, complex128

cdef extern from 'common.h':
    cdef void errclear()

cdef extern from 'refmaps.h':

    ctypedef struct Mapping:
        int32 nEl
        int32 nQP
        int32 dim
        int32 nEP
        FMField *bf
        FMField *bfGM # Volume or SurfaceExtra only.
        FMField *det # detJMR or detJSR.
        FMField *normal # Surface only.
        FMField *volume
        float64 totalVolume

cdef class CMapping:
    cdef Mapping geo[1]
    cdef FMField[1] _bf, _bfg, _det, _normal, _volume

    cdef public np.ndarray bf
    cdef public np.ndarray bfg
    cdef public np.ndarray det
    cdef public np.ndarray normal
    cdef public np.ndarray volume

    cdef public int n_el, n_qp, dim, n_ep

    # cdef public object integral
    # cdef public object qp
    # cdef public object ps
    # cdef public object mtx_t
