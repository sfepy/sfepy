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

    ctypedef enum MappingMode:
        MM_Volume
        MM_Surface
        MM_SurfaceExtra

    ctypedef struct Mapping:
        MappingMode mode
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

    cdef int32 map_print(Mapping *obj, FILE *file, int32 mode)

    cdef int32 map_describe(Mapping *obj,
                            float64 *coorIn, int32 nNod, int32 dim,
                            int32 *conn, int32 nEl, int32 nEP,
                            FMField *bfGR, FMField *ebfGR, FMField *weight)

    cdef int32 map_integrate(Mapping *obj, FMField *out, FMField *in_,
                             int32 mode)

    cdef int32 map_getElementDiameters(Mapping *obj, FMField *out,
                                       int32 *edges,
                                       int32 edges_nRow, int32 edges_nCol,
                                       float64 *coorIn, int32 nNod, int32 dim,
                                       int32 *conn, int32 nEl, int32 nEP,
                                       int32 *elList, int32 elList_nRow,
                                       int32 mode)
    cdef int32 map_evaluateBFBGM(Mapping *obj, FMField *bfBGR,
                                 FMField *ebfBGR,
                                 float64 *coorIn, int32 nNod, int32 dim,
                                 int32 *fis, int32 nFa, int32 nFP,
                                 int32 *conn, int32 nEl, int32 nEP)

cdef class CMapping:
    cdef Mapping geo[1]

    cdef FMField _bf[1], _bfg[1], _det[1], _normal[1], _volume[1]

    cdef public np.ndarray bf
    cdef public np.ndarray bfg
    cdef public np.ndarray det
    cdef public np.ndarray normal
    cdef public np.ndarray volume

    cdef public str mode
    cdef public tuple shape
    cdef public int n_el, n_qp, dim, n_ep

    # Auxiliary attributes to be assigned from Python.
    cdef public object integral
    cdef public object qp
    cdef public object ps
