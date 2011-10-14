# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
cimport cython

cimport numpy as np
import numpy as np

from _fmfield cimport (FMField, array2fmfield4, array2fmfield3,
                       array2fmfield2, array2fmfield1)

from types cimport int32, float64, complex128

cdef extern from 'common.h':
    cdef void errclear()

cdef extern from 'geometry.h':

    ctypedef enum GeometryMode:
        GM_Material
        GM_Spatial

    ctypedef struct VolumeGeometry:
        GeometryMode mode
        int32 nEl
        int32 nQP
        int32 dim
        int32 nEP
        FMField *bfGM
        FMField *det # detJMR or detJSR.
        FMField *volume
        float64 totalVolume

    cdef int32 vg_describe(VolumeGeometry *obj,
                           float64 *coorIn, int32 nNod, int32 dim,
                           int32 *conn, int32 nEl, int32 nEP,
                           FMField *bfGR, FMField *ebfGR, FMField *weight)

    cdef int32 vg_integrate(VolumeGeometry *obj, FMField *out, FMField *in_,
                            int32 mode)

    cdef int32 vg_getElementDiameters(VolumeGeometry *obj, FMField *out,
                                      int32 *edges,
                                      int32 edges_nRow, int32 edges_nCol,
                                      float64 *coorIn, int32 nNod, int32 dim,
                                      int32 *conn, int32 nEl, int32 nEP,
                                      int32 *elList, int32 elList_nRow,
                                      int32 mode)

    ctypedef struct SurfaceGeometry:
        GeometryMode mode
        int32 nFa
        int32 nQP
        int32 dim
        int32 nFP
        FMField *normal
        FMField *det # detJMR.
        FMField *bfBGM

        FMField *area
        float64 totalArea

    cdef int32 sg_describe(SurfaceGeometry *obj,
                           float64 *coorIn, int32 nNod, int32 dim,
                           int32 *fconn, int32 nFa, int32 nSP,
                           FMField *bfGR, FMField *weight)

    cdef int32 sg_integrate(SurfaceGeometry *obj, FMField *out, FMField *in_,
                            int32 mode)

    cdef int32 sg_evaluateBFBGM(SurfaceGeometry *obj, FMField *bfBGR,
                                FMField *ebfBGR,
                                float64 *coorIn, int32 nNod, int32 dim,
                                int32 *fis, int32 nFa, int32 nFP,
                                int32 *conn, int32 nEl, int32 nEP)

cdef class CVolumeMapping:
    cdef VolumeGeometry geo[1]

    cdef FMField _bfg[1], _det[1], _volume[1]

    cdef public np.ndarray bfg
    cdef public np.ndarray det
    cdef public np.ndarray volume

    cdef public tuple shape
    cdef public int n_el, n_qp, dim, n_ep

    # Auxiliary attributes to be assigned from Python.
    cdef public object integral
    cdef public object qp
    cdef public object ps
    cdef public object bf

cdef class CSurfaceMapping:
    cdef SurfaceGeometry geo[1]

    cdef FMField _normal[1], _det[1], _area[1], _bfbg[1]

    cdef public np.ndarray normal
    cdef public np.ndarray det
    cdef public np.ndarray area
    cdef public np.ndarray bfbg

    cdef public tuple shape
    cdef public int n_fa, n_qp, dim, n_fp

    # Auxiliary attributes to be assigned from Python.
    cdef public object integral
    cdef public object qp
    cdef public object ps
    cdef public object bf
