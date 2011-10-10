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
