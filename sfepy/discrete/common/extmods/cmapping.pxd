# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
cimport cython

from _fmfield cimport (FMField, array2fmfield4)

from types cimport int32, float64

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
