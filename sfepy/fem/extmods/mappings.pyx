# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
cimport cython

cimport numpy as np
import numpy as np

cimport _fmfield as _f

from types cimport int32, float64, complex128

cdef extern from "fmfield.h":
    ctypedef struct FMField:
        int32 nCell
        int32 nLev
        int32 nRow
        int32 nCol
        float64 *val0
        float64 *val
        int32 nAlloc
        int32 cellSize

        int32 offset
        int32 nColFull

        int32 stride

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

cdef class CVolumeMapping:
    cdef VolumeGeometry geo[1]

    cdef _f.FMField _bfg[1], _det[1], _volume[1]

    cdef public np.ndarray bfg
    cdef public np.ndarray det
    cdef public np.ndarray volume

    cdef public tuple shape

    def __cinit__(self, n_el, n_qp, dim, n_ep):
        self.bfg = np.empty((n_el, n_qp, dim, n_ep), dtype=np.float64)
        _f.array2fmfield4(self._bfg, self.bfg)
        self.geo.bfGM = self._bfg

        self.det = np.empty((n_el, n_qp, 1, 1), dtype=np.float64)
        _f.array2fmfield4(self._det, self.det)
        self.geo.det = self._det

        self.volume = np.empty((n_el, 1, 1, 1), dtype=np.float64)
        _f.array2fmfield4(self._volume, self.volume)
        self.geo.volume = self._volume

        self.shape = (n_el, n_qp, dim, n_ep)

        self.geo.nEl = n_el
        self.geo.nQP = n_qp
        self.geo.dim = dim
        self.geo.nEP = n_ep
        self.geo.mode = GM_Material

    def __str__(self):
        return 'CVolumeMapping: n_el %d, n_qp %d, dim: %d, n_ep: %d' \
               % self.shape

    def describe(self,
                 np.ndarray[float64, mode='c', ndim=2] coors not None,
                 np.ndarray[int32, mode='c', ndim=2] conn not None,
                 np.ndarray[float64, mode='c', ndim=3] bfgr not None,
                 np.ndarray[float64, mode='c', ndim=3] ebfgr not None,
                 np.ndarray[float64, mode='c', ndim=1] weights not None):
        """
        Describe the element geometry - compute the reference element
        mapping.
        """
        cdef _f.FMField _bfgr[1], _ebfgr[1], _weights[1]
        cdef float64 *_coors = &coors[0, 0]
        cdef int32 n_nod = coors.shape[0]
        cdef int32 dim = coors.shape[1]
        cdef int32 *_conn = &conn[0, 0]
        cdef int32 n_el = conn.shape[0]
        cdef int32 n_ep = conn.shape[1]

        # _f.array2fmfield2(_coors, coors)
        _f.array2fmfield3(_bfgr, bfgr)
        _f.array2fmfield3(_ebfgr, ebfgr)
        _f.array2fmfield1(_weights, weights)

        vg_describe(self.geo, _coors, n_nod, dim, _conn, n_el, n_ep,
                    _bfgr, _ebfgr, _weights)
