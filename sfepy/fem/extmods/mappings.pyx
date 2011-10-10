# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
import numpy as np

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

    def __cinit__(self, n_el, n_qp, dim, n_ep):
        self.bfg = np.empty((n_el, n_qp, dim, n_ep), dtype=np.float64)
        array2fmfield4(self._bfg, self.bfg)
        self.geo.bfGM = self._bfg

        self.det = np.empty((n_el, n_qp, 1, 1), dtype=np.float64)
        array2fmfield4(self._det, self.det)
        self.geo.det = self._det

        self.volume = np.empty((n_el, 1, 1, 1), dtype=np.float64)
        array2fmfield4(self._volume, self.volume)
        self.geo.volume = self._volume

        self.shape = (n_el, n_qp, dim, n_ep)

        self.geo.nEl = self.n_el = n_el
        self.geo.nQP = self.n_qp = n_qp
        self.geo.dim = self.dim = dim
        self.geo.nEP = self.n_ep = n_ep
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
        cdef int32 ret = 0
        cdef FMField _bfgr[1], _ebfgr[1], _weights[1]
        cdef float64 *_coors = &coors[0, 0]
        cdef int32 n_nod = coors.shape[0]
        cdef int32 dim = coors.shape[1]
        cdef int32 *_conn = &conn[0, 0]
        cdef int32 n_el = conn.shape[0]
        cdef int32 n_ep = conn.shape[1]

        # array2fmfield2(_coors, coors)
        array2fmfield3(_bfgr, bfgr)
        array2fmfield3(_ebfgr, ebfgr)
        array2fmfield1(_weights, weights)

        ret = vg_describe(self.geo, _coors, n_nod, dim, _conn, n_el, n_ep,
                          _bfgr, _ebfgr, _weights)

        if ret:
            errclear()
            raise ValueError('ccore error (see above)')

    def integrate(self,
                  np.ndarray[float64, mode='c', ndim=4] out not None,
                  np.ndarray[float64, mode='c', ndim=4] arr not None,
                  int32 mode=0):
        """
        Integrate `arr` over the domain of the mapping into `out`.
        """
        cdef int32 ret = 0
        cdef FMField _out[1], _arr[1]

        if not ((out.shape[0] == arr.shape[0])
                and (out.shape[1] == 1)
                and (out.shape[2] == arr.shape[2])
                and (out.shape[3] == arr.shape[3])
                and (out.shape[0] == self.n_el)
                and (arr.shape[1] == self.n_qp)):
            raise ValueError('incompatible shapes! (%s, out: (%d, %d, %d, %d)'
                             ', arr: (%d, %d, %d, %d))'
                             % (self.shape,
                                out.shape[0], out.shape[1],
                                out.shape[2], out.shape[3],
                                arr.shape[0], arr.shape[1],
                                arr.shape[2], arr.shape[3]))

        array2fmfield4(_out, out)
        array2fmfield4(_arr, arr)

        ret = vg_integrate(self.geo, _out, _arr, mode)

        if ret:
            errclear()
            raise ValueError('ccore error (see above)')
