# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
import numpy as np

cdef class CVolumeMapping:

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

    def cprint(self, int32 mode=0):
        vg_print(self.geo, stdout, mode)

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

        return ret

    def get_element_diameters(self,
                              np.ndarray[float64, mode='c', ndim=4]
                              out not None,
                              np.ndarray[int32, mode='c', ndim=2]
                              edges not None,
                              np.ndarray[float64, mode='c', ndim=2]
                              coors not None,
                              np.ndarray[int32, mode='c', ndim=2]
                              conn not None,
                              np.ndarray[int32, mode='c', ndim=1]
                              el_list not None,
                              int32 mode):
        """
        Compute diameters of selected elements.
        """
        cdef int32 ret = 0
        cdef FMField _out[1]
        cdef int32 *_edges = &edges[0, 0]
        cdef int32 edges_nr = edges.shape[0]
        cdef int32 edges_nc = edges.shape[1]
        cdef float64 *_coors = &coors[0, 0]
        cdef int32 n_nod = coors.shape[0]
        cdef int32 dim = coors.shape[1]
        cdef int32 *_conn = &conn[0, 0]
        cdef int32 n_el = conn.shape[0]
        cdef int32 n_ep = conn.shape[1]
        cdef int32 *_el_list = &el_list[0]
        cdef int32 n_el2 = el_list.shape[0]

        array2fmfield4(_out, out)

        ret = vg_getElementDiameters(self.geo, _out,
                                     _edges, edges_nr, edges_nc,
                                     _coors, n_nod, dim,
                                     _conn, n_el, n_ep, _el_list, n_el2, mode)

        if ret:
            errclear()
            raise ValueError('ccore error (see above)')

        return ret

cdef class CSurfaceMapping:

    def __cinit__(self, n_fa, n_qp, dim, n_fp):
        self.normal = np.empty((n_fa, n_qp, dim, 1), dtype=np.float64)
        array2fmfield4(self._normal, self.normal)
        self.geo.normal = self._normal

        self.det = np.empty((n_fa, n_qp, 1, 1), dtype=np.float64)
        array2fmfield4(self._det, self.det)
        self.geo.det = self._det

        self.area = np.empty((n_fa, 1, 1, 1), dtype=np.float64)
        array2fmfield4(self._area, self.area)
        self.geo.area = self._area

        self.bfbg = None
        self.geo.bfBGM = NULL

        self.shape = (n_fa, n_qp, dim, n_fp)

        self.geo.nFa = self.n_fa = n_fa
        self.geo.nQP = self.n_qp = n_qp
        self.geo.dim = self.dim = dim
        self.geo.nFP = self.n_fp = n_fp
        self.geo.mode = GM_Material

    def alloc_extra_data(self, int n_ep):
        """
        Allocate boundary base function gradient array.

        This must be called before calling CSurfaceMapping.evaluate_bfbgm()!
        """
        self.bfbg = np.empty((self.n_fa, self.n_qp, self.dim, n_ep),
                             dtype=np.float64)
        array2fmfield4(self._bfbg, self.bfbg)
        self.geo.bfBGM = self._bfbg

    def __str__(self):
        return 'CSurfaceMapping: n_fa %d, n_qp %d, dim: %d, n_fp: %d' \
               % self.shape

    def cprint(self, int32 mode=0):
        sg_print(self.geo, stdout, mode)

    def describe(self,
                 np.ndarray[float64, mode='c', ndim=2] coors not None,
                 np.ndarray[int32, mode='c', ndim=2] conn not None,
                 np.ndarray[float64, mode='c', ndim=3] bfgr not None,
                 np.ndarray[float64, mode='c', ndim=1] weights not None):
        """
        Describe the elememt surface geometry - compute the reference
        element surface mapping.
        """
        cdef int32 ret = 0
        cdef FMField _bfgr[1], _weights[1]
        cdef float64 *_coors = &coors[0, 0]
        cdef int32 n_nod = coors.shape[0]
        cdef int32 dim = coors.shape[1]
        cdef int32 *_conn = &conn[0, 0]
        cdef int32 n_fa = conn.shape[0]
        cdef int32 n_fp = conn.shape[1]

        array2fmfield3(_bfgr, bfgr)
        array2fmfield1(_weights, weights)

        ret = sg_describe(self.geo, _coors, n_nod, dim, _conn, n_fa, n_fp,
                          _bfgr, _weights)

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
                and (out.shape[0] == self.n_fa)
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

        ret = sg_integrate(self.geo, _out, _arr, mode)

        if ret:
            errclear()
            raise ValueError('ccore error (see above)')

        return ret

    def evaluate_bfbgm(self,
                       np.ndarray[float64, mode='c', ndim=4] bfbgr not None,
                       np.ndarray[float64, mode='c', ndim=4] ebfbgr not None,
                       np.ndarray[float64, mode='c', ndim=2] coors not None,
                       np.ndarray[int32, mode='c', ndim=2] fis not None,
                       np.ndarray[int32, mode='c', ndim=2] conn not None):
        """
        Evaluate volume base function gradients in surface quadrature
        points.
        """
        cdef int32 ret = 0
        cdef FMField _bfbgr[1], _ebfbgr[1]
        cdef float64 *_coors = &coors[0, 0]
        cdef int32 n_nod = coors.shape[0]
        cdef int32 dim = coors.shape[1]
        cdef int32 *_fis = &fis[0, 0]
        cdef int32 n_fa = fis.shape[0]
        cdef int32 n_fp = fis.shape[1]
        cdef int32 *_conn = &conn[0, 0]
        cdef int32 n_el = conn.shape[0]
        cdef int32 n_ep = conn.shape[1]

        if self.bfbg is None:
            raise ValueError('CSurfaceMapping.alloc_extra_data() not called!')

        array2fmfield4(_bfbgr, bfbgr)
        array2fmfield4(_ebfbgr, ebfbgr)

        ret = sg_evaluateBFBGM(self.geo, _bfbgr, _ebfbgr, _coors, n_nod, dim,
                               _fis, n_fa, n_fp, _conn, n_el, n_ep)

        if ret:
            errclear()
            raise ValueError('ccore error (see above)')
