# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""
import numpy as np

cdef class CMapping:

    def __cinit__(self, n_el, n_qp, dim, n_ep, mode='volume', flag=0):
        if flag:
            self.bf = np.empty((n_el, n_qp, 1, n_ep), dtype=np.float64)

        else:
            self.bf = np.empty((1, n_qp, 1, n_ep), dtype=np.float64)

        array2fmfield4(self._bf, self.bf)
        self.geo.bf = self._bf

        self.det = np.empty((n_el, n_qp, 1, 1), dtype=np.float64)
        array2fmfield4(self._det, self.det)
        self.geo.det = self._det

        self.volume = np.empty((n_el, 1, 1, 1), dtype=np.float64)
        array2fmfield4(self._volume, self.volume)
        self.geo.volume = self._volume

        if mode == 'volume':
            self.bfg = np.empty((n_el, n_qp, dim, n_ep), dtype=np.float64)
            array2fmfield4(self._bfg, self.bfg)
            self.geo.bfGM = self._bfg

            self.normal = None
            self.geo.normal = NULL

        else:
            self.bfg = None
            self.geo.bfGM = NULL

            self.normal = np.empty((n_el, n_qp, dim, 1), dtype=np.float64)
            array2fmfield4(self._normal, self.normal)
            self.geo.normal = self._normal

        self.mode = mode
        self.shape = (n_el, n_qp, dim, n_ep)

        self.geo.mode = {'volume' : MM_Volume,
                         'surface' : MM_Surface,
                         'surface_extra' : MM_SurfaceExtra}[mode]
        self.geo.nEl = self.n_el = n_el
        self.geo.nQP = self.n_qp = n_qp
        self.geo.dim = self.dim = dim
        self.geo.nEP = self.n_ep = n_ep

    def alloc_extra_data(self, n_epv):
        if self.mode != 'surface_extra':
            raise ValueError('CMapping not in surface_extra mode!')

        self.bfg = np.empty((self.n_el, self.n_qp, self.dim, n_epv),
                            dtype=np.float64)
        array2fmfield4(self._bfg, self.bfg)
        self.geo.bfGM = self._bfg

    def __str__(self):
        return 'CMapping: mode: %s, n_el %d, n_qp %d, dim: %d, n_ep: %d' \
               % ((self.mode,) + self.shape)

    def cprint(self, int32 mode=0):
        map_print(self.geo, stdout, mode)

    def describe(self,
                 np.ndarray[float64, mode='c', ndim=2] coors not None,
                 np.ndarray[int32, mode='c', ndim=2] conn not None,
                 np.ndarray[float64, mode='c', ndim=3] bfgr not None,
                 np.ndarray[float64, mode='c', ndim=4] ebfgr,
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

        if ebfgr is None:
            ebfgr = np.array([0], ndmin=4, dtype=np.float64)

        array2fmfield3(_bfgr, bfgr)
        array2fmfield4(_ebfgr, ebfgr)
        array2fmfield1(_weights, weights)

        ret = map_describe(self.geo, _coors, n_nod, dim, _conn, n_el, n_ep,
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

        if self.mode == 'volume':
            n_row_ok = out.shape[2] == arr.shape[2]

        else:
            n_row_ok = (((mode < 3) and (out.shape[2] == arr.shape[2]))
                        or
                        ((mode >= 3) and (out.shape[2] == 1)
                         and (arr.shape[2] == self.dim)))

        if not ((out.shape[0] == arr.shape[0])
                and (out.shape[1] == 1)
                and n_row_ok
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

        ret = map_integrate(self.geo, _out, _arr, mode)

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

        ret = map_getElementDiameters(self.geo, _out,
                                      _edges, edges_nr, edges_nc,
                                      _coors, n_nod, dim,
                                      _conn, n_el, n_ep, _el_list, n_el2, mode)

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

        if self.bfg is None:
            raise ValueError('CMapping.alloc_extra_data() not called!')

        array2fmfield4(_bfbgr, bfbgr)
        array2fmfield4(_ebfbgr, ebfbgr)

        ret = map_evaluateBFBGM(self.geo, _bfbgr, _ebfbgr, _coors, n_nod, dim,
                                _fis, n_fa, n_fp, _conn, n_el, n_ep)

        if ret:
            errclear()
            raise ValueError('ccore error (see above)')
