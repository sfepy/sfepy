# -*- Mode: Python -*-
"""
Low level reference mapping functionality.
"""

cdef class CMapping:
    def __cinit__(self, bf, det, volume, bfg, normal, dim):
        self.bf = bf
        array2fmfield4(self._bf, bf)
        self.geo.bf = self._bf

        self.det = det
        array2fmfield4(self._det, det)
        self.geo.det = self._det

        self.volume = volume
        array2fmfield4(self._volume, volume)
        self.geo.volume = self._volume

        self.bfg = bfg
        if bfg is not None:
            array2fmfield4(self._bfg, bfg)
            self.geo.bfGM = self._bfg
        else:
            self.geo.bfGM = NULL

        self.normal = normal
        if normal is not None:
            array2fmfield4(self._normal, normal)
            self.geo.normal = self._normal
        else:
            self.geo.normal = NULL

        n_el, n_qp = det.shape[:2]
        n_ep = bf.shape[3]
        self.geo.nEl = self.n_el = n_el
        self.geo.nQP = self.n_qp = n_qp
        self.geo.dim = self.dim = dim
        self.geo.nEP = self.n_ep = n_ep

    def __str__(self):
        return 'CMapping: n_el %d, n_qp %d, dim: %d, n_ep: %d' \
               % (self.n_el, self.n_qp, self.dim, self.n_ep)
