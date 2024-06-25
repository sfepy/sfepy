# -*- Mode: Python -*-
# cython: language_level=3
"""
Low level reference mapping functionality.
"""

cdef class CMapping:
    def __cinit__(self, bf, det, volume, bfg, normal, dim):
        array2fmfield4(self._bf, bf)
        self.geo.bf = self._bf

        array2fmfield4(self._det, det)
        self.geo.det = self._det

        array2fmfield4(self._volume, volume)
        self.geo.volume = self._volume

        if bfg is not None:
            array2fmfield4(self._bfg, bfg)
            self.geo.bfGM = self._bfg
        else:
            self.geo.bfGM = NULL

        if normal is not None:
            array2fmfield4(self._normal, normal)
            self.geo.normal = self._normal
        else:
            self.geo.normal = NULL

        n_el, n_qp = det.shape[:2]
        n_ep = bf.shape[3]
        self.geo.nEl = n_el
        self.geo.nQP = n_qp
        self.geo.dim = dim
        self.geo.nEP = n_ep

    def __str__(self):
        return 'CMapping: n_el %d, n_qp %d, dim: %d, n_ep: %d' \
               % (self.geo.nEl, self.geo.nQP, self.geo.dim, self.geo.nEP)
