# -*- Mode: Python -*-
# cython: language_level=3
"""
Basic types.
"""
cimport numpy as np
import numpy as np

cdef extern from "numpy/arrayobject.h":
    cdef int NPY_NO_DEPRECATED_API = 0

ctypedef np.complex128_t complex128
ctypedef np.float64_t float64
ctypedef np.int32_t int32
ctypedef np.uint32_t uint32
