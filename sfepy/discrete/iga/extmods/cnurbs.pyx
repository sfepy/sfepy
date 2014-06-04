# -*- Mode: Python -*-
cimport cython

import numpy as np
cimport numpy as np

from sfepy.discrete.fem.extmods.types cimport int32, uint32, float64

from sfepy.discrete.fem.extmods._fmfield cimport (FMField, array2fmfield1)

cdef extern from 'nurbs.h':
    cdef int32 _eval_bernstein_basis \
         'eval_bernstein_basis'(FMField *funs, FMField *ders,
                                float64 x, uint32 degree)

def eval_bernstein_basis(np.ndarray funs not None,
                         np.ndarray ders not None,
                         float64 x,
                         uint32 degree):
    cdef int32 ret
    cdef FMField _funs[1], _ders[1]

    array2fmfield1(_funs, funs)
    array2fmfield1(_ders, ders)

    ret = _eval_bernstein_basis(_funs, _ders, x, degree)
    return ret
