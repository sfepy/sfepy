# -*- Mode: Python -*-
"""
Low level functions.
"""
cimport cython

cimport sfepy.discrete.fem.extmods._geommech as gmch
from types cimport int32

from sfepy.discrete.fem.extmods._fmfield cimport (FMField, array2fmfield4)

cimport numpy as np
import numpy as np

def geme_mulAVSB3py(np.ndarray vs not None,
                   np.ndarray inp not None):
    cdef int32 ret
    cdef FMField _out[1], _vs[1], _inp[1]

    out = np.zeros_like(inp)
    array2fmfield4(_out, out)
    array2fmfield4(_vs, vs)
    array2fmfield4(_inp, inp)

    ret = gmch.geme_mulAVSB3(_out, _vs, _inp)

    return out, ret
