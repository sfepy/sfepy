# -*- Mode: Python -*-
# cython: language_level=3
"""
Low level functions.
"""
cimport cython

cimport sfepy.discrete.common.extmods._geommech as gmch
from sfepy.discrete.common.extmods.types cimport int32

from sfepy.discrete.common.extmods._fmfield cimport (
    FMField, FMF_SetCell, array2fmfield4,
)

cimport numpy as np
import numpy as np

def geme_mulAVSB3py(np.ndarray vs not None,
                   np.ndarray inp not None):
    cdef int32 ret
    cdef FMField[1] _out, _vs, _inp

    out = np.zeros_like(inp)
    array2fmfield4(_out, out)
    array2fmfield4(_vs, vs)
    array2fmfield4(_inp, inp)

    ret = 0
    for ii in range(out.shape[0]):
        FMF_SetCell(_out, ii);
        FMF_SetCell(_vs, ii);
        FMF_SetCell(_inp, ii);
        _ret = gmch.geme_mulAVSB3(_out, _vs, _inp)
        ret = ret and _ret

    return out, ret
