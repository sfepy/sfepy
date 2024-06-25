# -*- Mode: Python -*-
# cython: language_level=3
"""
Low level functions.
"""
from sfepy.discrete.common.extmods._fmfield cimport FMField
from sfepy.discrete.common.extmods.types cimport int32

cdef extern from 'geommech.h':
     int32 geme_mulAVSB3(FMField *out, FMField *vs, FMField *inp)
