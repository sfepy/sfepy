# -*- Mode: Python -*-
"""
Low level functions.
"""
from sfepy.discrete.common.extmods._fmfield cimport FMField
from types cimport int32

cdef extern from 'geommech.h':
     int32 geme_mulAVSB3(FMField *out, FMField *vs, FMField *inp)
