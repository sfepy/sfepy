# -*- Mode: Python -*-
"""
Template file for 'lobatto.pyx'.
"""
cimport cython

from types cimport int32, float64, complex128

cdef extern from 'math.h':
    cdef float64 sqrt(float64 x)
    cdef float64 pow(float64 x, float64 y)

ctypedef float64 (*fun)(float64 x)

# Start of generated code.
%s
# End of generated code.
