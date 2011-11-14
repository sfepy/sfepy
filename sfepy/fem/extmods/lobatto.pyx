# -*- Mode: Python -*-
"""
Template file for 'lobatto.pyx'.
"""
cimport cython

import numpy as np
cimport numpy as np

from types cimport int32, float64, complex128

cdef extern from 'math.h':
    cdef float64 sqrt(float64 x)
    cdef float64 pow(float64 x, float64 y)

ctypedef float64 (*fun)(float64 x)

# Start of generated code.

# Lobatto functions.

cdef float64 lobatto_000(float64 x):
    return -x/2 + 1.0/2.0

cdef float64 lobatto_001(float64 x):
    return x/2 + 1.0/2.0

cdef float64 lobatto_002(float64 x):
    return sqrt(6)*pow(x, 2)/4 - sqrt(6)/4

cdef float64 lobatto_003(float64 x):
    return x*(sqrt(10)*pow(x, 2) - sqrt(10))/4

cdef float64 lobatto_004(float64 x):
    return pow(x, 2)*(5*sqrt(14)*pow(x, 2) - 6*sqrt(14))/16 + sqrt(14)/16

cdef float64 lobatto_005(float64 x):
    return x*(pow(x, 2)*(21*sqrt(2)*pow(x, 2) - 30*sqrt(2)) + 9*sqrt(2))/16

cdef float64 lobatto_006(float64 x):
    return pow(x, 2)*(pow(x, 2)*(21*sqrt(22)*pow(x, 2) - 35*sqrt(22)) + 15*sqrt(22))/32 - sqrt(22)/32

cdef float64 lobatto_007(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(33*sqrt(26)*pow(x, 2) - 63*sqrt(26)) + 35*sqrt(26)) - 5*sqrt(26))/32

cdef float64 lobatto_008(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(429*sqrt(30)*pow(x, 2) - 924*sqrt(30)) + 630*sqrt(30)) - 140*sqrt(30))/256 + 5*sqrt(30)/256

cdef float64 lobatto_009(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(715*sqrt(34)*pow(x, 2) - 1716*sqrt(34)) + 1386*sqrt(34)) - 420*sqrt(34)) + 35*sqrt(34))/256

cdef float64 lobatto_010(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(2431*sqrt(38)*pow(x, 2) - 6435*sqrt(38)) + 6006*sqrt(38)) - 2310*sqrt(38)) + 315*sqrt(38))/512 - 7*sqrt(38)/512

# Kernel functions.

cdef float64 kernel_002(float64 x):
    return -sqrt(6)

cdef float64 kernel_003(float64 x):
    return -sqrt(10)*x

cdef float64 kernel_004(float64 x):
    return -5*sqrt(14)*pow(x, 2)/4 + sqrt(14)/4

cdef float64 kernel_005(float64 x):
    return x*(-21*sqrt(2)*pow(x, 2) + 9*sqrt(2))/4

cdef float64 kernel_006(float64 x):
    return pow(x, 2)*(-21*sqrt(22)*pow(x, 2) + 14*sqrt(22))/8 - sqrt(22)/8

cdef float64 kernel_007(float64 x):
    return x*(pow(x, 2)*(-33*sqrt(26)*pow(x, 2) + 30*sqrt(26)) - 5*sqrt(26))/8

cdef float64 kernel_008(float64 x):
    return pow(x, 2)*(pow(x, 2)*(-429*sqrt(30)*pow(x, 2) + 495*sqrt(30)) - 135*sqrt(30))/64 + 5*sqrt(30)/64

cdef float64 kernel_009(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(-715*sqrt(34)*pow(x, 2) + 1001*sqrt(34)) - 385*sqrt(34)) + 35*sqrt(34))/64

cdef float64 kernel_010(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(-2431*sqrt(38)*pow(x, 2) + 4004*sqrt(38)) - 2002*sqrt(38)) + 308*sqrt(38))/128 - 7*sqrt(38)/128

# Legendre functions.

cdef float64 legendre_000(float64 x):
    return 1

cdef float64 legendre_001(float64 x):
    return x

cdef float64 legendre_002(float64 x):
    return 3*pow(x, 2)/2 - 1.0/2.0

cdef float64 legendre_003(float64 x):
    return x*(5*pow(x, 2) - 3)/2

cdef float64 legendre_004(float64 x):
    return pow(x, 2)*(35*pow(x, 2) - 30)/8 + 3.0/8.0

cdef float64 legendre_005(float64 x):
    return x*(pow(x, 2)*(63*pow(x, 2) - 70) + 15)/8

cdef float64 legendre_006(float64 x):
    return pow(x, 2)*(pow(x, 2)*(231*pow(x, 2) - 315) + 105)/16 - 5.0/16.0

cdef float64 legendre_007(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(429*pow(x, 2) - 693) + 315) - 35)/16

cdef float64 legendre_008(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(6435*pow(x, 2) - 12012) + 6930) - 1260)/128 + 35.0/128.0

cdef float64 legendre_009(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(12155*pow(x, 2) - 25740) + 18018) - 4620) + 315)/128

cdef float64 legendre_010(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(46189*pow(x, 2) - 109395) + 90090) - 30030) + 3465)/256 - 63.0/256.0

# Lists of functions.

cdef int32 max_order = 10

cdef fun lobatto[11]

lobatto[:] = [&lobatto_000, &lobatto_001, &lobatto_002, &lobatto_003, &lobatto_004, &lobatto_005, &lobatto_006, &lobatto_007, &lobatto_008, &lobatto_009, &lobatto_010]

cdef fun kernel[9]

kernel[:] = [&kernel_002, &kernel_003, &kernel_004, &kernel_005, &kernel_006, &kernel_007, &kernel_008, &kernel_009, &kernel_010]

cdef fun legendre[11]

legendre[:] = [&legendre_000, &legendre_001, &legendre_002, &legendre_003, &legendre_004, &legendre_005, &legendre_006, &legendre_007, &legendre_008, &legendre_009, &legendre_010]

# End of generated code.

@cython.boundscheck(False)
def eval_lobatto(np.ndarray[float64, mode='c', ndim=1] coors not None,
                 int32 order):
    """
    Evaluate Lobatto function of the given order in given points.
    """
    cdef int32 ii
    cdef fun eval_fun
    cdef int32 n_coor = coors.shape[0]
    cdef np.ndarray[float64, ndim=1] out = np.zeros(n_coor, dtype=np.float64)
    cdef float64 *_coors = &coors[0]
    cdef float64 *_out = &out[0]

    if (order < 0) or (order > max_order):
        raise ValueError('order must be in [0, %d]! (was %d)'
                         % (max_order, order))

    eval_fun = lobatto[order]
    for ii in range(0, n_coor):
        _out[ii] = eval_fun(_coors[ii])

    return out
