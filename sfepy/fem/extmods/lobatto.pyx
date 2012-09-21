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

# Derivatives of Lobatto functions.

cdef float64 d_lobatto_000(float64 x):
    return -1.0/2.0

cdef float64 d_lobatto_001(float64 x):
    return 1.0/2.0

cdef float64 d_lobatto_002(float64 x):
    return sqrt(6)*x/2

cdef float64 d_lobatto_003(float64 x):
    return 3*sqrt(10)*pow(x, 2)/4 - sqrt(10)/4

cdef float64 d_lobatto_004(float64 x):
    return x*(20*sqrt(14)*pow(x, 2) - 12*sqrt(14))/16

cdef float64 d_lobatto_005(float64 x):
    return pow(x, 2)*(105*sqrt(2)*pow(x, 2) - 90*sqrt(2))/16 + 9*sqrt(2)/16

cdef float64 d_lobatto_006(float64 x):
    return x*(pow(x, 2)*(126*sqrt(22)*pow(x, 2) - 140*sqrt(22)) + 30*sqrt(22))/32

cdef float64 d_lobatto_007(float64 x):
    return pow(x, 2)*(pow(x, 2)*(231*sqrt(26)*pow(x, 2) - 315*sqrt(26)) + 105*sqrt(26))/32 - 5*sqrt(26)/32

cdef float64 d_lobatto_008(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(3432*sqrt(30)*pow(x, 2) - 5544*sqrt(30)) + 2520*sqrt(30)) - 280*sqrt(30))/256

cdef float64 d_lobatto_009(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(6435*sqrt(34)*pow(x, 2) - 12012*sqrt(34)) + 6930*sqrt(34)) - 1260*sqrt(34))/256 + 35*sqrt(34)/256

cdef float64 d_lobatto_010(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(24310*sqrt(38)*pow(x, 2) - 51480*sqrt(38)) + 36036*sqrt(38)) - 9240*sqrt(38)) + 630*sqrt(38))/512

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

# Derivatives of kernel functions.

cdef float64 d_kernel_002(float64 x):
    return 0

cdef float64 d_kernel_003(float64 x):
    return -sqrt(10)

cdef float64 d_kernel_004(float64 x):
    return -5*sqrt(14)*x/2

cdef float64 d_kernel_005(float64 x):
    return -63*sqrt(2)*pow(x, 2)/4 + 9*sqrt(2)/4

cdef float64 d_kernel_006(float64 x):
    return x*(-84*sqrt(22)*pow(x, 2) + 28*sqrt(22))/8

cdef float64 d_kernel_007(float64 x):
    return pow(x, 2)*(-165*sqrt(26)*pow(x, 2) + 90*sqrt(26))/8 - 5*sqrt(26)/8

cdef float64 d_kernel_008(float64 x):
    return x*(pow(x, 2)*(-2574*sqrt(30)*pow(x, 2) + 1980*sqrt(30)) - 270*sqrt(30))/64

cdef float64 d_kernel_009(float64 x):
    return pow(x, 2)*(pow(x, 2)*(-5005*sqrt(34)*pow(x, 2) + 5005*sqrt(34)) - 1155*sqrt(34))/64 + 35*sqrt(34)/64

cdef float64 d_kernel_010(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(-19448*sqrt(38)*pow(x, 2) + 24024*sqrt(38)) - 8008*sqrt(38)) + 616*sqrt(38))/128

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

# Derivatives of Legendre functions.

cdef float64 d_legendre_000(float64 x):
    return 0

cdef float64 d_legendre_001(float64 x):
    return 1

cdef float64 d_legendre_002(float64 x):
    return 3*x

cdef float64 d_legendre_003(float64 x):
    return 15*pow(x, 2)/2 - 3.0/2.0

cdef float64 d_legendre_004(float64 x):
    return x*(140*pow(x, 2) - 60)/8

cdef float64 d_legendre_005(float64 x):
    return pow(x, 2)*(315*pow(x, 2) - 210)/8 + 15.0/8.0

cdef float64 d_legendre_006(float64 x):
    return x*(pow(x, 2)*(1386*pow(x, 2) - 1260) + 210)/16

cdef float64 d_legendre_007(float64 x):
    return pow(x, 2)*(pow(x, 2)*(3003*pow(x, 2) - 3465) + 945)/16 - 35.0/16.0

cdef float64 d_legendre_008(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(51480*pow(x, 2) - 72072) + 27720) - 2520)/128

cdef float64 d_legendre_009(float64 x):
    return pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(109395*pow(x, 2) - 180180) + 90090) - 13860)/128 + 315.0/128.0

cdef float64 d_legendre_010(float64 x):
    return x*(pow(x, 2)*(pow(x, 2)*(pow(x, 2)*(461890*pow(x, 2) - 875160) + 540540) - 120120) + 6930)/256

# Lists of functions.

cdef int32 max_order = 10

cdef fun lobatto[11]

lobatto[:] = [&lobatto_000, &lobatto_001, &lobatto_002, &lobatto_003, &lobatto_004, &lobatto_005, &lobatto_006, &lobatto_007, &lobatto_008, &lobatto_009, &lobatto_010]

cdef fun d_lobatto[11]

d_lobatto[:] = [&d_lobatto_000, &d_lobatto_001, &d_lobatto_002, &d_lobatto_003, &d_lobatto_004, &d_lobatto_005, &d_lobatto_006, &d_lobatto_007, &d_lobatto_008, &d_lobatto_009, &d_lobatto_010]

cdef fun kernel[9]

kernel[:] = [&kernel_002, &kernel_003, &kernel_004, &kernel_005, &kernel_006, &kernel_007, &kernel_008, &kernel_009, &kernel_010]

cdef fun d_kernel[9]

d_kernel[:] = [&d_kernel_002, &d_kernel_003, &d_kernel_004, &d_kernel_005, &d_kernel_006, &d_kernel_007, &d_kernel_008, &d_kernel_009, &d_kernel_010]

cdef fun legendre[11]

legendre[:] = [&legendre_000, &legendre_001, &legendre_002, &legendre_003, &legendre_004, &legendre_005, &legendre_006, &legendre_007, &legendre_008, &legendre_009, &legendre_010]

cdef fun d_legendre[11]

d_legendre[:] = [&d_legendre_000, &d_legendre_001, &d_legendre_002, &d_legendre_003, &d_legendre_004, &d_legendre_005, &d_legendre_006, &d_legendre_007, &d_legendre_008, &d_legendre_009, &d_legendre_010]

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

@cython.boundscheck(False)
def eval_lobatto_tensor_product(np.ndarray[float64, mode='c', ndim=2]
                                coors not None,
                                np.ndarray[int32, mode='c', ndim=2]
                                nodes not None,
                                float64 cmin, float64 cmax, int32 order,
                                int32 diff=False):
    cdef np.ndarray[float64, ndim=3] out
    cdef np.ndarray[float64, ndim=2] lambdas
    cdef np.ndarray[float64, ndim=2] xis
    cdef float64 dx = cmax - cmin
    cdef int32 ii, ifun, ic, ir, io, nr
    cdef int32 n_coor = coors.shape[0]
    cdef int32 dim = coors.shape[1]
    cdef int32 n_fun = nodes.shape[0]
    cdef float64 *_xis, *_out
    cdef int32 *_nodes = &nodes[0, 0]
    cdef fun eval_fun

    if (order < 1) or (order > max_order):
        raise ValueError('order must be in [1, %d]! (was %d)'
                         % (max_order, order))

    nr = 1 if not diff else dim
    out = np.ones((n_coor, nr, n_fun), dtype=np.float64)

    # Transform coordinates via affine coordinates lambda to be in [-1, 1].
    lambdas = (coors - cmin) / dx
    xis = 2.0 * lambdas - 1.0
    _xis = &xis[0, 0]
    _out = &out[0, 0, 0]
    if not diff:
        for ii in range(0, dim):
            for ifun in range(0, n_fun):
                eval_fun = lobatto[_nodes[dim * ifun + ii]]
                for ic in range(0, n_coor):
                    _out[n_fun * ic + ifun] *= eval_fun(_xis[dim * ic + ii])

    else:
        for ii in range(0, dim):
            for ifun in range(0, n_fun):
                for ir in range(0, dim):
                    if ir == ii:
                        eval_fun = d_lobatto[_nodes[dim * ifun + ii]]

                    else:
                        eval_fun = lobatto[_nodes[dim * ifun + ii]]

                    for ic in range(0, n_coor):
                        io = n_fun * (nr * ic + ir) + ifun
                        _out[io] *= eval_fun(_xis[dim * ic + ii])

        # Multiply by 2 due to the transformation of coordinates.
        for ii in range(0, n_coor * nr * n_fun):
            _out[ii] *= 2.0

    return out
