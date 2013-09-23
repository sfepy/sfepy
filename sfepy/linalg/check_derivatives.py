"""
Utilities for checking derivatives of functions.
"""

import numpy as nm

def check_fx(x0, fx, fx_args, dfx, dfx_args=None, delta=1e-5):
    """
    Check derivatives of a (vectorized) scalar function of a scalar variable.
    """
    if dfx_args is None:
        dfx_args = fx_args

    dfx_a = dfx(x0, *dfx_args)

    x = x0 + delta
    f1 = fx(x, *fx_args)

    x = x0 - delta
    f2 = fx(x, *fx_args)

    dfx_d = 0.5 * (f1 - f2) / delta

    error = nm.linalg.norm(dfx_a - dfx_d, nm.inf)

    print 'analytical:', dfx_a
    print 'difference:', dfx_d
    print 'error:', error

    return dfx_a, dfx_d, error
