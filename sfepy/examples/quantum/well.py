"""
Quantum potential well.

See :ref:`quantum-quantum_common`.
"""
from __future__ import absolute_import

from sfepy.examples.quantum.quantum_common import common

def get_exact(n_eigs, box_size, dim):
    from numpy import pi

    if dim == 2:
        eigs = [pi**2/(2*box_size**2)*x
                for x in [2, 5, 5, 8, 10, 10, 13, 13, 17, 17, 18, 20, 20]]

    elif dim == 3:
        eigs = [pi**2/(2*box_size**2)*x
                for x in [3, 6, 6, 6, 9, 9, 9, 11, 11, 11,
                          12, 14, 14, 14, 14, 14, 14, 17, 17, 17]]

    return eigs

def fun_v(ts, coor, mode=None, **kwargs):
    from numpy import zeros_like

    if not mode == 'qp': return

    out = {}
    val = zeros_like(coor[:,0])

    val.shape = (val.shape[0], 1, 1)
    out['V'] = val
    return out

def define(n_eigs=10, tau=0.0):
    l = common(fun_v, get_exact=get_exact, n_eigs=n_eigs, tau=tau)
    return l
