"""
Quantum oscillator.

See :ref:`quantum-quantum_common`.
"""
from sfepy.linalg import norm_l2_along_axis

from sfepy.examples.quantum.quantum_common import common

def get_exact(n_eigs, box_size, dim):
    if dim == 2:
        eigs = [1] + [2]*2 + [3]*3 + [4]*4 + [5]*5 + [6]*6

    elif dim == 3:
        eigs = [float(1)/2 + x for x in [1] + [2]*3 + [3]*6 + [4]*10]

    return eigs

def fun_v(ts, coor, mode=None, **kwargs):
    if not mode == 'qp': return

    out = {}
    C = 0.5
    val = C * norm_l2_along_axis(coor, axis=1, squared=True)

    val.shape = (val.shape[0], 1, 1)
    out['V'] = val
    return out

def define(n_eigs=20, tau=0.0):
    l = common(fun_v, get_exact=get_exact, n_eigs=n_eigs, tau=tau)
    return l
