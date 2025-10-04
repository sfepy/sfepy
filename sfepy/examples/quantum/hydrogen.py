"""
Hydrogen atom.

See :ref:`quantum-quantum_common`.
"""
from sfepy.linalg import norm_l2_along_axis

from sfepy.examples.quantum.quantum_common import common

def get_exact(n_eigs, box_size, dim):
    Z = 1
    if dim == 2:
        eigs = [-float(Z)**2/2/(n-0.5)**2/4
                for n in [1] + [2]*3 + [3]*5 + [4]*8 + [5]*15]

    elif dim == 3:
        eigs = [-float(Z)**2/2/n**2 for n in [1] + [2]*2**2 + [3]*3**2]

    return eigs

def fun_v(ts, coor, mode=None, **kwargs):
    if not mode == 'qp': return

    out = {}
    C = 0.5
    r = norm_l2_along_axis(coor, axis=1)
    V = - C * 1.0 / r

    V.shape = (V.shape[0], 1, 1)
    out['V'] = V
    return out

def define(n_eigs=5, tau=-1.0):
    l = common(fun_v, get_exact=get_exact, n_eigs=n_eigs, tau=tau)
    return l
