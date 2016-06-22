from __future__ import absolute_import
from sfepy.linalg import norm_l2_along_axis

from examples.quantum.quantum_common import common

def fun_v(ts, coor, mode=None, **kwargs):
    if not mode == 'qp': return

    out = {}
    C = 0.5
    r = norm_l2_along_axis(coor, axis=1)
    V = - C * 5.0 / r

    V.shape = (V.shape[0], 1, 1)
    out['V'] = V
    return out

def define(n_eigs=10, tau=-15):
    l = common(fun_v, n_eigs=n_eigs, tau=tau)
    return l
