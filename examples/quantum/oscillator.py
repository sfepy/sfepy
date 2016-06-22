from __future__ import absolute_import
from sfepy.linalg import norm_l2_along_axis

from examples.quantum.quantum_common import common

def fun_v(ts, coor, mode=None, **kwargs):
    if not mode == 'qp': return

    out = {}
    C = 0.5
    val = C * norm_l2_along_axis(coor, axis=1, squared=True)

    val.shape = (val.shape[0], 1, 1)
    out['V'] = val
    return out

def define(n_eigs=20, tau=0.0):
    l = common(fun_v, n_eigs=n_eigs, tau=tau)
    return l
