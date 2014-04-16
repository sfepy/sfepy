from quantum_common import common

def fun_v(ts, coor, mode=None, **kwargs):
    from numpy import zeros_like

    if not mode == 'qp': return

    out = {}
    val = zeros_like(coor[:,0])

    val.shape = (val.shape[0], 1, 1)
    out['V'] = val
    return out

def define(n_eigs=10, tau=0.0):
    l = common(fun_v, n_eigs=n_eigs, tau=tau)
    return l
