from quantum_common import common

def fun_v(ts, coor, mode=None, region=None, ig=None):
    from numpy import zeros_like

    if not mode == 'qp': return

    out = {}
    val = zeros_like(coor[:,0])

    val.shape = (val.shape[0], 1, 1)
    out['V'] = val
    return out

def define():
    l = common(fun_v, n_eigs=10, tau=0.0)
    return l
