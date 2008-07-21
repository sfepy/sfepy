from dft import common

def define():
    l = common(dim=2, n_eigs=7, tau=-13)
    return l

def fun_v( ts, coor, region, ig, mode = None, vhxc = None ):
    import numpy as nm

    if vhxc is None:
        vhxc = 0.0

    out = {}
    C = 0.5
    r = nm.sqrt( coor[:,0]**2 + coor[:,1]**2 )
    vc = - C * 5.0 / r
    V = vhxc + vc

    out['V'] = V
    return out
