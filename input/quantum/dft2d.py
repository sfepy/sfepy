from dft import common

n_electron = 7
n_eigs = None

def define():
    l = common(dim=2, n_eigs=n_eigs, n_electron=n_electron, tau=-10)
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

def core_pot(bc, ts, coor):
    from sfepy.base.la import norm_l2_along_axis
    
    r = norm_l2_along_axis(coor)
    return n_electron / r
