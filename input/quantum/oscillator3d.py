from hydrogen import common

def define():
    l = common(dim=3, n_eigs=20, tau=0.0)
    return l

def fun_v( ts, coor, region, ig, mode = None ):
    out = {}
    C = 0.5
    val = C * (coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2)
    out['V'] = val
    return out
