from hydrogen import common

def define():
    l = common(dim=2, n_eigs=10, tau=0.0)
    return l

def funV( ts, coor, region, ig, mode = None ):
    from numpy import zeros_like, where

    out = {}
    val = zeros_like( coor[:,0] )
    out['V'] = val
    return out
