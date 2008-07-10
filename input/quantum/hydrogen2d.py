from hydrogen import common

def define():
    l = common('tmp/mesh.vtk', 2)
    return l

def funV( ts, coor, region, ig, mode = None ):
    from numpy import sqrt

    out = {}
    C = 0.5
    r = sqrt( coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2 )
    V = - C * 1.0 / r
    out['V'] = V
    return out
