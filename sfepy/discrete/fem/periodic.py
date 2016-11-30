from __future__ import print_function
import numpy as nm

from sfepy.discrete.fem.mesh import find_map

periodic_cache = {}

##
# c: 05.05.2008, r: 05.05.2008
eps = 1e-12
def set_accuracy(eps):
    globals()['eps'] = eps

##
# c: 18.10.2006, r: 05.05.2008
def match_grid_line(coors1, coors2, which, get_saved=True):
    """
    Match coordinates `coors1` with `coors2` along the axis `which`.
    """
    if coors1.shape != coors2.shape:
        raise ValueError('incompatible shapes: %s == %s'\
              % (coors1.shape, coors2.shape))

    key = coors1.shape + ('line', which)
    if key in periodic_cache and get_saved:
        return periodic_cache[key]
    else:
        c1 = coors1[:,which]
        c2 = coors2[:,which]
        i1 = nm.argsort(c1)
        i2 = nm.argsort(c2)

        if not nm.all(nm.abs(c1[i1] - c2[i2]) < eps):
            print(c1[i1])
            print(c2[i2])
            print(nm.abs(c1[i1] - c2[i2]).max())
            raise ValueError('cannot match nodes!')

        periodic_cache[key] = (i1, i2)

        return i1, i2

##
# 18.10.2006, c
# last revision: 18.10.2006
def match_x_line(coors1, coors2, get_saved=True):
    return match_grid_line(coors1, coors2, 0, get_saved)
def match_y_line(coors1, coors2, get_saved=True):
    return match_grid_line(coors1, coors2, 1, get_saved)
def match_z_line(coors1, coors2, get_saved=True):
    return match_grid_line(coors1, coors2, 2, get_saved)

##
# 01.06.2007, c
# last revision: 01.06.2007
def match_grid_plane(coors1, coors2, which, get_saved=True):
    """
    Match coordinates `coors1` with `coors2` along the plane with normal axis
    `which`.
    """
    if coors1.shape != coors2.shape:
        raise ValueError('incompatible shapes: %s == %s'\
              % (coors1.shape, coors2.shape))

    key = coors1.shape + ('plane', which)
    if key in periodic_cache and get_saved:
        return periodic_cache[key]
    else:
        offset = coors1[0,which] - coors2[0,which]
        aux = coors2.copy()
        aux[:,which] += offset
        i1, i2 = find_map(coors1, aux, join = False)

        if i1.shape[0] != coors1.shape[0]:
            print(coors1[i1])
            print(coors2[i2])
            print(nm.abs(coors1[i1] - coors2[i2]).max(0))
            ii = nm.setdiff1d(nm.arange(coors1.shape[0]), i1)
            print(coors1[ii])
            print(coors2[ii])
            raise ValueError('cannot match nodes!')

        periodic_cache[key] = (i1, i2)

        return i1, i2

##
# 01.06.2007, c
# last revision: 01.06.2007
def match_x_plane(coors1, coors2, get_saved=True):
    return match_grid_plane(coors1, coors2, 0, get_saved)
def match_y_plane(coors1, coors2, get_saved=True):
    return match_grid_plane(coors1, coors2, 1, get_saved)
def match_z_plane(coors1, coors2, get_saved=True):
    return match_grid_plane(coors1, coors2, 2, get_saved)

def match_coors(coors1, coors2, get_saved=True):
    """
    Match coordinates `coors1` with `coors2`.
    """
    if coors1.shape != coors2.shape:
        raise ValueError('incompatible shapes: %s == %s'\
                         % (coors1.shape, coors2.shape))

    key = coors1.shape + ('coors',)
    if key in periodic_cache and get_saved:
        return periodic_cache[key]
    else:
        i1, i2 = find_map(coors1, coors2, join=False)

        if i1.shape[0] != coors1.shape[0]:
            print(coors1[i1])
            print(coors2[i2])
            print(nm.abs(coors1[i1] - coors2[i2]).max(0))
            ii = nm.setdiff1d(nm.arange(coors1.shape[0]), i1)
            print(coors1[ii])
            print(coors2[ii])
            raise ValueError('cannot match nodes!')

        periodic_cache[key] = (i1, i2)

        return i1, i2
