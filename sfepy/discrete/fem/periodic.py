import numpy as nm

from sfepy.discrete.fem.mesh import find_map

periodic_cache = {}

##
# c: 05.05.2008, r: 05.05.2008
eps = 1e-9
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


def match_plane_by_dir(coors1, coors2, direction, get_saved=True):
    """
    Match coordinates `coors1` with `coors2` in a given direction.
    """
    if coors1.shape != coors2.shape:
        raise ValueError('incompatible shapes: %s == %s'\
                         % (coors1.shape, coors2.shape))

    key_dir = None if direction is None else tuple(direction)
    key = coors1.shape + ('dir', key_dir)
    if key in periodic_cache and get_saved:
        return periodic_cache[key]
    else:
        aux = coors2.copy()
        if direction is not None:
            direction = nm.asarray(direction, dtype=nm.float64).reshape((3, 1))
            direction = direction[:coors1.shape[1]]
            direction /= nm.linalg.norm(direction)
            x0p = coors1[0] - coors2
            dist = nm.linalg.norm(x0p - nm.dot(x0p, direction) * direction.T,
                                  axis=1)
            idx = nm.where(dist < eps)[0]
            if len(idx) < 1:
                print(direction)
                print(dist)
                raise ValueError('cannot match nodes!')

            aux += coors1[0] - coors2[idx[0]]

        i1, i2 = find_map(coors1, aux, join=False)

        if i1.shape[0] != coors1.shape[0]:
            print(direction)
            print(coors1[i1])
            print(coors2[i2])
            print(nm.abs(coors1[i1] - coors2[i2]).max(0))
            ii = nm.setdiff1d(nm.arange(coors1.shape[0]), i1)
            print(coors1[ii])
            print(coors2[ii])
            raise ValueError('cannot match nodes!')

        periodic_cache[key] = (i1, i2)

        return i1, i2

def get_grid_plane(idim):
    return nm.eye(3)[idim]

def match_x_plane(coors1, coors2, get_saved=True):
    return match_plane_by_dir(coors1, coors2, get_grid_plane(0), get_saved)
def match_y_plane(coors1, coors2, get_saved=True):
    return match_plane_by_dir(coors1, coors2, get_grid_plane(1), get_saved)
def match_z_plane(coors1, coors2, get_saved=True):
    return match_plane_by_dir(coors1, coors2, get_grid_plane(2), get_saved)
def match_grid_plane(coors1, coors2, idim, get_saved=True):
    return match_plane_by_dir(coors1, coors2, get_grid_plane(idim), get_saved)
def match_coors(coors1, coors2, get_saved=True):
    return match_plane_by_dir(coors1, coors2, None, get_saved)