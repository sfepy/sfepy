import numpy as nm

from sfepy.base.base import assert_, output
from sfepy.base.la import norm_l2_along_axis as norm

def transform_bar_to_space_coors(bar_coors, coors):
    """
    Transform barycentric coordinates `bar_coors` within simplices with
    vertex coordinates `coors` to space coordinates.
    """
    space_coors = nm.zeros((coors.shape[0], coors.shape[2]), dtype=coors.dtype)
    for iv in range(bar_coors.shape[1]):
        space_coors += bar_coors[:,iv:iv+1] * coors[:,iv,:]

    return space_coors

def get_simplex_circumcentres(coors, force_inside_eps=None):
    """
    Compute the circumcentres of `n_s` simplices in 1D, 2D and 3D.

    Parameters
    ----------
    coors : array
        The coordinates of the simplices with `n_v` vertices given in an
        array of shape `(n_s, n_v, dim)`, where `dim` is the space
        dimension and `2 <= n_v <= (dim + 1)`.
    force_inside_eps : float, optional
        If not None, move the circumcentres that are outside of their
        simplices or closer to their boundary then `force_inside_eps` so
        that they are inside the simplices at the distance given by
        `force_inside_eps`. It is ignored for edges.

    Returns
    -------
    centres : array
        The circumcentre coordinates as an array of shape `(n_s, dim)`.
    """
    n_s, n_v, dim = coors.shape

    assert_(2 <= n_v <= (dim + 1))
    assert_(1 <= dim <= 3)

    if n_v == 2: # Edges.
        centres = 0.5 * nm.sum(coors, axis=1)

    else:
        if n_v == 3: # Triangles.
            a2 = norm(coors[:,1,:] - coors[:,2,:], squared=True)
            b2 = norm(coors[:,0,:] - coors[:,2,:], squared=True)
            c2 = norm(coors[:,0,:] - coors[:,1,:], squared=True)

            bar_coors = nm.c_[a2 * (-a2 + b2 + c2),
                              b2 * (a2 - b2 + c2),
                              c2 * (a2 + b2 - c2)]

        elif n_v == 4: # Tetrahedrons.
            a2 = norm(coors[:,2,:] - coors[:,1,:], squared=True)
            b2 = norm(coors[:,2,:] - coors[:,0,:], squared=True)
            c2 = norm(coors[:,1,:] - coors[:,0,:], squared=True)
            d2 = norm(coors[:,3,:] - coors[:,0,:], squared=True)
            e2 = norm(coors[:,3,:] - coors[:,1,:], squared=True)
            f2 = norm(coors[:,3,:] - coors[:,2,:], squared=True)
            bar_coors = nm.c_[(d2 * a2 * (f2 + e2 - a2)
                               + b2 * e2 * (a2 + f2 - e2)
                               + c2 * f2 * (e2 + a2 - f2)
                               - 2 * a2 * e2 * f2),
                              (e2 * b2 * (f2 + d2 - b2)
                               +  c2 * f2 * (d2 + b2 - f2)
                               +  a2 * d2 * (b2 + f2 - d2)
                               - 2 * b2 * d2 * f2),
                              (f2 * c2 * (e2 + d2 - c2)
                               +  b2 * e2 * (d2 + c2 - e2)
                               +  a2 * d2 * (c2 + e2 - d2)
                               - 2 * c2 * e2 * d2),
                              (d2 * a2 * (b2 + c2 - a2)
                               +  e2 * b2 * (c2 + a2 - b2)
                               +  f2 * c2 * (a2 + b2 - c2)
                               - 2 * a2 * b2 * c2)]

        else:
            raise ValueError('unsupported simplex! (%d vertices)' % n_v)

        bar_coors /= nm.sum(bar_coors, axis=1)[:,None]
        if force_inside_eps is not None:
            bc = 1.0 / n_v
            limit = 0.9 * bc
            bar_centre = nm.array([bc] * n_v, dtype=nm.float64)

            eps = float(force_inside_eps)
            if eps > limit:
                output('force_inside_eps is too big, adjusting! (%e -> %e)'
                       % (eps, limit))
                eps = limit

            # Flag is True where the barycentre is closer to the simplex
            # boundary then eps, or outside of the simplex.
            mb = nm.min(bar_coors, axis=1)
            flag = nm.where(mb < eps)[0]

            # Move the bar_coors[flag] towards bar_centre so that it is
            # inside at the eps distance.
            mb = mb[flag]
            alpha = ((eps - mb) / (bar_centre[0] - mb))[:,None]

            bar_coors[flag] = (1.0 - alpha) * bar_coors[flag] \
                              + alpha * bar_centre[None,:]

        centres = transform_bar_to_space_coors(bar_coors, coors)

    return centres
