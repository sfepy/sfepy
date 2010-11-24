import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import assert_, output
from sfepy.linalg.utils import norm_l2_along_axis as norm
from sfepy.linalg.utils import mini_newton

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

def barycentric_coors(coors, s_coors):
    """
    Get barycentric (area in 2D, volume in 3D) coordinates of points
    with coordinates `coors` w.r.t. the simplex given by `s_coors`.

    Returns
    -------
    bc : array
        The barycentric coordinates. Then reference element coordinates
        `xi = dot(bc.T, ref_coors)`.
    """
    n_v, dim = s_coors.shape
    n_c, dim2 = coors.shape
    assert_(dim == dim2)
    assert_(n_v == (dim + 1))

    mtx = nm.ones((n_v, n_v), nm.float64)
    mtx[0:n_v-1,:] = s_coors.T

    rhs = nm.empty((n_v,n_c), nm.float64)
    rhs[0:n_v-1,:] = coors.T
    rhs[n_v-1,:] = 1.0

    bc = nla.solve(mtx, rhs)

    return bc

def points_in_simplex(coors, s_coors, eps=1e-8):
    """
    Test if points with coordinates `coors` are in the simplex given by
    `s_coors`.
    """
    n_c, dim = coors.shape
    bc = barycentric_coors(coors, s_coors)
    flag = nm.ones((n_c,), dtype=nm.bool)
    for idim in xrange(dim + 1):
        flag &= nm.where((bc[idim,:] > -eps)
                         & (bc[idim,:] < (1.0 + eps)), True, False)
    return flag


def inverse_element_mapping(coors, e_coors, eval_base, ref_coors,
                            suppress_errors=False):
    """
    Given spatial element coordinates, find the inverse mapping for
    points with coordinats X = X(xi), i.e. xi = xi(X).

    Returns
    -------
    xi : array
        The reference element coordinates.
    """
    n_v, dim = e_coors.shape
    if coors.ndim == 2:
        n_c, dim2 = coors.shape
    else:
        n_c, dim2 = 1, coors.shape[0]

    assert_(dim == dim2)

    if n_v == (dim + 1): # Simplex.
        bc = barycentric_coors(coors, e_coors)
        xi = nm.dot(bc.T, ref_coors)

    else: # Tensor-product and other.
        def residual(xi):
            bf = eval_base(xi[nm.newaxis,:],
                           suppress_errors=suppress_errors).squeeze()
            res = coors - nm.dot(bf, e_coors)
            return res.squeeze()

        def matrix(xi):
            bfg = eval_base(xi[nm.newaxis,:], diff=True,
                            suppress_errors=suppress_errors).squeeze()
            mtx = - nm.dot(bfg, e_coors)
            return mtx

        xi0 = nm.array([0.0, 0.0, 0.0])
        xi = mini_newton(residual, xi0, matrix)

    return xi

def rotation_matrix2d(angle):
    """
    Construct a 2D (plane) rotation matrix corresponding to `angle`.
    """
    angle *= nm.pi / 180.0
    mtx = nm.array([[nm.cos(angle), -nm.sin(angle)],
                    [nm.sin(angle), nm.cos(angle)]], dtype=nm.float64)
    return mtx

def make_axis_rotation_matrix(direction, angle):
    r"""
    Create a rotation matrix :math:`\ull{R}` corresponding to the
    rotation around a general axis :math:`\ul{d}` by a specified angle
    :math:`\alpha`.

    .. math::
        \ull{R} = \ul{d}\ul{d}^T + \cos(\alpha) (I - \ul{d}\ul{d}^T) +
        \sin(\alpha) \skewop(\ul{d})

    Parameters
    ----------
    direction : array
        The rotation axis direction vector :math:`\ul{d}`.
    angle : float
        The rotation angle :math:`\alpha`.

    Returns
    -------
    mtx : array
        The rotation matrix :math:`\ull{R}`.
    """
    d = nm.array(direction, dtype=nm.float64)
    d /= nm.linalg.norm(d)

    eye = nm.eye(3, dtype=nm.float64)
    ddt = nm.outer(d, d)
    skew = nm.array([[    0,  d[2],  -d[1]],
                     [-d[2],     0,  d[0]],
                     [d[1], -d[0],    0]], dtype=nm.float64)

    mtx = ddt + nm.cos(angle) * (eye - ddt) + nm.sin(angle) * skew
    return mtx

def get_coors_in_cylinder(coors, centre, axis, radius, length, inside=True):
    """
    Return indices of coordinates inside or outside a cylinder given by
    centre, axis, radius and length.
    """
    vec = coors.T - centre

    drv = nm.cross(axis, vec, axisb=0)
    dr = nm.sqrt(nm.sum(drv * drv, 1))
    dl = nm.dot(axis, vec)

    if inside:
        out = nm.where((dl >= 0.0) & (dl <= length) & (dr <= radius))[0]
    else:
        out = nm.where((dl >= 0.0) & (dl <= length) & (dr >= radius))[0]

    return out
