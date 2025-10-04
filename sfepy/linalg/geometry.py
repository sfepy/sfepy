import numpy as nm
import numpy.linalg as nla

from scipy.special import factorial

from sfepy.base.base import assert_, output
from sfepy.linalg.utils import norm_l2_along_axis as norm
from sfepy.linalg.utils import mini_newton, dets_fast

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

def get_simplex_volumes(cells, coors):
    """
    Get volumes of simplices in nD.

    Parameters
    ----------
    cells : array, shape (n, d)
        The indices of `n` simplices with `d` vertices into `coors`.
    coors : array
        The coordinates of simplex vertices.

    Returns
    -------
    volumes : array
        The volumes of the simplices.
    """
    scoors = coors[cells]
    deltas = scoors[:, 1:] - scoors[:, :1]
    dim = coors.shape[1]
    volumes = dets_fast(deltas) / factorial(dim)

    return volumes

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
    flag = nm.ones((n_c,), dtype=bool)
    for idim in range(dim + 1):
        flag &= nm.where((bc[idim,:] > -eps)
                         & (bc[idim,:] < (1.0 + eps)), True, False)
    return flag

def flag_points_in_polygon2d(polygon, coors):
    """
    Test if points are in a 2D polygon.

    Parameters
    ----------
    polygon : array, (:, 2)
        The polygon coordinates.
    coors: array, (:, 2)
        The coordinates of points.

    Returns
    -------
    flag : bool array
        The flag that is True for points that are in the polygon.

    Notes
    -----
    This is a semi-vectorized version of [1].

    [1] PNPOLY - Point Inclusion in Polygon Test, W. Randolph Franklin (WRF)
    """
    flag = nm.zeros(coors.shape[0], dtype=bool)
    nv = polygon.shape[0]
    px, py = coors[:, 0], coors[:, 1]
    for ii in range(nv):
        vix, viy = polygon[ii, 0], polygon[ii, 1]
        vjx, vjy = polygon[ii-1, 0], polygon[ii-1, 1]
        aux = nm.where((viy > py) != (vjy > py))
        flag[aux] = nm.where((px[aux] < (vjx - vix)
                              * (py[aux] - viy) / (vjy - viy) + vix),
                             ~flag[aux], flag[aux])
    return flag

def inverse_element_mapping(coors, e_coors, eval_basis, ref_coors,
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
            bf = eval_basis(xi[nm.newaxis,:].copy(),
                            suppress_errors=suppress_errors).squeeze()
            res = coors - nm.dot(bf, e_coors)
            return res.squeeze()

        def matrix(xi):
            bfg = eval_basis(xi[nm.newaxis,:].copy(), diff=True,
                             suppress_errors=suppress_errors).squeeze()
            mtx = - nm.dot(bfg, e_coors)
            return mtx

        xi0 = nm.zeros(dim, dtype=nm.float64)
        xi = mini_newton(residual, xi0, matrix)

    return xi

def get_perpendiculars(vec):
    """
    For a given vector, get a unit vector perpendicular to it in 2D, or get two
    mutually perpendicular unit vectors perpendicular to it in 3D.
    """
    nvec = nm.linalg.norm(vec)
    vec /= nvec

    if vec.shape[0] == 2:
        out = nm.array([vec[1], -vec[0]], dtype=nm.float64)

    else:
        aux = nm.array([0.0, 0.0, 1.0], dtype=nm.float64)

        v1 = nm.cross(vec, aux)
        if nm.linalg.norm(v1) < 0.1:
            # vec and aux close to being co-linear.
            aux = nm.array([0.0, 1.0, 0.0], dtype=nm.float64)

            v1 = nm.cross(vec, aux)

        v1 /= nm.linalg.norm(v1)

        v2 = nm.cross(vec, v1)
        v2 /= nm.linalg.norm(v2)

        out = (v1, v2)

    return out

def get_face_areas(faces, coors):
    """
    Get areas of planar convex faces in 2D and 3D.

    Parameters
    ----------
    faces : array, shape (n, m)
        The indices of `n` faces with `m` vertices into `coors`.
    coors : array
        The coordinates of face vertices.

    Returns
    -------
    areas : array
        The areas of the faces.
    """
    faces = nm.asarray(faces)
    coors = nm.asarray(coors)

    n_v = faces.shape[1]

    if n_v == 3:
        aux = coors[faces]
        if aux.shape[-1] == 2:
            zz = nm.zeros(aux.shape[:-1] + (1,))
            aux = nm.concatenate((aux, zz), axis=-1)

        v1 = aux[:, 1, :] - aux[:, 0, :]
        v2 = aux[:, 2, :] - aux[:, 0, :]
        areas = 0.5 * norm(nm.cross(v1, v2))

    elif n_v == 4:
        areas1 = get_face_areas(faces[:, [0, 1, 2]], coors)
        areas2 = get_face_areas(faces[:, [0, 2, 3]], coors)

        areas = areas1 + areas2

    else:
        raise ValueError('unsupported faces! (%d vertices)' % n_v)

    return areas

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

    Notes
    -----
    The matrix follows the right hand rule: if the right hand thumb
    points along the axis vector :math:`\ul{d}` the fingers show the
    positive angle rotation direction.

    Examples
    --------
    Make transformation matrix for rotation of coordinate system by 90
    degrees around 'z' axis.

    >>> mtx = make_axis_rotation_matrix([0., 0., 1.], nm.pi/2)
    >>> mtx
    array([[ 0.,  1.,  0.],
           [-1.,  0.,  0.],
           [ 0.,  0.,  1.]])

    Coordinates of vector :math:`[1, 0, 0]^T` w.r.t. the original system
    in the rotated system. (Or rotation of the vector by -90 degrees in
    the original system.)

    >>> nm.dot(mtx, [1., 0., 0.])
    >>> array([ 0., -1.,  0.])

    Coordinates of vector :math:`[1, 0, 0]^T` w.r.t. the rotated system
    in the original system. (Or rotation of the vector by +90 degrees in
    the original system.)

    >>> nm.dot(mtx.T, [1., 0., 0.])
    >>> array([ 0.,  1.,  0.])
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

def get_coors_in_tube(coors, centre, axis, radius_in, radius_out, length,
                      inside_radii=True):
    """
    Return indices of coordinates inside a tube given by
    centre, axis vector, inner and outer radii and length.

    Parameters
    ----------
    inside_radii : bool, optional
        If False, select points outside the radii, but within the tube
        length.

    Notes
    -----
    All float comparisons are done using `<=` or `>=` operators,
    i.e. the points on the boundaries are taken into account.
    """
    coors = nm.asarray(coors)
    centre = nm.asarray(centre)

    vec = coors - centre[None, :]

    drv = nm.cross(axis, vec, axisb=1)
    dr = nm.sqrt(nm.sum(drv * drv, 1))
    dl = nm.dot(vec, axis)

    l2 = 0.5 * length

    if inside_radii:
        out = nm.where((dl >= -l2) & (dl <= l2) &
                       (dr >= radius_in) & (dr <= radius_out))[0]

    else:
        out = nm.where((dl >= -l2) & (dl <= l2) &
                       (dr <= radius_in) & (dr >= radius_out))[0]

    return out

def get_coors_in_ball(coors, centre, radius, radius2=None, inside=True):
    """
    Return indices of coordinates inside or outside a ball given by
    centre and radius::

      inside radius radius2 condition
      True   r      None          |x - c| <= r
      True   r      r2      r2 <= |x - c| <= r
      False  r      None    |x - c| >= r
      False  r      r2      |x - c| >= r | |x - c| <= r2

    Notes
    -----
    All float comparisons are done using `<=` or `>=` operators,
    i.e. the points on the boundaries are taken into account.
    """
    coors = nm.asarray(coors)
    centre = nm.asarray(centre)

    vec = coors - centre[None, :]
    nvec = norm(vec)
    if radius2 is None:
        if inside:
            out = nm.where(nvec <= radius)[0]

        else:
            out = nm.where(nvec >= radius)[0]

    else:
        if inside:
            out = nm.where((nvec <= radius) & (nvec >= radius2))[0]

        else:
            out = nm.where((nvec >= radius) & (nvec <= radius2))[0]

    return out
