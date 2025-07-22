import numpy as nm

from sfepy.base.base import assert_, Struct
import sfepy.linalg as la

class ContactPlane(Struct):

    def __init__(self, anchor, normal, bounds, eps=1e-14):
        Struct.__init__(self, anchor=nm.array(anchor, dtype=nm.float64),
                        bounds=nm.asarray(bounds, dtype=nm.float64),
                        eps=eps)
        self.normal = nm.asarray(normal, dtype=nm.float64)

        norm = nm.linalg.norm
        self.normal /= norm(self.normal)

        e3 = [0.0, 0.0, 1.0]
        dd = nm.dot(e3, self.normal)
        rot_angle = nm.arccos(dd)

        if nm.abs(rot_angle) < eps:
            mtx = nm.eye(3, dtype=nm.float64)
            bounds2d = self.bounds[:, :2]

        else:
            if nm.abs(rot_angle - nm.pi) < eps:
                rot_axis = [1.0, 0.0, 0.0]

            else:
                rot_axis = nm.cross([0.0, 0.0, 1.0], self.normal)

            mtx = la.make_axis_rotation_matrix(rot_axis, rot_angle)

            mm = la.insert_strided_axis(mtx, 0, self.bounds.shape[0])
            rbounds = la.dot_sequences(mm, self.bounds)
            bounds2d = rbounds[:, :2]

        assert_(nm.allclose(nm.dot(mtx, self.normal), e3,
                            rtol=0.0, atol=1e-12))

        self.adotn = nm.dot(self.anchor, self.normal)

        self.rot_angle = rot_angle
        self.mtx = mtx
        self.bounds2d = bounds2d

    def mask_points(self, points):
        mm = la.insert_strided_axis(self.mtx, 0, points.shape[0])
        points2d = la.dot_sequences(mm, points)[:, :2]

        return la.flag_points_in_polygon2d(self.bounds2d, points2d)

    def get_distance(self, points):
        dist = la.dot_sequences(points, self.normal) - self.adotn

        return dist

class ContactSphere(Struct):

    def __init__(self, centre, radius):
        self.centre = nm.asarray(centre)
        self.radius = radius

    def mask_points(self, points, eps):
        dist2 = la.norm_l2_along_axis(points - self.centre, squared=True)
        radius2 = self.radius**2
        mask = dist2 <= ((1 + eps)**2) * radius2
        return mask

    def get_distance(self, points):
        """
        Get the penetration distance and normals of points w.r.t. the sphere
        surface.

        Returns
        -------
        d : array
            The penetration distance.
        normals : array
            The normals from the points to the sphere centre.
        """
        vecs = self.centre - points
        dist = la.norm_l2_along_axis(vecs)
        # Prevent zero division.
        ii = dist > 1e-8
        normals = nm.where(ii[:, None], vecs[ii] / dist[ii][:, None],
                           vecs[ii])
        return self.radius - dist, normals

    def _get_derivatives(self, points):
        vecs = self.centre - points
        dist = la.norm_l2_along_axis(vecs)

        # Distance derivative w.r.t. point coordinates.
        dd = vecs / dist[:, None]

        normals = dd
        # Unit normal derivative w.r.t. point coordinates.
        dim = points.shape[1]
        ee = nm.eye(dim)[None, ...]
        nnt = normals[..., None] * normals[..., None, :]
        dn = - (ee - nnt) / dist[:, None, None]

        return dd, dn

def plot_polygon(ax, polygon):
    from sfepy.postprocess.plot_dofs import _get_axes

    dim = polygon.shape[1]
    ax = _get_axes(ax, dim)

    pp = nm.r_[polygon, polygon[:1]]
    px, py = pp[:, 0], pp[:, 1]
    if dim == 2:
        ax.plot(px, py)

    else:
        pz = pp[:, 2]
        ax.plot(px, py, pz)

    return ax

def plot_points(ax, points, marker, **kwargs):
    from sfepy.postprocess.plot_dofs import _get_axes

    dim = points.shape[1]
    ax = _get_axes(ax, dim)

    px, py = points[:, 0], points[:, 1]
    if dim == 2:
        ax.plot(px, py, marker, **kwargs)

    else:
        pz = points[:, 2]
        ax.plot(px, py, pz, marker, **kwargs)

    return ax

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # Test and draw the plane.
    anchor = [1, 1, 1]
    normal = [2, -1, 1]
    bounds = [[-2, 0, 0],
              [2, 1, 0],
              [4, 3, 1],
              [1, 3, 1],
              [2, 2, 1]]
    cp = ContactPlane(anchor, normal, bounds)

    pps = 2 * nm.random.rand(20, 3)
    mask = cp.mask_points(pps)

    dist = cp.get_distance(pps)

    v1, v2 = la.get_perpendiculars(cp.normal)

    ax = plot_polygon(None, cp.bounds)
    ax = plot_polygon(ax, nm.r_[cp.anchor[None, :],
                                cp.anchor[None, :] + cp.normal[None, :]])
    ax = plot_polygon(ax, nm.r_[cp.anchor[None, :],
                                cp.anchor[None, :] + v1])
    ax = plot_polygon(ax, nm.r_[cp.anchor[None, :],
                                cp.anchor[None, :] + v2])
    ax = plot_points(ax, cp.anchor[None, :], 'r*')
    ax = plot_points(ax, pps[mask], 'bs', ms=10, mec='None')
    ax = plot_points(ax, pps[~mask], 'go', ms=10, mec='None')

    mask = dist >= 0.0
    ax = plot_points(ax, pps[mask], 'r^', mec='None')
    ax = plot_points(ax, pps[~mask], 'kv', mec='None')

    # Test and draw the sphere.
    pps = nm.random.rand(5000, 3)

    centre = [0, 0.5, 0.5]
    radius = 0.8
    cs = ContactSphere(centre, radius)
    mask = cs.mask_points(pps, 0.0)
    dist = cs.get_distance(pps)

    ax = plot_points(None, cs.centre[None, :], 'b*', ms=30)
    ax = plot_points(ax, pps[mask], 'kv')
    ax = plot_points(ax, pps[~mask], 'r.')

    plt.show()
