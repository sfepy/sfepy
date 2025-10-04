"""
Functions to visualize the geometry elements and numbering and orientation of
their facets (edges and faces).

The standard geometry elements can be plotted by running::

  $ python sfepy/postprocess/plot_facets.py
"""
import numpy as nm
import matplotlib.pyplot as plt

from sfepy.linalg import (get_perpendiculars, normalize_vectors,
                          make_axis_rotation_matrix)
from sfepy.postprocess.plot_dofs import _get_axes, plot_mesh, plot_global_dofs

def plot_geometry(ax, gel):
    """
    Plot a geometry element as a wireframe.
    """
    ax = plot_mesh(ax, gel.coors, [gel.conn], gel.edges)
    ax = plot_global_dofs(ax, gel.coors, [gel.conn])

    return ax

def plot_edges(ax, gel, length):
    """
    Plot edges of a geometry element as numbered arrows.
    """
    dim = gel.dim
    ax = _get_axes(ax, dim)

    if gel.edges is None: return ax

    l2 = 0.5 * length
    for ii, edge in enumerate(gel.edges):
        cc = gel.coors[edge]
        centre = 0.5 * cc.sum(axis=0)

        vdir = (cc - centre)
        normalize_vectors(vdir)

        cc = l2 * vdir + centre
        draw_arrow(ax, cc, length=0.3*length, linewidth=3, color='b')

        ax.text(*centre, s=ii,
                color='b', fontsize=10, weight='light')

    return ax

def plot_faces(ax, gel, radius, n_point):
    """
    Plot faces of a 3D geometry element as numbered oriented arcs. An arc
    centre corresponds to the first node of a face. It points from the first
    edge towards the last edge of the face.
    """
    dim = gel.dim
    ax = _get_axes(ax, dim)

    if dim < 3: return ax

    for ii, face in enumerate(gel.faces):
        cc = gel.coors[face]

        t1 = cc[1, :] - cc[0, :]
        t2 = cc[-1, :] - cc[0, :]
        n = nm.cross(t1, t2)

        nt1 = nm.linalg.norm(t1)
        nt2 = nm.linalg.norm(t2)
        angle = nm.arccos(nm.dot(t1, t2) / (nt1 * nt2))

        da = angle / (n_point - 1)

        mtx = make_axis_rotation_matrix(n, da)

        rt = cc[0] + radius * t1 / nt1
        coors = [rt]
        for ip in range(n_point - 1):
            rt = nm.dot(mtx.T, (rt - cc[0])) + cc[0]
            coors.append(rt)

        coors = nm.array(coors, dtype=nm.float64)
        centre = coors.sum(axis=0) / coors.shape[0]

        draw_arrow(ax, coors, length=0.3*radius, linewidth=3, color='r')

        ax.text(*centre, s=ii,
                color='r', fontsize=10, weight='light')

    return ax

def draw_arrow(ax, coors, angle=20.0, length=0.3, **kwargs):
    """
    Draw a line ended with an arrow head, in 2D or 3D.
    """
    color = kwargs.get('color', 'b')

    c0 = coors[-2]
    c1 = coors[-1]

    vd = c1 - c0
    nvd = nm.linalg.norm(vd)
    vd /= nvd

    c0 = c1 - length * vd

    ps = get_perpendiculars(vd)

    rangle = nm.deg2rad(min(angle, 60.0))
    plength = length * nm.arctan(rangle)

    if coors.shape[1] == 2:
        from matplotlib.patches import Polygon

        cx, cy = coors[:, 0], coors[:, 1]

        ax.plot(cx, cy, **kwargs)

        p0 = c0 + plength * ps
        p1 = c0 - plength * ps

        pol = Polygon([p0, p1, c1], color=color)
        ax.add_artist(pol)

    else:
        import mpl_toolkits.mplot3d as plt3

        cx, cy, cz = coors[:, 0], coors[:, 1], coors[:, 2]

        ax.plot(cx, cy, cz, **kwargs)

        p00 = c0 + plength * ps[0]
        p01 = c0 - plength * ps[0]

        p10 = c0 + plength * ps[1]
        p11 = c0 - plength * ps[1]

        arr = plt3.art3d.Poly3DCollection([[p00, p01, c1],
                                           [p10, p11, c1]], color=color)
        ax.add_collection3d(arr)

if __name__ == '__main__':
    from sfepy.discrete.fem.geometry_element import (GeometryElement,
                                                     geometry_data)

    for key, gd in geometry_data.items():
        if key == '1_2' : continue

        gel = GeometryElement(key)

        ax = plot_geometry(None, gel)
        ax = plot_edges(ax, gel, length=0.2)
        ax = plot_faces(ax, gel, radius=0.3, n_point=5)

        dd = 0.05
        ax.set_xlim([-dd, 1.0 + dd])
        ax.set_ylim([-dd, 1.0 + dd])
        if gel.dim == 3:
            ax.set_zlim([-dd, 1.0 + dd])

        plt.show()
