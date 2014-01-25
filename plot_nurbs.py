import numpy as nm
import matplotlib.pyplot as plt

from sfepy.fem.geometry_element import GeometryElement
from sfepy.mesh.mesh_generators import get_tensor_product_conn
import sfepy.postprocess.plot_dofs as pd
from sfepy.postprocess.plot_dofs import _get_axes

def _get_knots_tuple(knots):
    if isinstance(knots, nm.ndarray) and (knots.ndim == 1):
        knots = (knots,)

    elif not isinstance(knots, tuple):
        raise ValueError('knots must be 1D array or a tuple of 1D arrays!')

    return knots

def plot_parametric_mesh(ax, knots, show=False):
    """
    Plot the parametric mesh of a NURBS given by its knots.
    """
    knots = _get_knots_tuple(knots)
    dim = len(knots)

    ax = _get_axes(ax, dim)

    uknots = [nm.unique(ii) for ii in knots]
    shape = [len(ii) for ii in uknots]

    ngrid = nm.mgrid[[slice(ii) for ii in shape]]
    coors = nm.r_[[uknots[ii][ig].ravel() for ii, ig in enumerate(ngrid)]].T

    conn, desc = get_tensor_product_conn(nm.array(shape))
    gel = GeometryElement(desc)

    ax = pd.plot_mesh(ax, coors, conn, gel.edges)
    pd.plot_points(ax, coors)

    if show:
        plt.show()

    return ax

def plot_control_mesh(ax, control_points, show=False):
    """
    Plot the control mesh of a NURBS given by its control points.
    """
    dim = control_points.shape[-1]
    ax = _get_axes(ax, dim)

    shape = control_points.shape

    coors = control_points.copy()
    coors.shape = (nm.prod(shape[:-1]), shape[-1])

    conn, desc = get_tensor_product_conn(nm.array(shape[:-1]))
    gel = GeometryElement(desc)

    ax = pd.plot_mesh(ax, coors, conn, gel.edges)
    pd.plot_points(ax, coors)

    if show:
        plt.show()

    return ax

def plot_iso_lines(ax, nurbs, color='b', n_points=100, show=False):
    """
    Plot the NURBS <object using iso-lines in Greville abscissae coordinates.
    """
    dim = nurbs.dim
    ax = _get_axes(ax, dim)

    gas = nurbs.greville()

    if dim == 1:
        ga = gas[0]
        line = nm.linspace(ga[0], ga[-1], n_points)

        vals = nurbs(line)[:, 0]
        ax.plot(line, vals, color)

    elif dim == 2:
        ga0 = gas[0]
        ga1 = gas[1]

        x1 = nm.linspace(ga1[0], ga1[-1], n_points)
        for x0 in ga0:
            vals = nurbs(x0, x1)
            ax.plot(vals[:, 0], vals[:, 1], color)

        x0 = nm.linspace(ga0[0], ga0[-1], n_points)
        for x1 in ga0:
            vals = nurbs(x0, x1)
            ax.plot(vals[:, 0], vals[:, 1], color)

    else:
        ga0 = gas[0]
        ga1 = gas[1]
        ga2 = gas[2]

        x2 = nm.linspace(ga2[0], ga2[-1], n_points)
        for x0 in ga0:
            for x1 in ga1:
                vals = nurbs(x0, x1, x2)
                ax.plot(vals[:, 0], vals[:, 1], vals[:, 2], color)

        x1 = nm.linspace(ga1[0], ga1[-1], n_points)
        for x0 in ga0:
            for x2 in ga2:
                vals = nurbs(x0, x1, x2)
                ax.plot(vals[:, 0], vals[:, 1], vals[:, 2], color)

        x0 = nm.linspace(ga0[0], ga0[-1], n_points)
        for x1 in ga1:
            for x2 in ga2:
                vals = nurbs(x0, x1, x2)
                ax.plot(vals[:, 0], vals[:, 1], vals[:, 2], color)

    if show:
        plt.show()

    return ax
