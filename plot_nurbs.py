import numpy as nm
import matplotlib.pyplot as plt

from sfepy.fem.geometry_element import GeometryElement
from sfepy.mesh.mesh_generators import get_tensor_product_conn
import sfepy.postprocess.plot_dofs as pd
from sfepy.postprocess.plot_dofs import _get_axes

from igalib import _get_knots_tuple

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

def _get_edges(n_ep, shape):
    dim = len(shape)
    aux = nm.arange(n_ep).reshape(shape)

    edges = []
    if dim == 3:
        for ii in xrange(shape[2] - 1):
            edges.append(aux[0, 0, ii:ii+2])
            edges.append(aux[-1, 0, ii:ii+2])
            edges.append(aux[0, -1, ii:ii+2])
            edges.append(aux[-1, -1, ii:ii+2])

        for ii in xrange(shape[1] - 1):
            edges.append(aux[0, ii:ii+2, 0])
            edges.append(aux[-1, ii:ii+2, 0])
            edges.append(aux[0, ii:ii+2, -1])
            edges.append(aux[-1, ii:ii+2, -1])

        for ii in xrange(shape[0] - 1):
            edges.append(aux[ii:ii+2, 0, 0])
            edges.append(aux[ii:ii+2, -1, 0])
            edges.append(aux[ii:ii+2, 0, -1])
            edges.append(aux[ii:ii+2, -1, -1])

    elif dim == 2:
        for ii in xrange(shape[1] - 1):
            edges.append(aux[0, ii:ii+2])
            edges.append(aux[-1, ii:ii+2])

        for ii in xrange(shape[0] - 1):
            edges.append(aux[ii:ii+2, 0])
            edges.append(aux[ii:ii+2, -1])

    else:
        for ii in xrange(shape[0] - 1):
            edges.append(aux[ii:ii+2])

    return nm.array(edges)

def plot_bezier_mesh(ax, control_points, conn, degrees, show=False):
    """
    Plot the Bezier mesh of a NURBS given by its control points and
    connectivity.
    """
    dim = control_points.shape[-1]
    ax = _get_axes(ax, dim)

    edges = _get_edges(conn.shape[1], nm.asarray(degrees) + 1)
    ax = pd.plot_mesh(ax, control_points, conn, edges)
    pd.plot_points(ax, control_points)

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

def plot_nurbs_basis_1d(ax, nurbs, n_points=100, show=False):
    """
    Plot a 1D NURBS basis.
    """
    dim = nurbs.dim
    ax = _get_axes(ax, dim)

    ga = nurbs.greville()[0]

    n_fun = nurbs.weights.shape[0]
    line = nm.linspace(ga[0], ga[-1], n_points)
    for ii in xrange(n_fun):
        field = nm.zeros(n_fun)
        field[ii] = 1.0

        vals = nurbs.evaluate(fields=field, u=line)
        plt.plot(line, vals)

    if show:
        plt.show()

    return ax

def plot_bezier_nurbs_basis_1d(ax, control_points, weights, degrees, cs, conn,
                               n_points=20, show=False):
    """
    Plot a 1D NURBS basis using the Bezier extraction and local Bernstein
    basis.
    """
    from igalib import eval_variable_in_qp
    ax = _get_axes(ax, 2)

    n_fun = weights.shape[0]
    line = nm.linspace(0, 1, n_points)[:, None]
    for ii in xrange(n_fun):
        variable = nm.zeros((n_fun, 1))
        variable[ii] = 1.0

        coors, vals, dets = eval_variable_in_qp(variable, line,
                                                control_points,
                                                weights, degrees,
                                                cs, conn)
        plt.plot(coors[:, 0], vals)

    if show:
        plt.show()

    return ax
