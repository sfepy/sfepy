"""
Functions to visualize quadrature points in reference elements.
"""
import matplotlib.pyplot as plt

from sfepy.base.base import output

from sfepy.postprocess.plot_dofs import _get_axes
from sfepy.postprocess.plot_facets import plot_geometry

def _get_qp(geometry, order):
    from sfepy.fem import Integral
    from sfepy.fem.geometry_element import GeometryElement

    aux = Integral('aux', order=order)
    coors, weights = aux.get_qp(geometry)
    true_order = aux.qps[geometry].order

    output('geometry:', geometry, 'order:', order, 'num. points:',
           coors.shape[0], 'true_order:', true_order)
    output('min. weight:', weights.min())
    output('max. weight:', weights.max())

    return GeometryElement(geometry), coors, weights

def plot_weighted_points(ax, coors, weights, min_radius=10, max_radius=50,
                         show_colorbar=False, show=False):
    """
    Plot points with given coordinates as circles/spheres with radii given by
    weights.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    wmin, wmax = weights.min(), weights.max()
    if (wmax - wmin) < 1e-12:
        nweights = weights * max_radius / wmax

    else:
        nweights = ((weights - wmin) * (max_radius - min_radius)
                    / (wmax - wmin) + min_radius)

    if dim == 3:
        sc = ax.scatter(coors[:, 0], coors[:, 1], coors[:, 2],
                        s=nweights, c=weights, alpha=1)

    else:
        sc = ax.scatter(coors[:, 0], coors[:, 1],
                        s=nweights, c=weights, alpha=1)

    if show_colorbar:
        plt.colorbar(sc)

    if show:
        plt.show()

    return ax

def plot_quadrature(ax, geometry, order, min_radius=10, max_radius=50,
                    show_colorbar=False, show=False):
    """
    Plot quadrature points for the given geometry and integration order.

    The points are plotted as circles/spheres with radii given by quadrature
    weights - the weights are mapped to [`min_radius`, `max_radius`] interval.
    """
    gel, coors, weights = _get_qp(geometry, order)
    dim = coors.shape[1]

    ax = _get_axes(ax, dim)

    plot_geometry(ax, gel, show=False)
    plot_weighted_points(ax, coors, weights,
                         min_radius=min_radius, max_radius=max_radius,
                         show_colorbar=show_colorbar, show=show)

    return ax
