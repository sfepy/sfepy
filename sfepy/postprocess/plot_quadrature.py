"""
Functions to visualize quadrature points in reference elements.
"""
import numpy as nm
import matplotlib.pyplot as plt

from sfepy.base.base import output

from sfepy.postprocess.plot_dofs import _get_axes, _to2d
from sfepy.postprocess.plot_facets import plot_geometry

def _get_qp(geometry, order):
    from sfepy.discrete import Integral
    from sfepy.discrete.fem.geometry_element import GeometryElement

    aux = Integral('aux', order=order)
    coors, weights = aux.get_qp(geometry)
    true_order = aux.qps[geometry].order

    output('geometry:', geometry, 'order:', order, 'num. points:',
           coors.shape[0], 'true_order:', true_order)
    output('min. weight:', weights.min())
    output('max. weight:', weights.max())

    return GeometryElement(geometry), coors, weights

def _get_bqp(geometry, order):
    from sfepy.discrete import Integral
    from sfepy.discrete.fem.geometry_element import GeometryElement
    from sfepy.discrete.fem import Mesh, FEDomain, Field

    gel = GeometryElement(geometry)

    mesh = Mesh.from_data('aux', gel.coors, None,
                          [gel.conn[None, :]], [[0]], [geometry])
    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')
    surf = domain.create_region('Surf', 'vertices of surface', 'facet')
    field = Field.from_args('f', nm.float64, shape=1,
                            region=omega, approx_order=1)
    field.setup_surface_data(surf)

    integral = Integral('aux', order=order)
    field.create_bqp('Surf', integral)

    sd = field.extra_data['sd_Surf']
    qp = field.qp_coors[(integral.order, sd.bkey)]

    output('geometry:', geometry, 'order:', order, 'num. points:',
           qp.vals.shape[1], 'true_order:',
           integral.qps[gel.surface_facet_name].order)
    output('min. weight:', qp.weights.min())
    output('max. weight:', qp.weights.max())

    return (gel, qp.vals.reshape((-1, mesh.dim)),
            nm.tile(qp.weights, qp.vals.shape[0]))

def plot_weighted_points(ax, coors, weights, min_radius=10, max_radius=50,
                         show_colorbar=False):
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

    coors = _to2d(coors)
    sc = ax.scatter(*coors.T, s=nweights, c=weights, alpha=1)

    if show_colorbar:
        plt.colorbar(sc)

    return ax

def label_points(ax, coors):
    """
    Label points with their indices.
    """
    dim = coors.shape[1]
    ax = _get_axes(ax, dim)

    shift = 0.02 * (coors.max(0) - coors.min(0))
    ccs = coors + shift
    for ic, cc in enumerate(ccs):
        ax.text(*cc, s='%d' % ic, color='b')

def plot_quadrature(ax, geometry, order, boundary=False,
                    min_radius=10, max_radius=50,
                    show_colorbar=False, show_labels=False):
    """
    Plot quadrature points for the given geometry and integration order.

    The points are plotted as circles/spheres with radii given by quadrature
    weights - the weights are mapped to [`min_radius`, `max_radius`] interval.
    """
    if not boundary:
        gel, coors, weights = _get_qp(geometry, order)

    else:
        gel, coors, weights = _get_bqp(geometry, order)

    dim = coors.shape[1]

    ax = _get_axes(ax, dim)

    plot_geometry(ax, gel)
    plot_weighted_points(ax, coors, weights,
                         min_radius=min_radius, max_radius=max_radius,
                         show_colorbar=show_colorbar)
    if show_labels:
        label_points(ax, coors)

    return ax, coors, weights
