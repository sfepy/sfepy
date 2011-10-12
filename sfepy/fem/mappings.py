"""
Finite element reference mappings.
"""
import numpy as nm

from sfepy.base.base import output, get_default, Struct
from sfepy.fem.poly_spaces import PolySpace
from sfepy.fem.extmods.mappings import CVolumeMapping, CSurfaceMapping

class PhysicalQPs(Struct):
    """
    Physical quadrature points in a region.
    """

    def get_merged_values(self):
        qps = nm.concatenate([self.values[ig] for ig in self.igs], axis=0)

        return qps

    def get_shape(self, rshape, ig):
        """
        Get shape from raveled shape.
        """
        n_qp = self.shape[ig][1]
        shape = (rshape[0] / n_qp, n_qp, rshape[1], rshape[2])

        return shape

def get_physical_qps(region, integral):
    """
    Get physical quadrature points corresponding to the given region
    and integral.
    """
    phys_qps = PhysicalQPs(igs=region.igs,
                           indx={}, rindx={}, qp_indx={},
                           n_per_group={}, shape={}, values={},
                           is_uniform=True)

    ii = 0
    for ig in region.igs:
        gmap = region.create_mapping(integral.kind[0], ig)

        gel = gmap.get_geometry()
        qp_coors, _ = integral.get_qp(gel.name)

        qps = gmap.get_physical_qps(qp_coors)

        n_el, n_qp = qps.shape[0], qps.shape[1]

        phys_qps.n_per_group[ig] = n_per_group = n_el * n_qp
        phys_qps.shape[ig] = qps.shape

        phys_qps.indx[ig] = slice(ii, ii + n_el)
        phys_qps.rindx[ig] = slice(ii * n_qp, (ii + n_el) * n_qp)

        aux = nm.tile(nm.array(qps.shape[1], dtype=nm.int32), n_per_group + 1)
        aux[0] = 0
        phys_qps.qp_indx[ig] = nm.cumsum(aux)

        ii += qps.shape[0]

        qps.shape = (n_per_group, qps.shape[2])
        phys_qps.values[ig] = qps

    return phys_qps

def get_jacobian(field, integral, integration='volume'):
    """
    Get the jacobian of reference mapping corresponding to `field` in
    quadrature points of the given `integral`.

    Returns
    -------
    jac : array
       The jacobian merged for all element groups.

    Notes
    -----
    Assumes the same element geometry in all element groups of the field!
    """
    jac = None
    for ig in field.igs:
        geo, _ = field.get_mapping(ig, field.region, integral, integration)
        _jac = geo.det
        if jac is None:
            jac = _jac

        else:
            jac = nm.concatenate((jac, _jac), axis=0)

    return jac

class Mapping(Struct):
    """
    Base class for mappings.
    """

    def __init__(self, coors, conn, poly_space=None, gel=None, order=1):
        self.coors = coors
        self.conn = conn

        try:
            self.coors[self.conn]

        except IndexError:
            output('coordinates shape: %s' % list(coors.shape))
            output('connectivity: min: %d, max: %d' % (conn.min(), conn.max()))
            msg = 'incompatible connectivity and coordinates (see above)'
            raise IndexError(msg)

        self.n_el, self.n_ep = conn.shape
        self.dim = self.coors.shape[1]

        if poly_space is None:
            poly_space = PolySpace.any_from_args(None, gel, order,
                                                 base='lagrange',
                                                 force_bubble=False)

        self.poly_space = poly_space

    def get_geometry(self):
        """
        Return reference element geometry as a GeometryElement instance.
        """
        return self.poly_space.geometry

    def get_base(self, coors, diff=False):
        """
        Get base functions or their gradient evaluated in given
        coordinates.
        """
        bf = self.poly_space.eval_base(coors, diff=diff)
        return bf

    def get_physical_qps(self, qp_coors):
        """
        Get physical quadrature points corresponding to given reference
        element quadrature points.

        Returns
        -------
        qps : array
            The physical quadrature points ordered element by element,
            i.e. with shape (n_el, n_qp, dim).
        """
        bf = self.get_base(qp_coors)
        qps = nm.dot(nm.atleast_2d(bf.squeeze()), self.coors[self.conn])
        # Reorder so that qps are really element by element.
        qps = nm.ascontiguousarray(nm.swapaxes(qps, 0, 1))

        return qps

class VolumeMapping(Mapping):
    """
    Mapping from reference domain to physical domain of the same space
    dimension.
    """

    def get_mapping(self, qp_coors, weights, poly_space=None):
        """
        Get the mapping for given quadrature points, weights, and
        polynomial space.

        Returns
        -------
        cmap : CVolumeMapping instance
            The volume mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        bf_g = self.get_base(qp_coors, diff=True)
        ebf_g = poly_space.eval_base(qp_coors, diff=True)


        cmap = CVolumeMapping(self.n_el, qp_coors.shape[0], self.dim,
                              poly_space.n_nod)
        cmap.describe(self.coors, self.conn, bf_g, ebf_g, weights)

        return cmap

class SurfaceMapping(Mapping):
    """
    Mapping from reference domain to physical domain of the space
    dimension higher by one.
    """

    def get_mapping(self, qp_coors, weights, poly_space=None):
        """
        Get the mapping for given quadrature points, weights, and
        polynomial space.

        Returns
        -------
        cmap : CSurfaceMapping instance
            The surface mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        bf_g = self.get_base(qp_coors, diff=True)

        if nm.allclose(bf_g, 0.0):
            raise ValueError('zero base function gradient!')

        cmap = CSurfaceMapping(self.n_el, qp_coors.shape[0], self.dim,
                               poly_space.n_nod)
        cmap.describe(self.coors, self.conn, bf_g, weights)

        return cmap
