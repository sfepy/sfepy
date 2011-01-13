"""
Finite element reference mappings.
"""
import numpy as nm

from sfepy.base.base import output, get_default, Struct
from sfepy.fem.poly_spaces import PolySpace
import extmods.geometry as gm

def get_physical_qps(region, integral):
    """
    Get physical quadrature points corresponding to the given region
    and integral.
    """
    phys_qps = Struct(group_indx = {},
                      el_indx = {},
                      n_qp = {},
                      values = {},
                      is_uniform = True)

    ii = 0
    for ig in region.igs:
        gmap = region.create_mapping(integral.kind[0], ig)

        gel = gmap.get_geometry()
        qp_coors, _ = integral.get_qp(gel.name)

        qps = gmap.get_physical_qps(qp_coors)

        phys_qps.n_qp[ig] = n_qp = qps.shape[0] * qps.shape[1]

        phys_qps.group_indx[ig] = nm.arange(ii, ii + n_qp,
                                            dtype=nm.int32)

        aux = nm.tile(nm.array(qps.shape[1], dtype=nm.int32), n_qp + 1)
        aux[0] = 0
        phys_qps.el_indx[ig] = nm.cumsum(aux)

        ii += n_qp

        qps.shape = (n_qp, qps.shape[2])
        phys_qps.values[ig] = qps

    return phys_qps

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
        geo : VolumeGeometry instance
            The geometry object that describes the mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        geo = gm.VolumeGeometry(self.n_el, qp_coors.shape[0], self.dim,
                                poly_space.n_nod)
        geo.mode = gm.GM_Material

        bf_g = self.get_base(qp_coors, diff=True)
        ebf_g = poly_space.eval_base(qp_coors, diff=True)

        if nm.allclose(bf_g, 0.0):
            raise ValueError('zero base function gradient!')

        try:
            geo.describe(self.coors, self.conn, bf_g, ebf_g, weights)
        except:
            gm.errclear()
            raise

        return geo

        return geo

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
        geo : SurfaceGeometry instance
            The geometry object that describes the mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        geo = gm.SurfaceGeometry(self.n_el, qp_coors.shape[0], self.dim,
                                 poly_space.n_nod)
        geo.mode = gm.GM_Material

        bf_g = self.get_base(qp_coors, diff=True)

        if nm.allclose(bf_g, 0.0):
            raise ValueError('zero base function gradient!')

        try:
            geo.describe(self.coors, self.conn, bf_g, weights)
        except:
            gm.errclear()
            raise

        return geo
