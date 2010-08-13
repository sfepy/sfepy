"""
Finite element reference mappings.
"""
from sfepy.base.base import *

from sfepy.fem.poly_spaces import PolySpace
import extmods.geometry as gm

class Mapping(Struct):
    """
    Base class for mappings.
    """

    def __init__(self, coors, conn, poly_space=None, gel=None, order=1):
        self.coors = coors
        self.conn = conn

        self.n_el, self.n_ep = conn.shape
        self.dim = self.coors.shape[1]

        if poly_space is None:
            poly_space = PolySpace.any_from_args(None, gel, order,
                                                 base='lagrange',
                                                 force_bubble=False)

        self.poly_space = poly_space

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

    def _describe(self, geo, qp_coors, weights):
        bf_g = self.get_base(qp_coors, diff=True)

        if nm.allclose(bf_g, 0.0):
            raise ValueError('zero base function gradient!')

        try:
            geo.describe(self.coors, self.conn, bf_g, weights)
        except:
            gm.errclear()
            raise

class VolumeMapping(Mapping):
    """
    Mapping from reference domain to physical domain of the same space
    dimension.
    """

    def get_mapping(self, qp_coors, weights):
        """
        Get the mapping for given quadrature points and weights.

        Returns
        -------
        geo : VolumeGeometry instance
            The geometry object that describes the mapping.
        """
        geo = gm.VolumeGeometry(self.n_el, qp_coors.shape[0], self.dim,
                                self.n_ep)
        geo.mode = gm.GM_Material
        self._describe(geo, qp_coors, weights)

        return geo

class SurfaceMapping(Mapping):
    """
    Mapping from reference domain to physical domain of the space
    dimension higher by one.
    """

    def get_mapping(self, qp_coors, weights):
        """
        Get the mapping for given quadrature points and weights.

        Returns
        -------
        geo : SurfaceGeometry instance
            The geometry object that describes the mapping.
        """
        geo = gm.SurfaceGeometry(self.n_el, qp_coors.shape[0], self.dim,
                                 self.n_ep)
        geo.mode = gm.GM_Material
        self._describe(geo, qp_coors, weights)

        return geo
