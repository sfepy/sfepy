"""
Finite element reference mappings.
"""
import numpy as nm

from sfepy import Config
from sfepy.base.base import get_default, output
from sfepy.base.mem_usage import raise_if_too_large
from sfepy.discrete.common.mappings import Mapping
from sfepy.discrete.common.extmods.mappings import CMapping
from sfepy.discrete import PolySpace


class FEMapping(Mapping):
    """
    Base class for finite element mappings.
    """

    def __init__(self, coors, conn, poly_space=None, gel=None, order=1):
        self.coors = coors
        self.conn = conn

        try:
            nm.take(self.coors, self.conn)

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
        self.indices = None

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
        if self.indices is not None:
            ii = max(self.dim - 1, 1)
            bf = nm.ascontiguousarray(bf[..., :ii:, self.indices])

        return bf

    def set_basis_indices(self, indices):
        """
        Set indices to cell-based basis that give the facet-based basis.
        """
        self.indices = indices

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

    def get_mapping(self, qp_coors, weights, poly_space=None, ori=None,
                    transform=None, mode='surface', is_face=False):
        """
        Get the mapping for given quadrature points, weights, and
        polynomial space.

        Returns
        -------
        cmap : CMapping instance
            The domain mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        bf_g = self.get_base(qp_coors, diff=True)
        if nm.allclose(bf_g, 0.0) and self.dim > 1:
            raise ValueError('zero base function gradient!')

        if not is_face:
            ebf_g = poly_space.eval_base(qp_coors, diff=True, ori=ori,
                                        force_axis=True, transform=transform)
            size = ebf_g.nbytes * self.n_el
            site_config = Config()
            raise_if_too_large(size, site_config.refmap_memory_factor())

            flag = (ori is not None) or (ebf_g.shape[0] > 1)
            mode = 'volume'
        else:
            flag = 0
            ebf_g = None

        cmap = CMapping(self.n_el, qp_coors.shape[0], self.dim,
                        poly_space.n_nod, mode=mode, flag=flag)
        cmap.describe(self.coors, self.conn, bf_g, ebf_g, weights)

        if self.dim == 1 and cmap.normal is not None:
            # Fix normals.
            ii = nm.where(self.conn == 0)[0]
            cmap.normal[ii] *= -1.0

        return cmap


class VolumeMapping(FEMapping):
    """
    Mapping from reference domain to physical domain of the same space
    dimension.
    """

    def get_mapping(self, qp_coors, weights, poly_space=None, ori=None,
                    transform=None):
        print('vol_map')
        return FEMapping.get_mapping(self, qp_coors, weights,
                                     poly_space=poly_space, ori=ori,
                                     transform=transform, is_face=False)


class SurfaceMapping(FEMapping):
    """
    Mapping from reference domain to physical domain of the space
    dimension higher by one.
    """

    def get_mapping(self, qp_coors, weights, poly_space=None, mode='surface'):
        print('surf_map')
        return FEMapping.get_mapping(self, qp_coors, weights,
                                     poly_space=poly_space, ori=None,
                                     transform=None, mode=mode, is_face=True)
