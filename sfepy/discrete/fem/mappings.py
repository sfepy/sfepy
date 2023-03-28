"""
Finite element reference mappings.
"""
import numpy as nm

from sfepy import Config
from sfepy.base.base import get_default, output
from sfepy.base.mem_usage import raise_if_too_large
from sfepy.discrete.common.mappings import Mapping, PyCMapping
from sfepy.discrete import PolySpace
from sfepy.linalg.utils import invs_fast, dets_fast


def eval_mapping_data_in_qp(coors, conn, dim, bf_g, weights,
                            ebf_g=None, is_face=False, eps=1e-15,
                            se_conn=None, se_bf_bg=None):
    """
    Evaluate mapping data.

    Parameters
    ----------
    coors: numpy.ndarray
        The nodal coordinates.
    conn: numpy.ndarray
        The element connectivity.
    dim: int
        The space dimension.
    bf_g: numpy.ndarray
        The derivatives of the domain basis functions with respect to the
        reference coordinates.
    weights: numpy.ndarray
        The weights of the quadrature points.
    ebf_g: numpy.ndarray
        The derivatives of the field basis functions with respect to the
        reference coordinates.
    is_face: bool
        Is it the boundary of a region?
    eps: float
        The tolerance for the normal vectors calculation.
    se_conn: numpy.ndarray
        The connectivity for the calculation of surface derivatives.
    se_bf_bg: numpy.ndarray
        The surface basis function derivatives with respect to the reference
        coordinates.

    Returns
    -------
    det: numpy.ndarray
        The determinant of the mapping evaluated in integration points.
    volume: numpy.ndarray
        The element (volume or surface) volumes in integration points.
    bfg: numpy.ndarray
        The derivatives of the basis functions with respect to the spatial
        coordinates. Can be evaluated either for surface elements if `bf_g`,
        `se_conn`, and `se_bf_bg` are given.
    normal: numpy.ndarray
        The normal vectors for the surface elements in integration points.
    """
    mtxRM = nm.einsum('qij,cjk->cqik', bf_g, coors[conn, :dim], optimize=True)

    n_el, n_qp = mtxRM.shape[:2]

    if is_face:
        # outward unit normal vector
        normal = nm.ones((n_el, n_qp, dim, 1), dtype=nm.float64)
        if dim == 1:
            det = nm.tile(weights, (n_el, 1)).reshape(n_el, n_qp, 1, 1)
            ii = nm.where(conn == 0)[0]
            normal[ii] *= -1.0
        elif dim == 2:
            c1, c2 = mtxRM[..., 0], mtxRM[..., 1]
            det0 = nm.sqrt(c1**2 + c2**2).reshape(n_el, n_qp, 1, 1)
            det = det0 * weights[:, None, None]
            det0[nm.abs(det0) < eps] = 1.0
            normal[..., 0, :] = c2
            normal[..., 1, :] = -c1
            normal /= det0
        elif dim == 3:
            j012 = mtxRM[..., 0, :].T
            j345 = mtxRM[..., 1, :].T
            c1 = (j012[1] * j345[2] - j345[1] * j012[2]).T
            c2 = (j012[0] * j345[2] - j345[0] * j012[2]).T
            c3 = (j012[0] * j345[1] - j345[0] * j012[1]).T
            det0 = nm.sqrt(c1**2 + c2**2 + c3**2).reshape(n_el, n_qp, 1, 1)
            det = det0 * weights[:, None, None]
            det0[nm.abs(det0) < eps] = 1.0
            normal[..., 0, 0] = c1
            normal[..., 1, 0] = -c2
            normal[..., 2, 0] = c3
            normal /= det0
    else:
        # det0 = nm.linalg.det(mtxRM)
        det0 = dets_fast(mtxRM)
        if nm.any(det0 <= 0.0):
            raise ValueError('warp violation!')

        det = det0 * weights
        det = det.reshape(n_el, n_qp, 1, 1)
        normal = None

    if ebf_g is not None:
        if is_face and se_conn is not None and se_bf_bg is not None:
            mtxRM = nm.einsum('cqij,cjk->cqik', se_bf_bg, coors[se_conn, :dim],
                              optimize=True)
            mtxRMI = invs_fast(mtxRM)
            bfg = nm.einsum('cqij,cqjk->cqik', mtxRMI, ebf_g, optimize=True)
        else:
            mtxRMI = invs_fast(mtxRM, det0)
            es_arg2 = 'x' if ebf_g.shape[0] == 1 else 'c'
            bfg = nm.einsum(f'cqij,{es_arg2}qjk->cqik', mtxRMI, ebf_g,
                            optimize=True)
    else:
        bfg = None

    volume = nm.sum(det, axis=1).reshape(n_el, 1, 1, 1)

    return det, volume, bfg, normal


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
        Get basis functions or their gradient evaluated in given
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

    def get_mapping(self, qp_coors, weights, bf=None, poly_space=None,
                    ori=None, transform=None, is_face=False,
                    extra=(None, None, None)):
        """
        Get the mapping for given quadrature points, weights, and
        polynomial space.

        Parameters
        ----------
        qp_coors: numpy.ndarray
            The coordinates of the integration points.
        weights:
            The integration weights.
        bf: numpy.ndarray
            The basis functions.
        poly_space: PolySpace instance
            The PolySpace instance.
        ori: numpy.ndarray
            Element orientation, used by hierarchical basis.
        transform: numpy.ndarray
            The transformation matrix applied to the basis functions.
        is_face: bool
            Is it the boundary of a region?
        extra: tuple
            The extra data for surface derivatives:
            - the derivatives of the field boundary basis functions with
              respect to the reference coordinates
            - the boundary connectivity
            - the derivatives of the domain boundary basis functions with
              respect to the reference coordinates

        Returns
        -------
        pycmap: PyCMapping instance
            The domain mapping data.
        """
        poly_space = get_default(poly_space, self.poly_space)

        bf_g = self.get_base(qp_coors, diff=True)
        if nm.allclose(bf_g, 0.0) and self.dim > 1:
            raise ValueError('zero basis function gradient!')

        if not is_face:
            ebf_g = poly_space.eval_base(qp_coors, diff=True, ori=ori,
                                         force_axis=True, transform=transform)
            size = ebf_g.nbytes * self.n_el
            site_config = Config()
            raise_if_too_large(size, site_config.refmap_memory_factor())
            se_conn, se_bf_bg = None, None
        else:
            se_conn, se_bf_bg, ebf_g = extra

        margs = eval_mapping_data_in_qp(self.coors, self.conn, self.dim,
                                        bf_g, weights, ebf_g, is_face=is_face,
                                        se_conn=se_conn, se_bf_bg=se_bf_bg)

        if bf is None:
            bf = nm.array([[[[0.]]]])
        elif len(bf.shape) == 3:
            bf = bf[None, ...]

        margs = (nm.ascontiguousarray(bf),) + margs + (self.dim,)
        pycmap = PyCMapping(*margs)

        return pycmap
