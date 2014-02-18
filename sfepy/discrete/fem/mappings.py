"""
Finite element reference mappings.
"""
import numpy as nm

from sfepy.base.base import output, get_default, Struct
from sfepy.discrete.fem.poly_spaces import PolySpace
from sfepy.discrete.fem.fe_surface import FESurface
from sfepy.discrete.fem.extmods.mappings import CMapping

class PhysicalQPs(Struct):
    """
    Physical quadrature points in a region.
    """
    def __init__(self, igs, n_total=0, is_uniform=True):
        Struct.__init__(self, igs=igs, n_total=n_total, indx={}, rindx={},
                        n_per_group={}, shape={}, values={},
                        is_uniform=is_uniform)
        for ig in self.igs:
            self.indx[ig] = slice(None)
            self.rindx[ig] = slice(None)
            self.n_per_group[ig] = 0
            self.shape[ig] = (0, 0, 0)
            self.values[ig] = nm.empty(self.shape[ig], dtype=nm.float64)

    def get_merged_values(self):
        qps = nm.concatenate([self.values[ig] for ig in self.igs], axis=0)

        return qps

    def get_shape(self, rshape, ig=None):
        """
        Get shape from raveled shape.
        """
        if ig is None:
            if self.is_uniform:
                n_qp = self.shape[self.igs[0]][1]

            else:
                msg = 'ig argument must be given for non-uniform QPs!'
                raise ValueError(msg)

        else:
            n_qp = self.shape[ig][1]

        if (rshape[0] / n_qp) * n_qp != rshape[0]:
            raise ValueError('incompatible shapes! (n_qp: %d, %s)'
                             % (n_qp, rshape))

        shape = (rshape[0] / n_qp, n_qp) + rshape[1:]

        return shape

def create_mapping(region, kind, ig):
    """
    Create mapping from reference to physical entities in a given
    region, given the integration kind ('v' or 's').

    This mapping can be used to compute the physical quadrature
    points.

    Parameters
    ----------
    region : Region instance
        The region defining the entities.
    kind : 'v' or 's'
        The kind of the entities: 'v' - cells, 's' - facets.
    ig : int
        The group index.

    Returns
    -------
    mapping : VolumeMapping or SurfaceMapping instance
        The requested mapping.
    """
    coors = region.domain.get_mesh_coors()
    if kind == 's':
        coors = coors[region.vertices]

    gel = region.domain.groups[ig].gel
    conn = region.domain.groups[ig].conn

    if kind == 'v':
        cells = region.get_cells(ig)

        mapping = VolumeMapping(coors, conn[cells], gel=gel)

    elif kind == 's':
        aux = FESurface('aux', region, gel.get_surface_entities(),
                        conn , ig)
        mapping = SurfaceMapping(coors, aux.leconn, gel=gel.surface_facet)

    return mapping

def get_physical_qps(region, integral, map_kind=None):
    """
    Get physical quadrature points corresponding to the given region
    and integral.
    """
    phys_qps = PhysicalQPs(region.igs)

    if map_kind is None:
        map_kind = 'v' if region.can_cells else 's'

    ii = 0
    for ig in region.igs:
        gmap = create_mapping(region, map_kind, ig)

        gel = gmap.get_geometry()
        qp_coors, _ = integral.get_qp(gel.name)

        qps = gmap.get_physical_qps(qp_coors)
        n_el, n_qp = qps.shape[0], qps.shape[1]

        phys_qps.n_per_group[ig] = n_per_group = n_el * n_qp
        phys_qps.shape[ig] = qps.shape

        phys_qps.indx[ig] = slice(ii, ii + n_el)
        phys_qps.rindx[ig] = slice(ii * n_qp, (ii + n_el) * n_qp)

        ii += qps.shape[0]

        qps.shape = (n_per_group, qps.shape[2])
        phys_qps.values[ig] = qps
        phys_qps.n_total += n_el * n_qp

    return phys_qps

def get_mapping_data(name, field, integral, region=None, integration='volume'):
    """
    General helper function for accessing reference mapping data.

    Get data attribute `name` from reference mapping corresponding to
    `field` in `region` in quadrature points of the given `integral` and
    `integration` type.

    Parameters
    ----------
    name : str
        The reference mapping attribute name.
    field : Field instance
        The field defining the reference mapping.
    integral : Integral instance
        The integral defining quadrature points.
    region : Region instance, optional
        If given, use the given region instead of `field` region.
    integration : one of ('volume', 'surface', 'surface_extra')
        The integration type.

    Returns
    -------
    data : array
        The required data merged for all element groups.

    Notes
    -----
    Assumes the same element geometry in all element groups of the field!
    """
    data = None
    if region is None:
        region = field.region
    for ig in region.igs:
        geo, _ = field.get_mapping(ig, region, integral, integration)
        _data = getattr(geo, name)
        if data is None:
            data = _data

        else:
            data = nm.concatenate((data, _data), axis=0)

    return data

def get_jacobian(field, integral, region=None, integration='volume'):
    """
    Get the jacobian of reference mapping corresponding to `field`.

    Parameters
    ----------
    field : Field instance
        The field defining the reference mapping.
    integral : Integral instance
        The integral defining quadrature points.
    region : Region instance, optional
        If given, use the given region instead of `field` region.
    integration : one of ('volume', 'surface', 'surface_extra')
        The integration type.

    Returns
    -------
    jac : array
        The jacobian merged for all element groups.

    See Also
    --------
    get_mapping_data()

    Notes
    -----
    Assumes the same element geometry in all element groups of the field!
    """
    jac = get_mapping_data('det', field, integral, region=region,
                           integration=integration)
    return jac

def get_normals(field, integral, region):
    """
    Get the normals of element faces in `region`.

    Parameters
    ----------
    field : Field instance
        The field defining the reference mapping.
    integral : Integral instance
        The integral defining quadrature points.
    region : Region instance
        The given of the element faces.

    Returns
    -------
    normals : array
        The normals merged for all element groups.

    See Also
    --------
    get_mapping_data()

    Notes
    -----
    Assumes the same element geometry in all element groups of the field!
    """
    normals = get_mapping_data('normal', field, integral, region=region,
                               integration='surface')
    return normals

class Mapping(Struct):
    """
    Base class for mappings.
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

    def get_mapping(self, qp_coors, weights, poly_space=None, ori=None):
        """
        Get the mapping for given quadrature points, weights, and
        polynomial space.

        Returns
        -------
        cmap : CMapping instance
            The volume mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        bf_g = self.get_base(qp_coors, diff=True)

        ebf_g = poly_space.eval_base(qp_coors, diff=True, ori=ori,
                                     force_axis=True)
        flag = ori is not None

        cmap = CMapping(self.n_el, qp_coors.shape[0], self.dim,
                        poly_space.n_nod, mode='volume', flag=flag)
        cmap.describe(self.coors, self.conn, bf_g, ebf_g, weights)

        return cmap

class SurfaceMapping(Mapping):
    """
    Mapping from reference domain to physical domain of the space
    dimension higher by one.
    """

    def get_mapping(self, qp_coors, weights, poly_space=None, mode='surface'):
        """
        Get the mapping for given quadrature points, weights, and
        polynomial space.

        Returns
        -------
        cmap : CMapping instance
            The surface mapping.
        """
        poly_space = get_default(poly_space, self.poly_space)

        bf_g = self.get_base(qp_coors, diff=True)

        if nm.allclose(bf_g, 0.0):
            raise ValueError('zero base function gradient!')

        cmap = CMapping(self.n_el, qp_coors.shape[0], self.dim,
                        poly_space.n_nod, mode=mode)
        cmap.describe(self.coors, self.conn, bf_g, None, weights)

        return cmap
