"""
Reference-physical domain mappings.
"""
import numpy as nm

from sfepy.base.base import Struct

class PhysicalQPs(Struct):
    """
    Physical quadrature points in a region.
    """

    def __init__(self, num=0):
        Struct.__init__(self, num=num, shape=(0, 0, 0))
        self.values = nm.empty(self.shape, dtype=nm.float64)

    def get_shape(self, rshape):
        """
        Get shape from raveled shape.
        """
        n_qp = self.shape[1]

        if n_qp > 0:
            if (rshape[0] / n_qp) * n_qp != rshape[0]:
                raise ValueError('incompatible shapes! (n_qp: %d, %s)'
                                 % (n_qp, rshape))

            shape = (rshape[0] // n_qp, n_qp) + rshape[1:]

        else:
            shape = (rshape[0], 0, 0, 0)

        return shape

class Mapping(Struct):
    """
    Base class for mappings.
    """

    @staticmethod
    def from_args(region, kind='v'):
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

        Returns
        -------
        mapping : VolumeMapping or SurfaceMapping instance
            The requested mapping.
        """
        from sfepy.discrete.fem.domain import FEDomain
        from sfepy.discrete.iga.domain import IGDomain

        if isinstance(region.domain, FEDomain):
            import sfepy.discrete.fem.mappings as mm
            coors = region.domain.get_mesh_coors()
            if kind == 's':
                coors = coors[region.vertices]

            conn, gel = region.domain.get_conn(ret_gel=True)

            if kind == 'v':
                cells = region.get_cells()

                mapping = mm.VolumeMapping(coors, conn[cells], gel=gel)

            elif kind == 's':
                from sfepy.discrete.fem.fe_surface import FESurface

                aux = FESurface('aux', region, gel.get_surface_entities(),
                                conn)
                mapping = mm.SurfaceMapping(coors, aux.leconn,
                                            gel=gel.surface_facet)

        elif isinstance(region.domain, IGDomain):
            import sfepy.discrete.iga.mappings as mm
            mapping = mm.IGMapping(region.domain, region.cells)

        else:
            raise ValueError('unknown domain class! (%s)' % type(region.domain))

        return mapping

def get_physical_qps(region, integral, map_kind=None):
    """
    Get physical quadrature points corresponding to the given region
    and integral.
    """
    phys_qps = PhysicalQPs()

    if map_kind is None:
        map_kind = 'v' if region.can_cells else 's'

    gmap = Mapping.from_args(region, map_kind)

    gel = gmap.get_geometry()
    qp_coors, _ = integral.get_qp(gel.name)

    qps = gmap.get_physical_qps(qp_coors)
    n_el, n_qp = qps.shape[0], qps.shape[1]

    phys_qps.num = n_el * n_qp
    phys_qps.shape = qps.shape

    qps.shape = (phys_qps.num, qps.shape[2])
    phys_qps.values = qps

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

    geo, _ = field.get_mapping(region, integral, integration)
    data = getattr(geo, name)

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
    get_mapping_data

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
    get_mapping_data

    Notes
    -----
    Assumes the same element geometry in all element groups of the field!
    """
    normals = get_mapping_data('normal', field, integral, region=region,
                               integration='surface')
    return normals
