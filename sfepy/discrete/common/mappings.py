"""
Reference-physical domain mappings.
"""
import numpy as nm

from sfepy.base.base import Struct
from sfepy.discrete.common.extmods.cmapping import CMapping


class PyCMapping(Struct):
    """
    Class for storing mapping data. Primary data in numpy arrays.
    Data for C functions translated to FMFields and embedded in CMapping.
    """
    def __init__(self, bf, det, volume, bfg, normal, dim):
        self.bf = bf
        self.det = det
        self.volume = volume
        self.bfg = bfg
        self.normal = normal

        self.cmap = CMapping(bf, det, volume, bfg, normal, dim)

        n_el, n_qp = det.shape[:2]
        n_ep = bf.shape[3]
        self.n_el = n_el
        self.n_qp = n_qp
        self.dim = dim
        self.n_ep = n_ep

    def integrate(self, out, field, mode=0):
        dim = field.shape[2]
        if mode < 3 or dim == 1:
            out[:] = nm.sum(field * self.det, axis=1)[:, None, :, :]
            if mode == 1:
                out /= self.volume
        elif dim == (self.tdim + 1) and self.normal is not None:
            out[:] = nm.dot(field, self.normal) * self.det / self.volume

        return 0


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
            if rshape[0] == 1:
                shape = (0, n_qp) + rshape[1:] # Constant parameter.

            elif (rshape[0] // n_qp) * n_qp != rshape[0]:
                raise ValueError('incompatible shapes! (n_qp: %d, %s)'
                                 % (n_qp, rshape))

            else:
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
        mapping : FEMapping or IGMapping instance
            The requested mapping.
        """
        from sfepy.discrete.fem.domain import FEDomain
        from sfepy.discrete.iga.domain import IGDomain

        if isinstance(region.domain, FEDomain):
            from sfepy.discrete.fem.mappings import FEMapping
            coors = region.domain.get_mesh_coors()
            if kind == 's':
                coors = coors[region.vertices]

            if kind == 'v':
                cells = region.get_cells()

                conn, gel = region.domain.get_conn(ret_gel=True,
                                                   tdim=region.tdim,
                                                   cells=region.cells)

                mapping = FEMapping(coors, conn, gel=gel)

            elif kind == 's':
                from sfepy.discrete.fem.fe_surface import FESurface
                aux, gel = FESurface.from_region('aux', region, ret_gel=True)
                mapping = FEMapping(coors, aux.leconn, gel=gel)

        elif isinstance(region.domain, IGDomain):
            from sfepy.discrete.iga.mappings import IGMapping
            mapping = IGMapping(region.domain, region.cells)

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
    if isinstance(gel, dict):
        qp_coors = {}
        for k, v in gel.items():
            qp, _ = integral.get_qp(v.name)
            qp_coors[k] = qp
    else:
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
