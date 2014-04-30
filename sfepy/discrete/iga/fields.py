"""
Fields for isogeometric analysis.
"""
import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete.iga.mappings import IGMapping

class IGField(Field):
    """
    Bezier extraction based NURBS field for isogeometric analysis.

    Notes
    -----
    The field has to cover the whole IGA domain.
    """
    family_name = 'volume_H1_iga'

    def __init__(self, name, dtype, shape, region, **kwargs):
        """
        Create a Bezier element isogeometric analysis field.

        Parameters
        ----------
        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or (2,)
            or 3 or (3,)) or 'vector', or a tuple. The field shape determines
            the shape of the FE base functions and is related to the number of
            components of variables and to the DOF per node count, depending
            on the field kind.
        region : Region
            The region where the field is defined.
        **kwargs : dict
            Additional keyword arguments.
        """
        shape = parse_shape(shape, region.domain.shape.dim)
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)

        self.domain = self.region.domain
        self.nurbs = self.domain.nurbs

        self._setup_kind()

        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape
        self.n_nod = self.nurbs.weights.shape[0]
        self.n_efun = nm.prod(self.nurbs.degrees + 1)

        self.mappings = {}

        self.igs = self.region.igs

    def get_true_order(self):
        return nm.prod(self.nurbs.degrees)

    def is_higher_order(self):
        """
        Return True, if the field's approximation order is greater than one.
        """
        return (self.nurbs.degrees > 1).any()

    def get_econn(self, conn_type, region, ig, is_trace=False,
                  integration=None):
        """
        Get DOF connectivity of the given type in the given region.
        """
        ct = conn_type.type if isinstance(conn_type, Struct) else conn_type

        if (ig not in self.igs) or (ig not in region.igs):
            return None

        if ct == 'volume':
            conn = self.nurbs.conn

        else:
            raise ValueError('unsupported connectivity type! (%s)' % ct)

        return conn

    def get_data_shape(self, ig, integral,
                       integration='volume', region_name=None):
        """
        Get element data dimensions.

        Parameters
        ----------
        ig : int
            The element group index.
        integral : Integral instance
            The integral describing used numerical quadrature.
        integration : 'volume', 'plate', 'surface', 'surface_extra' or 'point'
            The term integration type.
        region_name : str
            The name of surface region, required when `shape_kind` is
            'surface'.

        Returns
        -------
        data_shape : 4 ints
            The `(n_el, n_qp, dim, n_en)` for volume shape kind.

        Notes
        -----
        - `n_el` = number of elements
        - `n_qp` = number of quadrature points per element/facet
        - `dim` = spatial dimension
        - `n_en` = number of element nodes
        """
        region = self.domain.regions[region_name]
        shape = region.shape[ig]
        dim = region.dim

        _, weights = integral.get_qp(self.domain.gel.name)
        n_qp = weights.shape[0]

        data_shape = (shape.n_cell, n_qp, dim, self.n_efun)

        return data_shape

    def get_dofs_in_region_group(self, region, ig, merge=True):
        """
        Return indices of DOFs that belong to the given region and group.
        """
        dofs = []

        for idim in xrange(region.tdim, 0, -1):
            if region.can[idim]:
                break

        else:
            raise ValueError('region "%s" has no facets or cells!'
                             % region.name)

        conn = self.nurbs.conn

        if idim == region.tdim: # Cells.
            cells = region.entities[idim]
            dofs.append(nm.unique(conn[cells]))

        else: # Facets.
            assert_(idim > 0)
            fis = region.get_facet_indices(ig, offset=False, force_ig=False)

            facets = self.domain.facets[2 - idim]

            aux = []
            for ii, fi in enumerate(fis):
                aux.append(conn[fi[0], facets[fi[1]]])

            dofs.append(nm.unique(nm.concatenate(aux)))

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Set the values of given DOFs using a function of space coordinates or
        value `fun`.

        Notes
        -----
        Works for a constant value over an entire patch side only.
        """
        if region is None:
            region = self.region

        if dpn is None:
            dpn = self.n_components

        nods = []
        vals = []

        aux = self.get_dofs_in_region(region, clean=True, warn=warn)
        nods = nm.unique(nm.hstack(aux))

        if nm.isscalar(fun):
            vals = nm.repeat([fun], nods.shape[0] * dpn)

        elif isinstance(fun, nm.ndarray):
            assert_(len(fun) == dpn)
            vals = nm.repeat(fun, nods.shape[0])

        else:
            raise ValueError('unknown function/value type! (%s)' % type(fun))

        return nods, vals

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if dct != 'volume':
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def create_mapping(self, ig, region, integral, integration):
        """
        Create a new reference mapping.
        """
        vals, weights = integral.get_qp(self.domain.gel.name)

        mapping = IGMapping(self.domain, self.region.cells)
        cmap = mapping.get_mapping(vals, weights)

        return cmap, mapping
