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
        self.is_surface = False

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

        nconn = self.nurbs.conn

        if ct == 'volume':
            if region.name == self.region.name:
                conn = nconn

            else:
                cells = region.get_cells(ig, true_cells_only=True)
                conn = nm.take(nconn, cells.astype(nm.int32), axis=0)

        elif ct == 'surface':
            fis = region.get_facet_indices(ig, offset=False, force_ig=False)
            tdim = region.kind_tdim
            facets = self.domain.facets[2 - tdim]

            conn = []
            for ii, fi in enumerate(fis):
                conn.append(nconn[fi[0], facets[fi[1]]])

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

        Notes
        -----
        `ig`, `merge` are not used.
        """
        idim = region.kind_tdim
        if idim < (self.domain.shape.tdim - 1):
            raise ValueError('region "%s" has no facets or cells!'
                             % region.name)

        if idim == region.tdim: # Cells.
            rconn = self.get_econn('volume', region, ig)
            dofs = nm.unique(rconn)

        else: # Facets.
            rconn = self.get_econn('surface', region, ig)
            dofs = nm.unique(nm.concatenate(rconn))

        if not merge:
            dofs = [dofs]

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
        mapping = IGMapping(self.domain, region.cells)
        cmap = mapping.get_mapping(vals, weights)

        return cmap, mapping

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, **kwargs):
        """
        Convert the DOFs corresponding to the field to a dictionary of
        output data usable by Mesh.write().

        Parameters
        ----------
        dofs : array, shape (n_nod, n_component)
            The array of DOFs reshaped so that each column corresponds
            to one component.
        var_name : str
            The variable name corresponding to `dofs`.
        dof_names : tuple of str
            The names of DOF components.
        key : str, optional
            The key to be used in the output dictionary instead of the
            variable name.

        Returns
        -------
        out : dict
            The output dictionary.
        """
        from sfepy.discrete.iga.utils import create_mesh_and_output

        key = key if key is not None else var_name

        num = 25 if self.region.dim == 3 else 101
        pars = (nm.linspace(0, 1, num),) * self.region.dim
        mesh, out = create_mesh_and_output(self.domain.nurbs, pars,
                                           **{key : dofs})
        out[key].var_name = var_name
        out[key].dofs = dof_names
        out['__mesh__'] = mesh

        return out
