import numpy as nm

from sfepy.base.base import Struct
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete import PolySpace
from sfepy.discrete.fem.mappings import FEMapping
from sfepy.discrete.fem.fields_base import _find_geometry
#from sfepy.discrete.fem.utils import get_min_value
from sfepy.discrete.fem.meshio import convert_complex_output

class L2ConstantVolumeField(Field):
    """
    The L2 constant-in-a-region approximation.
    """
    family_name = 'volume_L2_constant'

    def __init__(self, name, dtype, shape, region, approx_order=0):
        """
        Create a L2 constant field.

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
        approx_order : int or tuple
            The FE approximation order. Ignored here.
        """
        field_dim = region.field_dim if hasattr(region, 'field_dim')\
            else region.domain.shape.dim
        shape = parse_shape(shape, field_dim)
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region, approx_order=0)
        self.domain = self.region.domain
        self.cmesh = self.region.cmesh

        self.gel, self.is_surface = _find_geometry(self.region)
        self._setup_kind()
        self._create_interpolant()

        cells = self.region.get_cells(true_cells_only=False)
        self.econn = nm.zeros((len(cells), 1), dtype=nm.int32)

        self.domain = self.region.domain

        self.n_components = int(nm.prod(self.shape))
        self.val_shape = self.shape
        self.n_nod = 1

        self.extra_data = {}
        self.mappings = {}

    def _create_interpolant(self):
        name = '%s_%s_%s_%d' % (self.gel.name, self.space,
                                self.poly_space_basis, self.approx_order)
        ps = PolySpace.any_from_args(name, self.gel, self.approx_order,
                                     basis='lagrange')
        self.poly_space = ps

    def setup_extra_data(self, info):
        pass

    def get_coor(self, nods=None):
        """
        Returns the barycenter of the field region.

        Parameters
        ----------
        nods : array, optional
            Ignored.
        """
        coors = nm.sum(self.domain.mesh.coors[self.region.vertices],
                       axis=0, keepdims=True)
        coors /= self.region.shape.n_vertex
        return coors

    def get_econn(self, conn_type, region, trace_region=None, local=False):
        """
        Get extended connectivity of the given type in the given region.

        Parameters
        ----------
        conn_type: tuple or string
            DOF connectivity type, ignored.
        region: sfepy.discrete.common.region.Region
            The region for which the connectivity is required.
        trace_region: None or string
            Ignored.
        local: bool
            Ignored.

        Returns
        -------
        econn: numpy.ndarray
            The extended connectivity array.
        """
        if region.name == self.region.name:
            conn = self.econn

        else:
            cells = region.get_cells(true_cells_only=False)
            conn = nm.zeros((len(cells), 1), dtype=nm.int32)

        return conn

    def get_data_shape(self, integral, integration='cell', region_name=None):
        """
        Get element data dimensions.

        Parameters
        ----------
        integral : Integral instance
            The integral describing used numerical quadrature.
        integration : 'cell'
            The term integration type. Ignored.
        region_name : str
            The name of the region of the integral.

        Returns
        -------
        data_shape : 4 ints
            The `(n_el, n_qp, dim, n_en)` for volume shape kind,
            `(n_fa, n_qp, dim, n_fn)` for surface shape kind.

        Notes
        -----
        - `n_el`, `n_fa` = number of elements/facets
        - `n_qp` = number of quadrature points per element/facet
        - `dim` = spatial dimension
        - `n_en`, `n_fn` = number of element/facet nodes
        - `n_nod` = number of element nodes
        """
        if integration not in ('cell', 'facet', 'facet_extra'):
            raise NotImplementedError('unsupported integration type! (%s)'
                                      % integration)

        region = self.domain.regions[region_name]
        n_cell = len(region.entities[region.kind_tdim])
        _, weights = integral.get_qp(self.gel.name)
        n_qp = weights.shape[0]
        dim = region.field_dim if hasattr(region, 'field_dim') else region.dim
        return (n_cell, n_qp, dim, 1)

    def get_dofs_in_region(self, region, merge=True):
        """
        Return indices of DOFs that belong to the given region.
        """
        return nm.zeros(1, dtype=nm.int32)

    def eval_basis(self, key, derivative, integral, iels=None,
                   from_geometry=False, base_only=True):
        qp_coors, qp_weights = integral.get_qp(self.gel.name)
        ps = self.poly_space
        bf = ps.eval_basis(qp_coors)
        if key[0] == 'b': # BQP
            num = 6 if self.gel.n_vertex == 4 else 4
            bf = nm.tile(bf, (num, 1, 1, 1))

        if base_only:
            return bf

        else:
            return bf, qp_weights

    def create_mapping(self, region, integral, integration,
                       return_mapping=True):
        domain = self.domain
        coors = domain.get_mesh_coors(actual=True)
        iels = region.get_cells(true_cells_only=(region.kind == 'cell'))

        if integration == 'cell':
            ps = self.poly_space
            geo_ps = self.gel.poly_space
            dconn = domain.get_conn(tdim=region.tdim, cells=iels)
            qp_coors, qp_weights = integral.get_qp(self.gel.name)
            bf = ps.eval_basis(qp_coors)

        elif integration in ('facet', 'facet_extra'):
            if self.is_surface:
                gel = self.gel
                ps = self.poly_space

            else:
                gel = self.gel.surface_facet
                ps = PolySpace.any_from_args('aux', gel, self.approx_order,
                                             basis='lagrange')

            geo_ps = gel.poly_space
            domain.create_surface_group(region)
            sd = domain.surface_groups[region.name]
            dconn = sd.get_connectivity()
            qp_coors, qp_weights = integral.get_qp(gel.name)
            bf = ps.eval_basis(qp_coors)

        mapping = FEMapping(coors, dconn, poly_space=geo_ps)
        out = mapping.get_mapping(qp_coors, qp_weights, bf, poly_space=ps,
                                  is_face=self.is_surface)

        if return_mapping:
            out = (out, mapping)

        return out

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, extend=True, fill_value=None,
                      linearization=None):
        _new_dofs = nm.tile(dofs, (self.region.shape.n_vertex, 1))
        if extend:
            if fill_value is None:
                fill_value = 0.0

            n_nod = self.domain.shape.n_nod
            new_dofs = nm.full((n_nod, dofs.shape[1]), fill_value,
                               dtype=self.dtype)
            new_dofs[self.region.vertices] = _new_dofs

        else:
            new_dofs = _new_dofs

        out = {}
        out[key] = Struct(name='output_data', mode='vertex',
                          data=new_dofs, var_name=var_name,
                          dofs=dof_names, region_name=self.region.name)

        out = convert_complex_output(out)
        return out

class L2ConstantSurfaceField(L2ConstantVolumeField):
    family_name = 'surface_L2_constant'
