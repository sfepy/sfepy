from __future__ import absolute_import

import numpy as nm

from sfepy.base.base import output, iter_dict_of_lists, Struct, basestr
from sfepy.base.timing import Timer
import six


def parse_approx_order(approx_order):
    """
    Parse the uniform approximation order value (str or int).
    """
    ao_msg = 'unsupported approximation order! (%s)'
    force_bubble = False
    discontinuous = False

    if approx_order is None:
        return 'iga', force_bubble, discontinuous

    elif isinstance(approx_order, basestr):
        if approx_order.startswith('iga'):
            return approx_order, force_bubble, discontinuous

    try:
        ao = int(approx_order)
    except ValueError:
        mode = approx_order[-1].lower()

        if mode == 'b':
            ao = int(approx_order[:-1])
            force_bubble = True

        elif mode == 'd':
            ao = int(approx_order[:-1])
            discontinuous = True

        else:
            raise ValueError(ao_msg % approx_order)

    if ao < 0:
        raise ValueError(ao_msg % approx_order)

    elif ao == 0:
        discontinuous = True

    return ao, force_bubble, discontinuous

def parse_shape(shape, dim):
    if isinstance(shape, basestr):
        try:
            shape = {'scalar' : (1,),
                     'vector' : (dim,)}[shape]
        except KeyError:
            raise ValueError('unsupported field shape! (%s)', shape)

    elif isinstance(shape, six.integer_types):
        shape = (int(shape),)

    return shape

def setup_extra_data(conn_info):
    """
    Setup extra data required for non-volume integration.
    """
    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        for var in info.all_vars:
            field = var.get_field()
            if var == info.primary:
                field.setup_extra_data(info.ps_tg, info, info.is_trace)

def fields_from_conf(conf, regions):
    fields = {}
    for key, val in six.iteritems(conf):
        field = Field.from_conf(val, regions)
        fields[field.name] = field

    return fields

class Field(Struct):
    """
    Base class for fields.
    """
    _all = None

    @staticmethod
    def from_args(name, dtype, shape, region, approx_order=1,
                  space='H1', poly_space_base='lagrange'):
        """
        Create a Field subclass instance corresponding to a given space.

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
        approx_order : int/str
            The FE approximation order, e.g. 0, 1, 2, '1B' (1 with bubble).
        space : str
            The function space name.
        poly_space_base : str
            The name of polynomial space base.

        Notes
        -----
        Assumes one cell type for the whole region!
        """
        conf = Struct(name=name, dtype=dtype, shape=shape, region=region.name,
                      approx_order=approx_order, space=space,
                      poly_space_base=poly_space_base)
        return Field.from_conf(conf, {region.name : region})

    @staticmethod
    def from_conf(conf, regions):
        """
        Create a Field subclass instance based on the configuration.
        """
        if Field._all is None:
            from sfepy import get_paths
            from sfepy.base.base import load_classes

            field_files = [ii for ii
                           in get_paths('sfepy/discrete/fem/fields*.py')
                           if 'fields_base.py' not in ii]
            field_files += get_paths('sfepy/discrete/iga/fields*.py')
            field_files += get_paths('sfepy/discrete/structural/fields*.py')
            field_files += get_paths('sfepy/discrete/dg/fields.py')
            Field._all = load_classes(field_files, [Field], ignore_errors=True,
                                      name_attr='family_name')
        table = Field._all

        space = conf.get('space', 'H1')
        poly_space_base = conf.get('poly_space_base', 'lagrange')

        key = space + '_' + poly_space_base

        approx_order = parse_approx_order(conf.approx_order)
        ao, force_bubble, discontinuous = approx_order
        region = regions[conf.region]

        if region.kind == 'cell':
            # Volume fields.
            kind = 'volume'

            if discontinuous:
                cls = table[kind + '_' + key + '_discontinuous']

            else:
                cls = table[kind + '_' + key]

            obj = cls(conf.name, conf.dtype, conf.shape, region,
                      approx_order=approx_order[:2])

        else:
            # Surface fields.
            kind = 'surface'

            cls = table[kind + '_' + key]
            obj = cls(conf.name, conf.dtype, conf.shape, region,
                      approx_order=approx_order[:2])

        return obj

    def _setup_kind(self):
        name = self.get('family_name', None,
                        'An abstract Field method called!')
        aux = name.split('_')
        self.space = aux[1]
        self.poly_space_base = aux[2]

    def clear_mappings(self, clear_all=False):
        """
        Clear current reference mappings.
        """
        self.mappings = {}
        if clear_all:
            if hasattr(self, 'mappings0'):
                self.mappings0.clear()
            else:
                self.mappings0 = {}

    def save_mappings(self):
        """
        Save current reference mappings to `mappings0` attribute.
        """
        import sfepy.base.multiproc as multi

        if multi.is_remote_dict(self.mappings0):
            for k, v in six.iteritems(self.mappings):
                m, _ = self.mappings[k]
                nv = (m.bf, m.bfg, m.det, m.volume, m.normal)
                self.mappings0[k] = nv
        else:
            self.mappings0 = self.mappings.copy()

    def get_mapping(self, region, integral, integration,
                    get_saved=False, return_key=False):
        """
        For given region, integral and integration type, get a reference
        mapping, i.e. jacobians, element volumes and base function
        derivatives for Volume-type geometries, and jacobians, normals
        and base function derivatives for Surface-type geometries
        corresponding to the field approximation.

        The mappings are cached in the field instance in `mappings`
        attribute. The mappings can be saved to `mappings0` using
        `Field.save_mappings`. The saved mapping can be retrieved by
        passing `get_saved=True`. If the required (saved) mapping
        is not in cache, a new one is created.

        Returns
        -------
        geo : CMapping instance
            The reference mapping.
        mapping : VolumeMapping or SurfaceMapping instance
            The mapping.
        key : tuple
            The key of the mapping in `mappings` or `mappings0`.
        """
        import sfepy.base.multiproc as multi

        key = (region.name, integral.order, integration)

        if get_saved:
            out = self.mappings0.get(key, None)
            if multi.is_remote_dict(self.mappings0) and out is not None:
                m, i = self.create_mapping(region, integral, integration)
                m.bf[:], m.bfg[:], m.det[:], m.volume[:] = out[0:4]
                if m.normal is not None:
                    m.normal[:] = m[4]
                out = m, i
        else:
            out = self.mappings.get(key, None)

        if out is None:
            out = self.create_mapping(region, integral, integration)
            self.mappings[key] = out

        if return_key:
            out = out + (key,)

        return out

    def create_eval_mesh(self):
        """
        Create a mesh for evaluating the field. The default implementation
        returns None, because this mesh is for most fields the same as the one
        created by `Field.create_mesh()`.
        """

    def evaluate_at(self, coors, source_vals, mode='val', strategy='general',
                    close_limit=0.1, get_cells_fun=None, cache=None,
                    ret_cells=False, ret_status=False, ret_ref_coors=False,
                    verbose=False):
        """
        Evaluate source DOF values corresponding to the field in the given
        coordinates using the field interpolation.

        Parameters
        ----------
        coors : array, shape ``(n_coor, dim)``
            The coordinates the source values should be interpolated into.
        source_vals : array, shape ``(n_nod, n_components)``
            The source DOF values corresponding to the field.
        mode : {'val', 'grad'}, optional
            The evaluation mode: the field value (default) or the field value
            gradient.
        strategy : {'general', 'convex'}, optional
            The strategy for finding the elements that contain the
            coordinates. For convex meshes, the 'convex' strategy might be
            faster than the 'general' one.
        close_limit : float, optional
            The maximum limit distance of a point from the closest
            element allowed for extrapolation.
        get_cells_fun : callable, optional
            If given, a function with signature ``get_cells_fun(coors, cmesh,
            **kwargs)`` returning cells and offsets that potentially contain
            points with the coordinates `coors`. Applicable only when
            `strategy` is 'general'. When not given,
            :func:`get_potential_cells()
            <sfepy.discrete.common.global_interp.get_potential_cells>` is used.
        cache : Struct, optional
            To speed up a sequence of evaluations, the field mesh and other
            data can be cached. Optionally, the cache can also contain the
            reference element coordinates as `cache.ref_coors`, `cache.cells`
            and `cache.status`, if the evaluation occurs in the same
            coordinates repeatedly. In that case the mesh related data are
            ignored. See :func:`Field.get_evaluate_cache()
            <sfepy.discrete.fem.fields_base.FEField.get_evaluate_cache()>`.
        ret_ref_coors : bool, optional
            If True, return also the found reference element coordinates.
        ret_status : bool, optional
            If True, return also the enclosing cell status for each point.
        ret_cells : bool, optional
            If True, return also the cell indices the coordinates are in.
        verbose : bool
            If False, reduce verbosity.

        Returns
        -------
        vals : array
            The interpolated values with shape ``(n_coor, n_components)`` or
            gradients with shape ``(n_coor, n_components, dim)`` according to
            the `mode`. If `ret_status` is False, the values where the status
            is greater than one are set to ``numpy.nan``.
        ref_coors : array
            The found reference element coordinates, if `ret_ref_coors` is True.
        cells : array
            The cell indices, if `ret_ref_coors` or `ret_cells` or `ret_status`
            are True.
        status : array
            The status, if `ret_ref_coors` or `ret_status` are True, with the
            following meaning: 0 is success, 1 is extrapolation within
            `close_limit`, 2 is extrapolation outside `close_limit`, 3 is
            failure, 4 is failure due to non-convergence of the Newton
            iteration in tensor product cells. If close_limit is 0, then for
            the 'general' strategy the status 5 indicates points outside of the
            field domain that had no potential cells.
        """
        from sfepy.discrete.common.global_interp import get_ref_coors
        from sfepy.discrete.common.extmods.crefcoors import evaluate_in_rc
        from sfepy.base.base import complex_types

        output('evaluating in %d points...' % coors.shape[0], verbose=verbose)

        ref_coors, cells, status = get_ref_coors(self, coors,
                                                 strategy=strategy,
                                                 close_limit=close_limit,
                                                 get_cells_fun=get_cells_fun,
                                                 cache=cache,
                                                 verbose=verbose)

        timer = Timer(start=True)

        # Interpolate to the reference coordinates.
        source_dtype = nm.float64 if source_vals.dtype in complex_types\
            else source_vals.dtype
        if mode == 'val':
            vals = nm.empty((coors.shape[0], source_vals.shape[1], 1),
                            dtype=source_dtype)
            cmode = 0

        elif mode == 'grad':
            vals = nm.empty((coors.shape[0], source_vals.shape[1],
                             coors.shape[1]),
                            dtype=source_dtype)
            cmode = 1

        ctx = self.create_basis_context()

        if source_vals.dtype in complex_types:
            valsi = vals.copy()
            evaluate_in_rc(vals, ref_coors, cells, status,
                           nm.ascontiguousarray(source_vals.real),
                           self.get_econn('volume', self.region), cmode, ctx)
            evaluate_in_rc(valsi, ref_coors, cells, status,
                           nm.ascontiguousarray(source_vals.imag),
                           self.get_econn('volume', self.region), cmode, ctx)
            vals = vals + valsi * 1j
        else:
            evaluate_in_rc(vals, ref_coors, cells, status, source_vals,
                           self.get_econn('volume', self.region), cmode, ctx)

        output('interpolation: %f s' % timer.stop(),verbose=verbose)

        output('...done',verbose=verbose)

        if mode == 'val':
            vals.shape = (coors.shape[0], source_vals.shape[1])

        if not ret_status:
            ii = nm.where(status > 1)[0]
            vals[ii] = nm.nan

        if ret_ref_coors:
            return vals, ref_coors, cells, status

        elif ret_status:
            return vals, cells, status

        elif ret_cells:
            return vals, cells

        else:
            return vals
