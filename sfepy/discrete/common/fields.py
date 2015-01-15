import numpy as nm

from sfepy.base.base import output, iter_dict_of_lists, Struct, basestr

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

    elif isinstance(shape, int):
        shape = (shape,)

    return shape

def setup_extra_data(conn_info):
    """
    Setup extra data required for non-volume integration.
    """
    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        for var in info.all_vars:
            field = var.get_field()
            field.setup_extra_data(info.ps_tg, info, info.is_trace)

def fields_from_conf(conf, regions):
    fields = {}
    for key, val in conf.iteritems():
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
            import sfepy
            from sfepy.base.base import load_classes

            field_files = [ii for ii
                           in sfepy.get_paths('sfepy/discrete/fem/fields*.py')
                           if 'fields_base.py' not in ii] \
                           + sfepy.get_paths('sfepy/discrete/iga/fields*.py')
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

    def get_dofs_in_region(self, region, merge=False, clean=False,
                           warn=False, igs=None):
        """
        Return indices of DOFs that belong to the given region.
        """
        if igs is None:
            igs = region.igs

        nods = []
        for ig in self.igs:
            if not ig in igs:
                nods.append(None)
                continue

            nn = self.get_dofs_in_region_group(region, ig)
            nods.append(nn)

        if merge:
            nods = [nn for nn in nods if nn is not None]
            nods = nm.unique(nm.hstack(nods))

        elif clean:
            for nn in nods[:]:
                if nn is None:
                    nods.remove(nn)
                    if warn is not None:
                        output(warn + ('%s' % region.name))

        return nods

    def clear_mappings(self, clear_all=False):
        """
        Clear current reference mappings.
        """
        self.mappings = {}
        if clear_all:
            self.mappings0 = {}

    def save_mappings(self):
        """
        Save current reference mappings to `mappings0` attribute.
        """
        self.mappings0 = self.mappings.copy()

    def get_mapping(self, ig, region, integral, integration,
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
        key = (region.name, integral.order, ig, integration)

        if get_saved:
            out = self.mappings0.get(key, None)

        else:
            out = self.mappings.get(key, None)

        if out is None:
            out = self.create_mapping(ig, region, integral, integration)
            self.mappings[key] = out

        if return_key:
            out = out + (key,)

        return out
