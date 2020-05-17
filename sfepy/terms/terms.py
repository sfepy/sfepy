from __future__ import absolute_import
import re
from copy import copy

import numpy as nm

from sfepy.base.base import (as_float_or_complex, get_default, assert_,
                             Container, Struct, basestr, goptions)
from sfepy.base.compat import in1d

# Used for imports in term files.
from sfepy.terms.extmods import terms
import six
from six.moves import range
from functools import reduce

_match_args = re.compile('^([^\(\}]*)\((.*)\)$').match
_match_virtual = re.compile('^virtual$').match
_match_state = re.compile('^state(_[_a-zA-Z0-9]+)?$').match
_match_parameter = re.compile('^parameter(_[_a-zA-Z0-9]+)?$').match
_match_material = re.compile('^material(_[_a-zA-Z0-9]+)?$').match
_match_material_opt = re.compile('^opt_material(_[_a-zA-Z0-9]+)?$').match
_match_material_root = re.compile('(.+)\.(.*)').match
_match_ts = re.compile('^ts$').match

def get_arg_kinds(arg_types):
    """
    Translate `arg_types` of a Term to a canonical form.

    Parameters
    ----------
    arg_types : tuple of strings
        The term argument types, as given in the `arg_types` attribute.

    Returns
    -------
    arg_kinds : list of strings
        The argument kinds - one of 'virtual_variable', 'state_variable',
        'parameter_variable', 'opt_material', 'ts', 'user'.
    """
    arg_kinds = []
    for ii, arg_type in enumerate(arg_types):
        if _match_virtual(arg_type):
            arg_kinds.append('virtual_variable')

        elif _match_state(arg_type):
            arg_kinds.append('state_variable')

        elif _match_parameter(arg_type):
            arg_kinds.append('parameter_variable')

        elif _match_material(arg_type):
            arg_kinds.append('material')

        elif _match_material_opt(arg_type):
            arg_kinds.append('opt_material')
            if ii > 0:
                msg = 'opt_material at position %d, must be at 0!' % ii
                raise ValueError(msg)

        elif _match_ts(arg_type):
            arg_kinds.append('ts')

        else:
            arg_kinds.append('user')

    return arg_kinds

def get_shape_kind(integration):
    """
    Get data shape kind for given integration type.
    """
    if integration == 'surface':
        shape_kind = 'surface'

    elif integration in ('volume', 'plate', 'surface_extra'):
        shape_kind = 'volume'

    elif integration == 'point':
        shape_kind = 'point'

    else:
        raise NotImplementedError('unsupported term integration! (%s)'
                                  % integration)

    return shape_kind

def split_complex_args(args):
    """
    Split complex arguments to real and imaginary parts.

    Returns
    -------
    newargs : dictionary
        Dictionary with lists corresponding to `args` such that each
        argument of numpy.complex128 data type is split to its real and
        imaginary part. The output depends on the number of complex
        arguments in 'args':

          - 0: list (key 'r') identical to input one

          - 1: two lists with keys 'r', 'i' corresponding to real
            and imaginary parts

          - 2: output dictionary contains four lists:

            - 'r' - real(arg1), real(arg2)
            - 'i' - imag(arg1), imag(arg2)
            - 'ri' - real(arg1), imag(arg2)
            - 'ir' - imag(arg1), real(arg2)
    """
    newargs = {}
    cai = []

    for ii, arg in enumerate(args):
        if isinstance(arg, nm.ndarray) and (arg.dtype == nm.complex128):
            cai.append(ii)

    if len(cai) > 0:
        newargs['r'] = list(args[:])
        newargs['i'] = list(args[:])

        arg1 = cai[0]
        newargs['r'][arg1] = args[arg1].real.copy()
        newargs['i'][arg1] = args[arg1].imag.copy()

        if len(cai) == 2:
            arg2 = cai[1]
            newargs['r'][arg2] = args[arg2].real.copy()
            newargs['i'][arg2] = args[arg2].imag.copy()

            newargs['ri'] = list(args[:])
            newargs['ir'] = list(args[:])
            newargs['ri'][arg1] = newargs['r'][arg1]
            newargs['ri'][arg2] = newargs['i'][arg2]
            newargs['ir'][arg1] = newargs['i'][arg1]
            newargs['ir'][arg2] = newargs['r'][arg2]

        elif len(cai) > 2:
            raise NotImplementedError('more than 2 complex arguments! (%d)'
                                      % len(cai))

    else:
        newargs['r'] = args[:]

    return newargs

def create_arg_parser():
    from pyparsing import Literal, Word, delimitedList, Group, \
         StringStart, StringEnd, Optional, nums, alphas, alphanums

    ident = Word(alphas, alphanums + "_")
    inumber = Word("+-" + nums, nums)

    history = Optional(Literal('[').suppress() + inumber
                       + Literal(']').suppress(), default=0)("history")
    history.setParseAction(lambda str, loc, toks: int(toks[0]))

    variable = Group(Word(alphas, alphanums + '._') + history)

    derivative = Group(Literal('d') + variable\
                       + Literal('/').suppress() + Literal('dt'))

    trace = Group(Literal('tr')
                  + Literal('(').suppress()
                  + Optional(ident + Literal(',').suppress(), default=None)
                  + variable
                  + Literal(')').suppress())

    generalized_var = derivative | trace | variable

    args = StringStart() + delimitedList(generalized_var) + StringEnd()

    return args

class ConnInfo(Struct):

    def get_region(self, can_trace=True):
        if self.is_trace and can_trace:
            return self.region.get_mirror_region(self.trace_region)
        else:
            return self.region

    def get_region_name(self, can_trace=True):
        if self.is_trace and can_trace:
            reg = self.region.get_mirror_region(self.trace_region)
        else:
            reg = self.region

        if reg is not None:
            return reg.name
        else:
            return None

class Terms(Container):

    @staticmethod
    def from_desc(term_descs, regions, integrals=None):
        """
        Create terms, assign each term its region.
        """
        from sfepy.terms import term_table

        terms = Terms()
        for td in term_descs:
            try:
                constructor = term_table[td.name]
            except:
                msg = "term '%s' is not in %s" % (td.name,
                                                  sorted(term_table.keys()))
                raise ValueError(msg)

            try:
                region = regions[td.region]
            except IndexError:
                raise KeyError('region "%s" does not exist!' % td.region)

            term = Term.from_desc(constructor, td, region, integrals=integrals)
            terms.append(term)

        return terms

    def __init__(self, objs=None):
        Container.__init__(self, objs=objs)

        self.update_expression()

    def insert(self, ii, obj):
        Container.insert(self, ii, obj)
        self.update_expression()

    def append(self, obj):
        Container.append(self, obj)
        self.update_expression()

    def update_expression(self):
        self.expression = []
        for term in self:
            aux = [term.sign, term.name, term.arg_str,
                   term.integral_name, term.region.name]
            self.expression.append(aux)

    def __mul__(self, other):
        out = Terms()
        for name, term in self.iteritems():
            out.append(term * other)

        return out

    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        if isinstance(other, Term):
            out = self.copy()
            out.append(other)

        elif isinstance(other, Terms):
            out = Terms(self._objs + other._objs)

        else:
            raise ValueError('cannot add Terms with %s!' % other)

        return out

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, Term):
            out = self + (-other)

        elif isinstance(other, Terms):
            out = self + (-other)

        else:
            raise ValueError('cannot subtract Terms with %s!' % other)

        return out

    def __rsub__(self, other):
        return -self + other

    def __pos__(self):
        return self

    def __neg__(self):
        return -1.0 * self

    def setup(self):
        for term in self:
            term.setup()

    def assign_args(self, variables, materials, user=None):
        """
        Assign all term arguments.
        """
        for term in self:
            term.assign_args(variables, materials, user)

    def get_variable_names(self):
        out = []
        for term in self:
            out.extend(term.get_variable_names())
        return list(set(out))

    def get_material_names(self):
        out = []
        for term in self:
            out.extend(term.get_material_names())
        return list(set(out))

    def get_user_names(self):
        out = []
        for term in self:
            out.extend(term.get_user_names())
        return list(set(out))

class Term(Struct):
    name = ''
    arg_types = ()
    arg_shapes = {}
    integration = 'volume'
    geometries = ['1_2', '2_3', '2_4', '3_4', '3_8']

    @staticmethod
    def new(name, integral, region, **kwargs):
        from sfepy.terms import term_table

        arg_str = _match_args(name)
        if arg_str is not None:
            name, arg_str = arg_str.groups()

        else:
            raise ValueError('bad term syntax! (%s)' % name)

        if name in term_table:
            constructor = term_table[name]

        else:
            msg = "term '%s' is not in %s" % (name, sorted(term_table.keys()))
            raise ValueError(msg)

        obj = constructor(name, arg_str, integral, region, **kwargs)
        return obj

    @staticmethod
    def from_desc(constructor, desc, region, integrals=None):
        from sfepy.discrete import Integrals

        if integrals is None:
            integrals = Integrals()

        integral = integrals.get(desc.integral)
        obj = constructor(desc.name, desc.args, integral, region)
        obj.sign = desc.sign

        return obj

    def __init__(self, name, arg_str, integral, region, **kwargs):
        self.name = name
        self.arg_str = arg_str
        self.region = region
        self._kwargs = kwargs
        self._integration = self.integration
        self.sign = 1.0

        self.set_integral(integral)

    def __mul__(self, other):
        try:
            mul = as_float_or_complex(other)

        except ValueError:
            raise ValueError('cannot multiply Term with %s!' % other)

        out = self.copy(name=self.name)

        out.sign = mul * self.sign

        return out

    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        if isinstance(other, Term):
            out = Terms([self, other])

        else:
            out = NotImplemented

        return out

    def __sub__(self, other):
        if isinstance(other, Term):
            out = Terms([self, -1.0 * other])

        else:
            out = NotImplemented

        return out

    def __pos__(self):
        return self

    def __neg__(self):
        out = -1.0 * self
        return out

    def get_str(self):
        return '{:+} * {}.{}.{}({})'.format(
            self.sign, self.name, self.integral.order,
            self.region.name, self.arg_str)

    def set_integral(self, integral):
        """
        Set the term integral.
        """
        self.integral = integral
        if self.integral is not None:
            self.integral_name = self.integral.name

    def setup(self):
        self.function = Struct.get(self, 'function', None)

        self.step = 0
        self.dt = 1.0
        self.is_quasistatic = False
        self.has_region = True

        self.setup_formal_args()

        if self._kwargs:
            self.setup_args(**self._kwargs)

        else:
            self.args = []

    def setup_formal_args(self):
        self.arg_names = []
        self.arg_steps = {}
        self.arg_derivatives = {}
        self.arg_traces = {}
        self.arg_trace_regions = {}

        parser = create_arg_parser()
        self.arg_desc = parser.parseString(self.arg_str)
        for arg in self.arg_desc:
            derivative = None
            trace = False
            trace_region = None

            if isinstance(arg[1], int):
                name, step = arg

            else:
                kind = arg[0]
                if kind == 'd':
                    name, step = arg[1]
                    derivative = arg[2]
                elif kind == 'tr':
                    trace = True
                    trace_region = arg[1]
                    name, step = arg[2]

            match = _match_material_root(name)
            if match:
                name = (match.group(1), match.group(2))

            self.arg_names.append(name)
            self.arg_steps[name] = step
            self.arg_derivatives[name] = derivative
            self.arg_traces[name] = trace
            self.arg_trace_regions[name] = trace_region

    def setup_args(self, **kwargs):
        self._kwargs = kwargs

        self.args = []
        for arg_name in self.arg_names:
            if isinstance(arg_name, basestr):
                self.args.append(self._kwargs[arg_name])

            else:
                self.args.append((self._kwargs[arg_name[0]], arg_name[1]))

        self.classify_args()
        self.check_args()

    def assign_args(self, variables, materials, user=None):
        """
        Check term argument existence in variables, materials, user data
        and assign the arguments to terms. Also check compatibility of
        field and term regions.
        """
        if user is None:
            user = {}

        user.setdefault('ts', Struct())

        kwargs = {}
        for arg_name in self.arg_names:
            if isinstance(arg_name, basestr):
                if arg_name in variables.names:
                    kwargs[arg_name] = variables[arg_name]

                elif arg_name in user:
                    kwargs[arg_name] = user[arg_name]

                else:
                    raise ValueError('argument %s not found!' % arg_name)

            else:
                arg_name = arg_name[0]
                if arg_name in materials.names:
                    kwargs[arg_name] = materials[arg_name]

                else:
                    raise ValueError('material argument %s not found!'
                                     % arg_name)

        self.setup_args(**kwargs)

    def classify_args(self):
        """
        Classify types of the term arguments and find matching call
        signature.

        A state variable can be in place of a parameter variable and
        vice versa.
        """
        self.names = Struct(name='arg_names',
                            material=[], variable=[], user=[],
                            state=[], virtual=[], parameter=[])

        # Prepare for 'opt_material' - just prepend a None argument if needed.
        if isinstance(self.arg_types[0], tuple):
            arg_types = self.arg_types[0]

        else:
            arg_types = self.arg_types

        if len(arg_types) == (len(self.args) + 1):
            self.args.insert(0, (None, None))
            self.arg_names.insert(0, (None, None))

        if isinstance(self.arg_types[0], tuple):
            assert_(len(self.modes) == len(self.arg_types))
            # Find matching call signature using variable arguments - material
            # and user arguments are ignored!
            matched = []
            for it, arg_types in enumerate(self.arg_types):
                arg_kinds = get_arg_kinds(arg_types)
                if self._check_variables(arg_kinds):
                    matched.append((it, arg_kinds))

            if len(matched) == 1:
                i_match, arg_kinds = matched[0]
                arg_types = self.arg_types[i_match]
                self.mode = self.modes[i_match]

            elif len(matched) == 0:
                msg = 'cannot match arguments! (%s)' % self.arg_names
                raise ValueError(msg)

            else:
                msg = 'ambiguous arguments! (%s)' % self.arg_names
                raise ValueError(msg)

        else:
            arg_types = self.arg_types
            arg_kinds = get_arg_kinds(self.arg_types)
            self.mode = Struct.get(self, 'mode', None)

            if not self._check_variables(arg_kinds):
                raise ValueError('cannot match variables! (%s)'
                                 % self.arg_names)

        # Set actual argument types.
        self.ats = list(arg_types)

        for ii, arg_kind in enumerate(arg_kinds):
            name = self.arg_names[ii]
            if arg_kind.endswith('variable'):
                names = self.names.variable

                if arg_kind == 'virtual_variable':
                    self.names.virtual.append(name)

                elif arg_kind == 'state_variable':
                    self.names.state.append(name)

                elif arg_kind == 'parameter_variable':
                    self.names.parameter.append(name)

            elif arg_kind.endswith('material'):
                # This should be better checked already in create_arg_parser().
                if not isinstance(name, tuple):
                    raise ValueError('wrong material argument %s of term %s!'
                                     % (name, self.get_str()))
                names = self.names.material

            else:
                names = self.names.user

            names.append(name)

        self.n_virtual = len(self.names.virtual)
        if self.n_virtual > 1:
            raise ValueError('at most one virtual variable is allowed! (%d)'
                             % self.n_virtual)

        self.set_arg_types()

        self.setup_integration()

    def _check_variables(self, arg_kinds):
        for ii, arg_kind in enumerate(arg_kinds):
            if arg_kind.endswith('variable'):
                var = self.args[ii]
                check = {'virtual_variable' : var.is_virtual,
                         'state_variable' : var.is_state_or_parameter,
                         'parameter_variable' : var.is_state_or_parameter}
                if not check[arg_kind]():
                    return False

        else:
            return True

    def set_arg_types(self):
        pass

    def check_args(self):
        """
        Common checking to all terms.

        Check compatibility of field and term regions.
        """
        vns = self.get_variable_names()
        for name in vns:
            field = self._kwargs[name].get_field()
            if field is None:
                continue

            region = self.region
            if self.arg_traces[name]:
                mreg_name = self.arg_trace_regions[name]
                if mreg_name is None:
                    mreg_name = region.setup_mirror_region(mreg_name,
                                                           ret_name=True)
                    self.arg_trace_regions[name] = mreg_name
                else:
                    region.setup_mirror_region(mreg_name)

                region = region.get_mirror_region(mreg_name)

            if not nm.all(in1d(region.vertices,
                               field.region.vertices)):
                msg = ('%s: incompatible regions: (self, field %s)'
                       + '(%s in %s)') %\
                       (self.name, field.name,
                        self.region.vertices, field.region.vertices)
                raise ValueError(msg)

    def get_variable_names(self):
        return self.names.variable

    def get_material_names(self):
        out = []
        for aux in self.names.material:
            if aux[0] is not None:
                out.append(aux[0])
        return out

    def get_user_names(self):
        return self.names.user

    def get_virtual_name(self):
        if not self.names.virtual:
            return None

        var = self.get_virtual_variable()
        return var.name

    def get_state_names(self):
        """
        If variables are given, return only true unknowns whose data are of
        the current time step (0).
        """
        variables = self.get_state_variables()
        return [var.name for var in variables]

    def get_parameter_names(self):
        return copy(self.names.parameter)

    def get_conn_key(self):
        """The key to be used in DOF connectivity information."""
        key = (self.name,) + tuple(self.arg_names)
        arg_traces = [k for k, v in self.arg_traces.items() if v]
        if len(arg_traces) > 0:
            atr = arg_traces[-1]
            trace = True, self.arg_trace_regions[atr], atr
        else:
            trace = False, None, None

        key += (self.integral_name, self.region.name) + trace

        return key

    def get_conn_info(self):
        vvar = self.get_virtual_variable()
        svars = self.get_state_variables()
        pvars = self.get_parameter_variables()

        all_vars = self.get_variables()

        dc_type = self.get_dof_conn_type()
        tgs = self.get_geometry_types()

        v_tg = None
        if vvar is not None:
            field = vvar.get_field()
            if field is not None:
                if vvar.name in tgs:
                    v_tg = tgs[vvar.name]

                else:
                    v_tg = None

        else:
            # No virtual variable -> all unknowns are in fact known parameters.
            pvars += svars
            svars = []

        region = self.get_region()
        if region is not None:
            arg_traces = [k for k, v in self.arg_traces.items() if v]
            if len(arg_traces) > 0:
                aname = arg_traces[-1]
                mreg_name = self.arg_trace_regions[aname]
                if mreg_name is None:
                    mreg_name = region.setup_mirror_region(mreg_name,
                                                           ret_name=True)
                    self.arg_trace_regions[aname] = mreg_name
                else:
                    region.setup_mirror_region(mreg_name)

        vals = []
        aux_pvars = []
        for svar in svars:
            # Allow only true state variables.
            if not svar.is_state():
                aux_pvars.append(svar)
                continue

            field = svar.get_field()
            is_trace = self.arg_traces[svar.name]
            trace_region = self.arg_trace_regions[svar.name]

            if svar.name in tgs:
                ps_tg = tgs[svar.name]
            else:
                ps_tg = v_tg

            val = ConnInfo(virtual=vvar,
                           state=svar,
                           primary=svar,
                           has_virtual=True,
                           has_state=True,
                           is_trace=is_trace,
                           trace_region=trace_region,
                           dc_type=dc_type,
                           v_tg=v_tg,
                           ps_tg=ps_tg,
                           region=region,
                           all_vars=all_vars)
            vals.append(val)

        pvars += aux_pvars
        for pvar in pvars:
            field = pvar.get_field()
            is_trace = self.arg_traces[pvar.name]
            trace_region = self.arg_trace_regions[pvar.name]

            if pvar.name in tgs:
                ps_tg = tgs[pvar.name]
            else:
                ps_tg = v_tg

            val = ConnInfo(virtual=vvar,
                           state=None,
                           primary=pvar.get_primary(),
                           has_virtual=vvar is not None,
                           has_state=False,
                           is_trace=is_trace,
                           trace_region=trace_region,
                           dc_type=dc_type,
                           v_tg=v_tg,
                           ps_tg=ps_tg,
                           region=region,
                           all_vars=all_vars)
            vals.append(val)

        if vvar and (len(vals) == 0):
            # No state, parameter variables, just the virtual one.
            val = ConnInfo(virtual=vvar,
                           state=vvar.get_primary(),
                           primary=vvar.get_primary(),
                           has_virtual=True,
                           has_state=False,
                           is_trace=False,
                           trace_region=None,
                           dc_type=dc_type,
                           v_tg=v_tg,
                           ps_tg=v_tg,
                           region=region,
                           all_vars=all_vars)
            vals.append(val)

        return vals

    def get_args_by_name(self, arg_names):
        """
        Return arguments by name.
        """
        out = []
        for name in arg_names:
            try:
                ii = self.arg_names.index(name)

            except ValueError:
                raise ValueError('non-existing argument! (%s)' % name)

            out.append(self.args[ii])

        return out

    def get_args(self, arg_types=None, **kwargs):
        """
        Return arguments by type as specified in arg_types (or
        self.ats). Arguments in **kwargs can override the ones assigned
        at the term construction - this is useful for passing user data.
        """
        ats = self.ats
        if arg_types is None:
            arg_types = ats
        args = []

        region_name, iorder = self.region.name, self.integral.order
        for at in arg_types:
            ii = ats.index(at)
            arg_name = self.arg_names[ii]
            if isinstance(arg_name, basestr):
                if arg_name in kwargs:
                    args.append(kwargs[arg_name])

                else:
                    args.append(self.args[ii])

            else:
                mat, par_name = self.args[ii]
                if mat is not None:
                    mat_data = mat.get_data((region_name, iorder), par_name)

                else:
                    mat_data = None

                args.append(mat_data)

        return args

    def get_kwargs(self, keys, **kwargs):
        """Extract arguments from **kwargs listed in keys (default is
        None)."""
        return [kwargs.get(name) for name in keys]

    def get_arg_name(self, arg_type, full=False, join=None):
        """
        Get the name of the argument specified by `arg_type.`

        Parameters
        ----------
        arg_type : str
            The argument type string.
        full : bool
            If True, return the full name. For example, if the name of a
            variable argument is 'u' and its time derivative is
            requested, the full name is 'du/dt'.
        join : str, optional
            Optionally, the material argument name tuple can be joined
            to a single string using the `join` string.

        Returns
        -------
        name : str
            The argument name.
        """
        try:
            ii = self.ats.index(arg_type)
        except ValueError:
            return None

        name = self.arg_names[ii]
        if full:
            # Include derivatives.
            if self.arg_derivatives[name]:
                name = 'd%s/%s' % (name, self.arg_derivatives[name])

        if (join is not None) and isinstance(name, tuple):
            name = join.join(name)

        return name

    def setup_integration(self):
        self.has_geometry = True

        self.geometry_types = {}
        if isinstance(self.integration, basestr):
            for var in self.get_variables():
                self.geometry_types[var.name] = self.integration

        else:
            if self.mode is not None:
                self.integration = self._integration[self.mode]

            if self.integration is not None:
                for arg_type, gtype in six.iteritems(self.integration):
                    var = self.get_args(arg_types=[arg_type])[0]
                    self.geometry_types[var.name] = gtype

        gtypes = list(set(self.geometry_types.values()))

        if 'surface_extra' in gtypes:
            self.dof_conn_type = 'volume'

        elif len(gtypes):
            self.dof_conn_type = gtypes[0]

    def get_region(self):
        return self.region

    def get_geometry_types(self):
        """
        Returns
        -------
        out : dict
            The required geometry types for each variable argument.
        """
        return self.geometry_types

    def get_dof_conn_type(self):
        return Struct(name='dof_conn_info', type=self.dof_conn_type,
                      region_name=self.region.name)

    def get_assembling_cells(self, shape=None):
        """
        Return the assembling cell indices into a DOF connectivity.
        """
        cells = nm.arange(shape[0], dtype=nm.int32)

        return cells

    def time_update(self, ts):
        if ts is not None:
            self.step = ts.step
            self.dt = ts.dt
            self.is_quasistatic = ts.is_quasistatic

        if 'ts' in self._kwargs:
            self._kwargs['ts'].update(ts)

    def advance(self, ts):
        """
        Advance to the next time step. Implemented in subclasses.
        """

    def get_vector(self, variable):
        """Get the vector stored in `variable` according to self.arg_steps
        and self.arg_derivatives. Supports only the backward difference w.r.t.
        time."""

        name = variable.name
        return variable(step=self.arg_steps[name],
                        derivative=self.arg_derivatives[name])

    def get_variables(self, as_list=True):

        if as_list:
            variables = self.get_args_by_name(self.names.variable)

        else:
            variables = {}
            for var in self.get_args_by_name(self.names.variable):
                variables[var.name] = var

        return variables

    def get_virtual_variable(self):
        aux = self.get_args_by_name(self.names.virtual)
        if len(aux) == 1:
            var = aux[0]

        else:
            var = None

        return var

    def get_state_variables(self, unknown_only=False):
        variables = self.get_args_by_name(self.names.state)

        if unknown_only:
            variables = [var for var in variables
                         if (var.kind == 'unknown') and
                         (self.arg_steps[var.name] == 0)]

        return variables

    def get_parameter_variables(self):
        return self.get_args_by_name(self.names.parameter)

    def get_materials(self, join=False):
        materials = self.get_args_by_name(self.names.material)

        for mat in materials:
            if mat[0] is None:
                materials.remove(mat)

        if join:
            materials = list(set(mat[0] for mat in materials))

        return materials

    def get_qp_key(self):
        """
        Return a key identifying uniquely the term quadrature points.
        """
        return (self.region.name, self.integral.order)

    def get_physical_qps(self):
        """
        Get physical quadrature points corresponding to the term region
        and integral.
        """
        from sfepy.discrete.common.mappings import get_physical_qps, PhysicalQPs

        if self.integration == 'point':
            phys_qps = PhysicalQPs()

        else:
            phys_qps = get_physical_qps(self.region, self.integral)

        return phys_qps

    def get_mapping(self, variable, get_saved=False, return_key=False):
        """
        Get the reference mapping from a variable.

        Notes
        -----
        This is a convenience wrapper of Field.get_mapping() that
        initializes the arguments using the term data.
        """
        integration = self.geometry_types[variable.name]
        is_trace = self.arg_traces[variable.name]

        if is_trace:
            mreg_name = self.arg_trace_regions[variable.name]
            region = self.region.get_mirror_region(mreg_name)

        else:
            region = self.region

        out = variable.field.get_mapping(region,
                                         self.integral, integration,
                                         get_saved=get_saved,
                                         return_key=return_key)

        return out

    def get_data_shape(self, variable):
        """
        Get data shape information from variable.

        Notes
        -----
        This is a convenience wrapper of FieldVariable.get_data_shape() that
        initializes the arguments using the term data.
        """
        integration = self.geometry_types[variable.name]
        is_trace = self.arg_traces[variable.name]

        if is_trace:
            mreg_name = self.arg_trace_regions[variable.name]
            region = self.region.get_mirror_region(mreg_name)

        else:
            region = self.region

        out = variable.get_data_shape(self.integral, integration, region.name)
        return out

    def get(self, variable, quantity_name, bf=None, integration=None,
            step=None, time_derivative=None):
        """
        Get the named quantity related to the variable.

        Notes
        -----
        This is a convenience wrapper of Variable.evaluate() that
        initializes the arguments using the term data.
        """
        name = variable.name

        step = get_default(step, self.arg_steps[name])
        time_derivative = get_default(time_derivative,
                                      self.arg_derivatives[name])
        integration = get_default(integration, self.geometry_types[name])

        data = variable.evaluate(mode=quantity_name,
                                 region=self.region, integral=self.integral,
                                 integration=integration,
                                 step=step, time_derivative=time_derivative,
                                 is_trace=self.arg_traces[name], bf=bf,
                                 trace_region=self.arg_trace_regions[name])
        return data

    def check_shapes(self, *args, **kwargs):
        """
        Check term argument shapes at run-time.
        """
        from sfepy.base.base import output
        from sfepy.mechanics.tensors import dim2sym

        dim = self.region.dim
        sym = dim2sym(dim)

        def _parse_scalar_shape(sh):
            if isinstance(sh, basestr):
                if sh == 'D':
                    return dim

                elif sh == 'D2':
                    return dim**2

                elif sh == 'S':
                    return sym

                elif sh == 'N': # General number.
                    return nm.inf

                elif sh == 'str':
                    return 'str'

                else:
                    return int(sh)

            else:
                return sh

        def _parse_tuple_shape(sh):
            if isinstance(sh, basestr):
                return tuple((_parse_scalar_shape(ii.strip())
                              for ii in sh.split(',')))

            else:
                return (int(sh),)

        arg_kinds = get_arg_kinds(self.ats)

        arg_shapes_list = self.arg_shapes
        if not isinstance(arg_shapes_list, list):
            arg_shapes_list = [arg_shapes_list]

        # Loop allowed shapes until a match is found, else error.
        allowed_shapes = []
        prev_shapes = {}
        actual_shapes = {}
        for _arg_shapes in arg_shapes_list:
            # Unset shapes are taken from the previous iteration.
            arg_shapes = copy(prev_shapes)
            arg_shapes.update(_arg_shapes)
            prev_shapes = arg_shapes

            allowed_shapes.append(arg_shapes)

            n_ok = 0
            for ii, arg_kind in enumerate(arg_kinds):
                if arg_kind in ('user', 'ts'):
                    n_ok += 1
                    continue

                arg = args[ii]
                key = '%s:%s' % (self.ats[ii], self.arg_names[ii])

                if self.mode is not None:
                    extended_ats = self.ats[ii] + ('/%s' % self.mode)

                else:
                    extended_ats = self.ats[ii]

                try:
                    sh = arg_shapes[self.ats[ii]]

                except KeyError:
                    sh = arg_shapes[extended_ats]

                if arg_kind.endswith('variable'):
                    n_el, n_qp, _dim, n_en, n_c = self.get_data_shape(arg)
                    actual_shapes[key] = (n_c,)
                    shape = _parse_scalar_shape(sh[0] if isinstance(sh, tuple)
                                                else sh)
                    if nm.isinf(shape):
                        n_ok += 1

                    else:
                        n_ok += shape == n_c

                elif arg_kind.endswith('material'):
                    if arg is None: # Switched-off opt_material.
                        n_ok += sh is None
                        continue

                    if sh is None:
                        continue

                    prefix = ''
                    if isinstance(sh, basestr):
                        aux = tuple(ii.strip() for ii in sh.split(':'))
                        if len(aux) == 2:
                            prefix, sh = aux

                    if sh == 'str':
                        n_ok += isinstance(arg, basestr)
                        continue

                    shape = _parse_tuple_shape(sh)
                    ls = len(shape)

                    aarg = nm.array(arg, ndmin=1)
                    actual_shapes[key] = aarg.shape

                    # Substiture general dimension 'N' with actual value.
                    iinfs = nm.where(nm.isinf(shape))[0]
                    if len(iinfs):
                        shape = list(shape)
                        for iinf in iinfs:
                            shape[iinf] = aarg.shape[-ls+iinf]
                        shape = tuple(shape)

                    if (ls > 1) or (shape[0] > 1):
                        # Array.
                        n_ok += shape == aarg.shape[-ls:]
                        actual_shapes[key] = aarg.shape[-ls:]

                    elif (ls == 1) and (shape[0] == 1):
                        # Scalar constant or callable as term argument
                        from numbers import Number
                        n_ok += isinstance(arg, Number) or callable(arg)

                else:
                    n_ok += 1

            if n_ok == len(arg_kinds):
                break

        else:
            term_str = self.get_str()
            output('allowed argument shapes for term "%s":' % term_str)
            output(allowed_shapes)
            output('actual argument shapes:')
            output(actual_shapes)
            raise ValueError('wrong arguments shapes for "%s" term! (see above)'
                             % term_str)

    def standalone_setup(self):
        from sfepy.discrete import create_adof_conns, Variables

        conn_info = {'aux' : self.get_conn_info()}
        adcs = create_adof_conns(conn_info, None)

        variables = Variables(self.get_variables())
        variables.set_adof_conns(adcs)

        materials = self.get_materials(join=True)

        for mat in materials:
            mat.time_update(None, [Struct(terms=[self])])

    def call_get_fargs(self, args, kwargs):
        try:
            fargs = self.get_fargs(*args, **kwargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        return fargs

    def call_function(self, out, fargs):
        try:
            status = self.function(out, *fargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        if status:
            terms.errclear()
            raise ValueError('term evaluation failed! (%s)' % self.name)

        return status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        out = nm.empty(shape, dtype=nm.float64)

        if mode == 'eval':
            status = self.call_function(out, fargs)
            # Sum over elements but not over components.
            out1 = nm.sum(out, 0).squeeze()
            return out1, status

        else:
            status = self.call_function(out, fargs)

            return out, status

    def eval_complex(self, shape, fargs, mode='eval', term_mode=None,
                     diff_var=None, **kwargs):
        rout = nm.empty(shape, dtype=nm.float64)

        fargsd = split_complex_args(fargs)

        # Assuming linear forms. Then the matrix is the
        # same both for real and imaginary part.
        rstatus = self.call_function(rout, fargsd['r'])
        if (diff_var is None) and len(fargsd) >= 2:
            iout = nm.empty(shape, dtype=nm.float64)
            istatus = self.call_function(iout, fargsd['i'])

            if mode == 'eval' and len(fargsd) >= 4:
                irout = nm.empty(shape, dtype=nm.float64)
                irstatus = self.call_function(irout, fargsd['ir'])
                riout = nm.empty(shape, dtype=nm.float64)
                ristatus = self.call_function(riout, fargsd['ri'])

                out = (rout - iout) + (riout + irout) * 1j
                status = rstatus or istatus or ristatus or irstatus

            else:
                out = rout + 1j * iout
                status = rstatus or istatus

        else:
            out, status = rout + 0j, rstatus

        if mode == 'eval':
            out1 = nm.sum(out, 0).squeeze()
            return out1, status

        else:
            return out, status

    def evaluate(self, mode='eval', diff_var=None,
                 standalone=True, ret_status=False, **kwargs):
        """
        Evaluate the term.

        Parameters
        ----------
        mode : 'eval' (default), or 'weak'
            The term evaluation mode.

        Returns
        -------
        val : float or array
            In 'eval' mode, the term returns a single value (the
            integral, it does not need to be a scalar), while in 'weak'
            mode it returns an array for each element.
        status : int, optional
            The flag indicating evaluation success (0) or failure
            (nonzero). Only provided if `ret_status` is True.
        iels : array of ints, optional
            The local elements indices in 'weak' mode. Only provided in
            non-'eval' modes.
        """
        if standalone:
            self.standalone_setup()

        kwargs = kwargs.copy()
        term_mode = kwargs.pop('term_mode', None)

        if mode in ('eval', 'el_eval', 'el_avg', 'qp'):
            args = self.get_args(**kwargs)
            self.check_shapes(*args)

            emode = 'eval' if mode == 'el_eval' else mode
            _args = tuple(args) + (emode, term_mode, diff_var)

            shape, dtype = self.get_eval_shape(*_args, **kwargs)

            if shape[0] == 0:
                val = nm.zeros(shape, dtype=dtype)
                status = 0

            else:
                fargs = self.call_get_fargs(_args, kwargs)

                if dtype == nm.float64:
                    val, status = self.eval_real(shape, fargs, mode,
                                                 term_mode,
                                                 **kwargs)

                elif dtype == nm.complex128:
                    val, status = self.eval_complex(shape, fargs, mode,
                                                    term_mode,
                                                    **kwargs)

                else:
                    raise ValueError('unsupported term dtype! (%s)' % dtype)

            val *= self.sign
            out = (val,)

        elif mode == 'weak':
            varr = self.get_virtual_variable()
            if varr is None:
                raise ValueError('no virtual variable in weak mode! (in "%s")'
                                 % self.get_str())

            if diff_var is not None:
                varc = self.get_variables(as_list=False)[diff_var]

            args = self.get_args(**kwargs)
            self.check_shapes(*args)

            n_elr, n_qpr, dim, n_enr, n_cr = self.get_data_shape(varr)
            n_row = n_cr * n_enr

            if diff_var is None:
                shape = (n_elr, 1, n_row, 1)

            else:
                n_elc, n_qpc, dim, n_enc, n_cc = self.get_data_shape(varc)
                n_col = n_cc * n_enc

                shape = (n_elr, 1, n_row, n_col)

            if shape[0] == 0:
                vals = nm.zeros(shape, dtype=varr.dtype)
                status = 0

            else:
                _args = tuple(args) + (mode, term_mode, diff_var)
                fargs = self.call_get_fargs(_args, kwargs)

                if varr.dtype == nm.float64:
                    vals, status = self.eval_real(shape, fargs, mode,
                                                  term_mode,
                                                  diff_var, **kwargs)

                elif varr.dtype == nm.complex128:
                    vals, status = self.eval_complex(shape, fargs, mode,
                                                     term_mode,
                                                     diff_var, **kwargs)

                else:
                    raise ValueError('unsupported term dtype! (%s)'
                                     % varr.dtype)

            if not isinstance(vals, tuple):
                vals *= self.sign
                iels = self.get_assembling_cells(vals.shape)

            else:
                vals = (self.sign * vals[0],) + vals[1:]
                iels = None

            out = (vals, iels)

        if goptions['check_term_finiteness']:
            assert_(nm.isfinite(out[0]).all(),
                    msg='"%s" term values not finite!' % self.get_str())

        if ret_status:
            out = out + (status,)

        if len(out) == 1:
            out = out[0]

        return out

    def assemble_to(self, asm_obj, val, iels, mode='vector', diff_var=None):
        """
        Assemble the results of term evaluation.

        For standard terms, assemble the values in `val` corresponding to
        elements/cells `iels` into a vector or a CSR sparse matrix `asm_obj`,
        depending on `mode`.

        For terms with a dynamic connectivity (e.g. contact terms), in
        `'matrix'` mode, return the extra COO sparse matrix instead. The extra
        matrix has to be added to the global matrix by the caller. By default,
        this is done in :func:`Equations.evaluate()
        <sfepy.discrete.equations.Equations.evaluate()>`.
        """
        import sfepy.discrete.common.extmods.assemble as asm

        vvar = self.get_virtual_variable()
        dc_type = self.get_dof_conn_type()

        extra = None

        if mode == 'vector':
            if asm_obj.dtype == nm.float64:
                assemble = asm.assemble_vector

            else:
                assert_(asm_obj.dtype == nm.complex128)
                assemble = asm.assemble_vector_complex
                for ii in range(len(val)):
                    if not(val[ii].dtype == nm.complex128):
                        val[ii] = nm.complex128(val[ii])

            if not isinstance(val, tuple):
                dc = vvar.get_dof_conn(dc_type)
                assert_(val.shape[2] == dc.shape[1])

                assemble(asm_obj, val, iels, 1.0, dc)

            else:
                vals, rows, var = val
                if var.eq_map is not None:
                    eq = var.eq_map.eq

                    rows = eq[rows]
                    active = (rows >= 0)
                    vals, rows = vals[active], rows[active]

                # Assumes no repeated indices in rows!
                asm_obj[rows] += vals

        elif mode == 'matrix':
            if asm_obj.dtype == nm.float64:
                assemble = asm.assemble_matrix

            else:
                assert_(asm_obj.dtype == nm.complex128)
                assemble = asm.assemble_matrix_complex

            svar = diff_var
            tmd = (asm_obj.data, asm_obj.indptr, asm_obj.indices)

            if ((asm_obj.dtype == nm.complex128)
                and (val.dtype == nm.float64)):
                val = val.astype(nm.complex128)

            sign = 1.0
            if self.arg_derivatives[svar.name]:
                if not self.is_quasistatic or (self.step > 0):
                    sign *= 1.0 / self.dt

                else:
                    sign = 0.0

            if not isinstance(val, tuple):
                rdc = vvar.get_dof_conn(dc_type)

                is_trace = self.arg_traces[svar.name]
                trace_region = self.arg_trace_regions[svar.name]
                cdc = svar.get_dof_conn(dc_type, is_trace, trace_region)
                assert_(val.shape[2:] == (rdc.shape[1], cdc.shape[1]))

                assemble(tmd[0], tmd[1], tmd[2], val, iels, sign, rdc, cdc)

            else:
                from scipy.sparse import coo_matrix

                vals, rows, cols, rvar, cvar = val
                if rvar.eq_map is not None:
                    req, ceq = rvar.eq_map.eq, cvar.eq_map.eq

                    rows, cols = req[rows], ceq[cols]
                    active = (rows >= 0) & (cols >= 0)
                    vals, rows, cols = vals[active], rows[active], cols[active]

                extra = coo_matrix((sign * vals, (rows, cols)),
                                   shape=asm_obj.shape)

        else:
            raise ValueError('unknown assembling mode! (%s)' % mode)

        return extra
