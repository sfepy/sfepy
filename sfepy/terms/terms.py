import re
from copy import copy

import numpy as nm

from sfepy.base.base import (as_float_or_complex, get_default, assert_,
                             Container, Struct, basestr)
from sfepy.base.compat import in1d

# Used for imports in term files.
from sfepy.terms.extmods import terms

from sfepy.linalg import split_range
#from sfepy.base.ioutils import read_cache_data, write_cache_data

_match_args = re.compile('^([^\(\}]*)\((.*)\)$').match
_match_var = re.compile( '^virtual$|^state(_[_a-zA-Z0-9]+)?$'\
                        + '|^parameter(_[_a-zA-Z0-9]+)?$' ).match
_match_state = re.compile( '^state(_[_a-zA-Z0-9]+)?$' ).match
_match_parameter = re.compile( '^parameter(_[_a-zA-Z0-9]+)?$' ).match
_match_material = re.compile( '^material(_[_a-zA-Z0-9]+)?$' ).match
_match_material_opt = re.compile( '^opt_material(_[_a-zA-Z0-9]+)?$' ).match
_match_material_root = re.compile( '(.+)\.(.*)' ).match

def get_shape_kind(integration):
    """
    Get data shape kind for given integration type.
    """
    if integration == 'surface':
        shape_kind = 'surface'

    elif integration in ('volume', 'surface_extra'):
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

def vector_chunk_generator( total_size, chunk_size, shape_in,
                            zero = False, set_shape = True, dtype = nm.float64 ):
    if not chunk_size:
        chunk_size = total_size
    shape = list( shape_in )

    sizes = split_range( total_size, chunk_size )
    ii = nm.array( 0, dtype = nm.int32 )
    for size in sizes:
        chunk = nm.arange( size, dtype = nm.int32 ) + ii
        if set_shape:
            shape[0] = size
        if zero:
            out = nm.zeros( shape, dtype = dtype )
        else:
            out = nm.empty( shape, dtype = dtype )
        yield out, chunk
        ii += size

def create_arg_parser():
    from pyparsing import Combine, Literal, Word, delimitedList, Group, \
         StringStart, StringEnd, Optional, nums, alphas, alphanums

    inumber = Word("+-"+nums, nums)

    history = Optional(Literal('[').suppress() + inumber
                       + Literal(']').suppress(), default=0)("history")
    history.setParseAction(lambda str, loc, toks: int(toks[0]))

    variable = Group(Word(alphas, alphanums + '._') + history)

    derivative = Group(Literal('d') + variable\
                       + Literal('/').suppress() + Literal('dt'))

    trace = Group(Literal('tr') + Literal('(').suppress() + variable \
                  + Literal(')').suppress())
    
    generalized_var = derivative | trace | variable

    args = StringStart() + delimitedList(generalized_var) + StringEnd()

    return args

def reorder_dofs_on_mirror(adof, dc, mirror_dc):
    nadof = nm.zeros_like(adof)
    for ifc in range(len(dc)):
        fc = dc[ifc]
        mfc = mirror_dc[ifc]
        ndrange = range(len(fc));
        lmap = -nm.ones((len(fc),), dtype=nm.int32)
        for ind in ndrange:
            for j in ndrange:
                if mfc[ind] == fc[j]:
                    lmap[ind] = j
                    break;
        nadof[ifc,:] = adof[ifc,lmap]

    return nadof
##
# 22.01.2006, c
class CharacteristicFunction( Struct ):
    ##
    # c: 22.01.2006, r: 09.05.2008
    def __init__( self, region ):
        self.igs = region.igs
        self.region = region
        self.i_current = None
        self.local_chunk = None
        self.ig = None

    def __call__(self, chunk_size, shape_in, zero=False, set_shape=True,
                  ret_local_chunk=False, dtype=nm.float64):
        els = self.region.cells[self.ig]
        for out, chunk in vector_chunk_generator( els.shape[0], chunk_size,
                                                  shape_in, zero, set_shape,
                                                  dtype ):
            self.local_chunk = chunk

            if ret_local_chunk:
                yield out, chunk
            else:
                yield out, els[chunk]
                
        self.local_chunk = None

    ##
    # 11.08.2006, c
    # 27.02.2007
    def set_current_group( self, ig ):
        self.ig = ig
        
    ##
    # 05.09.2006, c
    def get_local_chunk( self ):
        return self.local_chunk

class ConnInfo(Struct):

    def get_region(self, can_trace=True):
        if self.is_trace and can_trace:
            return self.region.get_mirror_region()[0]
        else:
            return self.region

    def get_region_name(self, can_trace=True):
        if self.is_trace and can_trace:
            reg = self.region.get_mirror_region()[0]
        else:
            reg = self.region

        if reg is not None:
            return reg.name
        else:
            return None

    def iter_igs(self):
        if self.region is not None:
            for ig in self.region.igs:
                if self.virtual_igs is not None:
                    ir = self.virtual_igs.index(ig)
                    rig = self.virtual_igs[ir]
                else:
                    rig = None

                if not self.is_trace:
                    ii = ig
                else:
                    ig_map_i = self.region.get_mirror_region()[2]
                    ii = ig_map_i[ig]

                if self.state_igs is not None:
                    ic = self.state_igs.index(ii)
                    cig = self.state_igs[ic]
                else:
                    cig = None

                yield rig, cig

        else:
            yield None, None

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

    ##
    # 24.07.2006, c
    # 02.08.2006
    def get_variable_names( self ):
        out = []
        for term in self:
            out.extend( term.get_variable_names() )
        return list( set( out ) )

    ##
    # 24.07.2006, c
    # 02.08.2006
    def get_material_names( self ):
        out = []
        for term in self:
            out.extend( term.get_material_names() )
        return list( set( out ) )

    ##
    # 24.07.2006, c
    # 02.08.2006
    def get_user_names( self ):
        out = []
        for term in self:
            out.extend( term.get_user_names() )
        return list( set( out ) )

    ##
    # 11.08.2006, c
    def set_current_group( self, ig ):
        for term in self:
            term.char_fun.set_current_group( ig )

class Term(Struct):
    name = ''
    arg_types = ()
    integration = 'volume'

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
            msg = "term '%s' is not in %s" % (name,
                                              sorted(term_table.keys()))
            raise ValueError(msg)

        obj = constructor(name, arg_str, integral, region, **kwargs)
        return obj
        
    @staticmethod
    def from_desc(constructor, desc, region, integrals=None):
        from sfepy.fem import Integrals
        
        if integrals is None:
            integrals = Integrals()

        obj = constructor(desc.name, desc.args, None, region)
        obj.set_integral(integrals.get(desc.integral, obj.get_integral_info()))
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

    def set_integral(self, integral):
        """
        Set the term integral.
        """
        self.integral = integral
        if self.integral is not None:
            self.integral_name = self.integral.name

            kind = self.get_integral_info()
            if kind != integral.kind:
                msg = "integral kind for term %s must be '%s'! (is '%s')" \
                      % (self.name, kind, integral.kind)
                raise ValueError(msg)

    def setup(self):
        self.char_fun = CharacteristicFunction(self.region)
        self.function = Struct.get(self, 'function', None)

        self.step = 0
        self.dt = 1.0
        self.is_quasistatic = False
        self.has_integral = True
        self.has_region = True

        self.itype = itype = None
        aux = re.compile('([a-z]+)_.*').match(self.name)
        if aux:
            itype = aux.group(1)
        self.raw_itype = itype

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

        parser = create_arg_parser()
        self.arg_desc = parser.parseString(self.arg_str)

        for arg in self.arg_desc:
            trace = False
            derivative = None

            if isinstance(arg[1], int):
                name, step = arg

            else:
                kind = arg[0]
                name, step = arg[1]
                if kind == 'd':
                    derivative = arg[2]
                elif kind == 'tr':
                    trace = True

            match = _match_material_root(name)
            if match:
                name = (match.group(1), match.group(2))

            self.arg_names.append(name)
            self.arg_steps[name] = step
            self.arg_derivatives[name] = derivative
            self.arg_traces[name] = trace

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

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """Subclasses either implement __call__ or plug in a proper _call()."""
        return self._call( diff_var, chunk_size, **kwargs )

    def _call( self, diff_var = None, chunk_size = None, **kwargs ):
        msg = 'base class method "_call" called for %s' % self.__class__.__name__
        raise RuntimeError( msg )

    def assign_args(self, variables, materials, user=None):
        """
        Check term argument existence in variables, materials, user data
        and assign the arguments to terms. Also check compatibility of
        field and term subdomain lists (igs).
        """
        if user is None:
            user = {}

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
        self.names = Struct(name = 'arg_names',
                            material = [], variable = [], user = [],
                            state = [], virtual = [], parameter = [])

        msg = "variable '%s' requested by term '%s' does not exist!"

        # check for "opt_material"
        if isinstance(self.arg_types[0], tuple):
            arg_types = self.arg_types[0]
        else:
            arg_types = self.arg_types

        matched = 0
        for ii, arg_type in enumerate(arg_types):
            if _match_material_opt(arg_type):
                matched += 1
                if ii > 0:
                    msg = 'opt_material at position %d, must be at 0!' % ii
                    raise ValueError(msg)
                if not(isinstance(self.args[ii], tuple)):
                    self.args.insert(ii, (None, None))
                    self.arg_names.insert(ii, (None, None))

        if matched > 1:
            msg = 'only one opt_material allowed, %d given!' % matched
            raise ValueError(msg)

        if isinstance(self.arg_types[0], tuple):
            assert_(len(self.modes) == len(self.arg_types))
            # Find matching call signature.
            matched = []
            for it, arg_types in enumerate( self.arg_types ):
                failed = False
                for ii, arg_type in enumerate( arg_types ):
                    name = self.arg_names[ii]
                    if _match_var( arg_type ):
                        names = self.names.variable
                        var = self.args[ii]

                        if _match_state( arg_type ) and \
                               var.is_state_or_parameter():
                            pass
                        elif (arg_type == 'virtual') and var.is_virtual():
                            pass
                        elif _match_parameter( arg_type ) and \
                                 var.is_state_or_parameter():
                            pass
                        else:
                            failed = True
                            break

                if not failed:
                    matched.append( it )

            if len( matched ) == 1:
                i_match = matched[0]
                arg_types = self.arg_types[i_match]
                self.mode = self.modes[i_match]
            elif len( matched ) == 0:
                msg = 'cannot match arguments! (%s)' % self.arg_names
                raise ValueError( msg )
            else:
                msg = 'ambiguous arguments! (%s)' % self.arg_names
                raise ValueError( msg )
        else:
            arg_types = self.arg_types
            self.mode = None

        # Set actual argument types.
        self.ats = list( arg_types )

        for ii, arg_type in enumerate( arg_types ):
            name = self.arg_names[ii]
            if _match_var( arg_type ):
                names = self.names.variable
                var = self.args[ii]

                if _match_state( arg_type ) and \
                       var.is_state_or_parameter():
                    self.names.state.append( name )
                elif (arg_type == 'virtual') and var.is_virtual():
                    self.names.virtual.append( name )
                elif _match_parameter( arg_type ) and \
                         var.is_state_or_parameter():
                    self.names.parameter.append( name )

            elif _match_material( arg_type ) or \
                     _match_material_opt( arg_type ):
                names = self.names.material

            else:
                names = self.names.user
            names.append( name )

        self.n_virtual = len( self.names.virtual )
        if self.n_virtual > 1:
            raise ValueError( 'at most one virtial variable is allowed! (%d)'\
                              % self.n_virtual )

        self.set_arg_types()

        self.setup_integration()

        if (self.raw_itype == 'dw') and (self.mode == 'eval'):
            self.itype = 'd'
        else:
            self.itype = self.raw_itype

    def set_arg_types( self ):
        pass

    def check_args(self):
        """
        Common checking to all terms.

        Check compatibility of field and term subdomain lists (igs).
        """
        vns = self.get_variable_names()
        for name in vns:
            field = self._kwargs[name].get_field()
            if field is None:
                continue

            if not nm.all(in1d(self.region.all_vertices,
                               field.region.all_vertices)):
                msg = ('%s: incompatible regions: (self, field %s)'
                       + '(%s in %s)') %\
                       (self.name, field.name,
                        self.region.all_vertices, field.region.all_vertices)
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

    ##
    # 26.07.2007, c
    def get_parameter_names( self ):
        return copy( self.names.parameter )

    def get_conn_key(self):
        """The key to be used in DOF connectivity information."""
        key = (self.name,) + tuple(self.arg_names)
        key += (self.integral_name, self.region.name)

        return key

    def get_conn_info(self):
        vvar = self.get_virtual_variable()
        svars = self.get_state_variables()
        pvars = self.get_parameter_variables()

        all_vars = self.get_variables()

        dc_type = self.get_dof_conn_type()
        tgs = self.get_geometry_types()

        v_igs = v_tg = None
        if vvar is not None:
            field = vvar.get_field()
            if field is not None:
                v_igs = field.igs
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
            is_any_trace = reduce(lambda x, y: x or y,
                                  self.arg_traces.values())
            if is_any_trace:
                region.setup_mirror_region()

        vals = []
        aux_pvars = []
        for svar in svars:
            # Allow only true state variables.
            if not svar.is_state():
                aux_pvars.append(svar)
                continue

            field = svar.get_field()
            if field is not None:
                s_igs = field.igs
            else:
                s_igs = None
            is_trace = self.arg_traces[svar.name]

            if svar.name in tgs:
                ps_tg = tgs[svar.name]
            else:
                ps_tg = v_tg

            val = ConnInfo(virtual = vvar, virtual_igs = v_igs,
                           state = svar, state_igs = s_igs,
                           primary = svar, primary_igs = s_igs,
                           has_virtual = True,
                           has_state = True,
                           is_trace = is_trace,
                           dc_type = dc_type,
                           v_tg = v_tg,
                           ps_tg = ps_tg,
                           region = region,
                           all_vars = all_vars)
            vals.append(val)

        pvars += aux_pvars
        for pvar in pvars:
            field = pvar.get_field()
            if field is not None:
                p_igs = field.igs
            else:
                p_igs = None
            is_trace = self.arg_traces[pvar.name]

            if pvar.name in tgs:
                ps_tg = tgs[pvar.name]
            else:
                ps_tg = v_tg

            val = ConnInfo(virtual = vvar, virtual_igs = v_igs,
                           state = None, state_igs = [],
                           primary = pvar.get_primary(), primary_igs = p_igs,
                           has_virtual = vvar is not None,
                           has_state = False,
                           is_trace = is_trace,
                           dc_type = dc_type,
                           v_tg = v_tg,
                           ps_tg = ps_tg,
                           region = region,
                           all_vars = all_vars)
            vals.append(val)

        if vvar and (len(vals) == 0):
            # No state, parameter variables, just the virtual one.
            val = ConnInfo(virtual = vvar, virtual_igs = v_igs,
                           state = vvar.get_primary(), state_igs = v_igs,
                           primary = vvar.get_primary(), primary_igs = v_igs,
                           has_virtual = True,
                           has_state = False,
                           is_trace = False,
                           dc_type = dc_type,
                           v_tg = v_tg,
                           ps_tg = v_tg,
                           region = region,
                           all_vars = all_vars)
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

        iname, region_name, ig = self.get_current_group()
        for at in arg_types:
            ii = ats.index(at)
            arg_name = self.arg_names[ii]
            ## print at, ii, arg_name
            if isinstance(arg_name, basestr):
                if arg_name in kwargs:
                    args.append(kwargs[arg_name])

                else:
                    args.append(self.args[ii])

            else:
                ## print self.names.material
                mat, par_name = self.args[ii]
                if mat is not None:
                    mat_data = mat.get_data((region_name, self.integral_name),
                                            ig, par_name)
                else:
                    mat_data = None

                args.append(mat_data)

        return args

    def get_kwargs( self, keys, **kwargs ):
        """Extract arguments from **kwargs listed in keys (default is
        None)."""
        return [kwargs.get( name ) for name in keys]

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

    def get_integral_info(self):
        """
        Get information on the term integral.

        Returns
        -------
        kind : 'v' or 's'
            The integral kind.
        """
        if self.integration:
            if self.integration == 'volume':
                kind = 'v'

            elif 'surface' in self.integration:
                kind = 's'

            elif self.integration == 'point':
                kind = None

            else:
                raise ValueError('unsupported term integration! (%s)'
                                 % self.integration)

        else:
            kind = None

        return kind

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
                for arg_type, gtype in self.integration.iteritems():
                    var = self.get_args(arg_types=[arg_type])[0]
                    self.geometry_types[var.name] = gtype

        gtypes = list(set(self.geometry_types.itervalues()))

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

    def get_current_group(self):
        return (self.integral_name, self.region.name, self.char_fun.ig)


    def get_dof_conn_type(self):
        return Struct(name='dof_conn_info', type=self.dof_conn_type,
                      region_name=self.region.name)

    ##
    # c: 16.02.2007, r: 15.01.2008
    def set_current_group( self, ig ):
        self.char_fun.set_current_group( ig )

    ##
    # 02.03.2007, c
    def igs( self ):
        return self.char_fun.igs

    def get_assembling_cells(self, shape=None):
        """
        According to the term integration type, return either the term
        region cell indices or local index sequence.
        """
        shape_kind = get_shape_kind(self.integration)
        ig = self.char_fun.ig

        cells = self.region.cells[ig]
        if shape_kind == 'surface':
            cells = nm.arange(cells.shape[0], dtype=nm.int32)

        elif shape_kind == 'point':
            cells = nm.arange(shape[0], dtype=nm.int32)

        return cells

    ##
    # c: 05.12.2007, r: 15.01.2008
    def iter_groups( self ):
        if self.dof_conn_type == 'point':
            igs = self.igs()[0:1]
        else:
            igs = self.igs()

        for ig in igs:
            if self.integration == 'volume':
                if not len(self.region.cells[ig]): continue
            self.set_current_group( ig )
            yield ig

    def time_update( self, ts ):
        if ts is not None:
            self.step = ts.step
            self.dt = ts.dt
            self.is_quasistatic = ts.is_quasistatic

    def advance(self, ts):
        """Advance to the next time step. Implemented in subclasses."""
        pass

    def get_vector( self, variable ):
        """Get the vector stored in `variable` according to self.arg_steps
        and self.arg_derivatives. Supports only the backward difference w.r.t.
        time."""

        name = variable.name
        return variable( step = self.arg_steps[name],
                         derivative = self.arg_derivatives[name] )

    def get_approximation(self, variable, get_saved=False):
        """
        Return approximation corresponding to `variable`. Also return
        the corresponding geometry (actual or saved, according to
        `get_saved`).
        """
        geo, _, key = self.get_mapping(variable, get_saved=get_saved,
                                       return_key=True)
        ig = key[2]
        ap = variable.get_approximation(ig)

        return ap, geo

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
        return (self.region.name, self.integral.name)

    def get_physical_qps(self):
        """
        Get physical quadrature points corresponding to the term region
        and integral.
        """
        from sfepy.fem.mappings import get_physical_qps

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
            region, ig_map, ig_map_i = self.region.get_mirror_region()
            ig = ig_map_i[self.char_fun.ig]

        else:
            region = self.region
            ig = self.char_fun.ig

        out = variable.field.get_mapping(ig, region,
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
            region, ig_map, ig_map_i = self.region.get_mirror_region()
            ig = ig_map_i[self.char_fun.ig]

        else:
            region = self.region
            ig = self.char_fun.ig

        out = variable.get_data_shape(ig, self.integral,
                                      integration, region.name)
        return out

    def get(self, variable, quantity_name, bf=None,
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

        data = variable.evaluate(self.char_fun.ig, mode=quantity_name,
                                 region=self.region, integral=self.integral,
                                 integration=self.geometry_types[name],
                                 step=step, time_derivative=time_derivative,
                                 is_trace=self.arg_traces[name], bf=bf)
        return data

    def check_shapes(self, *args, **kwargs):
        """
        Default implementation of function to check term argument shapes
        at run-time.
        """
        pass

    def standalone_setup(self):
        from sfepy.fem import setup_dof_conns

        conn_info = {'aux' : self.get_conn_info()}

        setup_dof_conns(conn_info)

        materials = self.get_materials(join=True)
        ## print materials

        for mat in materials:
            mat.time_update(None, [Struct(terms=[self])])

    def call_get_fargs(self, args, kwargs):
        try:
            fargs = self.get_fargs(*args, **kwargs)

        except RuntimeError:
            terms.errclear()
            raise ValueError

        return fargs

    def call_function(self, out, fargs):
        try:
            status = self.function(out, *fargs)

        except RuntimeError:
            terms.errclear()
            raise ValueError

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
            out, status = rout, rstatus

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

        if mode == 'eval':
            val = 0.0
            status = 0
            for ig in self.iter_groups():
                args = self.get_args(**kwargs)
                self.check_shapes(*args)

                _args = tuple(args) + (mode, term_mode, diff_var)
                fargs = self.call_get_fargs(_args, kwargs)

                shape, dtype = self.get_eval_shape(*_args, **kwargs)

                if dtype == nm.float64:
                    _v, stat = self.eval_real(shape, fargs, mode, term_mode,
                                               **kwargs)

                elif dtype == nm.complex128:
                    _v, stat = self.eval_complex(shape, fargs, mode, term_mode,
                                                 **kwargs)

                else:
                    raise ValueError('unsupported term dtype! (%s)' % dtype)

                val += _v
                status += stat

            val *= self.sign

        elif mode in ('el_avg', 'el', 'qp'):
            vals = None
            iels = nm.empty((0, 2), dtype=nm.int32)
            status = 0
            for ig in self.iter_groups():
                args = self.get_args(**kwargs)
                self.check_shapes(*args)

                _args = tuple(args) + (mode, term_mode, diff_var)
                fargs = self.call_get_fargs(_args, kwargs)

                shape, dtype = self.get_eval_shape(*_args, **kwargs)

                if dtype == nm.float64:
                    val, stat = self.eval_real(shape, fargs, mode, term_mode,
                                               **kwargs)

                elif dtype == nm.complex128:
                    val, stat = self.eval_complex(shape, fargs, mode, term_mode,
                                                  **kwargs)

                if vals is None:
                    vals = val

                else:
                    vals = nm.r_[vals, val]

                _iels = self.get_assembling_cells(val.shape)
                aux = nm.c_[nm.repeat(ig, _iels.shape[0])[:,None],
                            _iels[:,None]]
                iels = nm.r_[iels, aux]
                status += stat

            vals *= self.sign

        elif mode == 'weak':
            vals = []
            iels = []
            status = 0

            varr = self.get_virtual_variable()
            if diff_var is not None:
                varc = self.get_variables(as_list=False)[diff_var]

            for ig in self.iter_groups():
                args = self.get_args(**kwargs)
                self.check_shapes(*args)

                _args = tuple(args) + (mode, term_mode, diff_var)
                fargs = self.call_get_fargs(_args, kwargs)

                n_elr, n_qpr, dim, n_enr, n_cr = self.get_data_shape(varr)
                n_row = n_cr * n_enr

                if diff_var is None:
                    shape = (n_elr, 1, n_row, 1)

                else:
                    n_elc, n_qpc, dim, n_enc, n_cc = self.get_data_shape(varc)
                    n_col = n_cc * n_enc

                    shape = (n_elr, 1, n_row, n_col)

                if varr.dtype == nm.float64:
                    val, stat = self.eval_real(shape, fargs, mode, term_mode,
                                               diff_var, **kwargs)

                elif varr.dtype == nm.complex128:
                    val, stat = self.eval_complex(shape, fargs, mode, term_mode,
                                                  diff_var, **kwargs)

                else:
                    raise ValueError('unsupported term dtype! (%s)'
                                     % varr.dtype)

                vals.append(self.sign * val)
                iels.append((ig, self.get_assembling_cells(val.shape)))
                status += stat

        # Setup return value.
        if mode == 'eval':
            out = (val,)

        else:
            out = (vals, iels)

        if ret_status:
            out = out + (status,)

        if len(out) == 1:
            out = out[0]

        return out

    def assemble_to(self, asm_obj, val, iels, mode='vector', diff_var=None):
        import sfepy.fem.extmods.assemble as asm

        vvar = self.get_virtual_variable()
        dc_type = self.get_dof_conn_type()

        if mode == 'vector':
            if asm_obj.dtype == nm.float64:
                assemble = asm.assemble_vector

            else:
                assert_(asm_obj.dtype == nm.complex128)
                assemble = asm.assemble_vector_complex
                for ii in range(len(val)):
                    if not(val[ii].dtype == nm.complex128):
                        val[ii] = nm.complex128(val[ii])

            for ii, (ig, _iels) in enumerate(iels):
                vec_in_els = val[ii]
                dc = vvar.get_dof_conn(dc_type, ig, active=True)
                assert_(vec_in_els.shape[2] == dc.shape[1])

                assemble(asm_obj, vec_in_els, _iels, 1.0, dc)

        elif mode == 'matrix':
            if asm_obj.dtype == nm.float64:
                assemble = asm.assemble_matrix

            else:
                assert_(asm_obj.dtype == nm.complex128)
                assemble = asm.assemble_matrix_complex

            svar = diff_var
            tmd = (asm_obj.data, asm_obj.indptr, asm_obj.indices)
            for ii, (ig, _iels) in enumerate(iels):
                mtx_in_els = val[ii]
                if ((asm_obj.dtype == nm.complex128)
                    and (mtx_in_els.dtype == nm.float64)):
                    mtx_in_els = mtx_in_els.astype(nm.complex128)

                rdc = vvar.get_dof_conn(dc_type, ig, active=True)

                is_trace = self.arg_traces[svar.name]
                ## print dc_type, ig, is_trace
                cdc = svar.get_dof_conn(dc_type, ig, active=True,
                                        is_trace=is_trace)
                ## print svar.name, cdc.shape
                assert_(mtx_in_els.shape[2:] == (rdc.shape[1], cdc.shape[1]))

                sign = 1.0
                if self.arg_derivatives[svar.name]:
                    if not self.is_quasistatic or (self.step > 0):
                        sign *= 1.0 / self.dt

                    else:
                        sign = 0.0

                assemble(tmd[0], tmd[1], tmd[2], mtx_in_els,
                         _iels, sign, rdc, cdc)

        else:
            raise ValueError('unknown assembling mode! (%s)' % mode)

##     def get_quadrature_orders(self):
##         """Curently, it just takes the order of the term integral."""
##         return self.integral.order
