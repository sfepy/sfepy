from sfepy.base.base import *
try:
    from sfepy.terms.extmods import terms
except (ImportError, AttributeError):
    from sfepy.base.base import output
    msg = 'sfepy extension modules are not compiled!\ntype "make"'
    raise ImportError( msg )
from sfepy.base.la import split_range, combine
#from sfepy.base.ioutils import read_cache_data, write_cache_data

_match_var = re.compile( '^virtual$|^state(_[_a-zA-Z0-9]+)?$'\
                        + '|^parameter(_[_a-zA-Z0-9]+)?$' ).match
_match_state = re.compile( '^state(_[_a-zA-Z0-9]+)?$' ).match
_match_parameter = re.compile( '^parameter(_[_a-zA-Z0-9]+)?$' ).match
_match_material = re.compile( '^material(_[_a-zA-Z0-9]+)?$' ).match
_match_material_root = re.compile( '(.+)\.(.*)' ).match

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

##
# 24.07.2006, c
class Terms( Container ):

    ##
    # 24.07.2006, c
    def classify_args( self, variables ):
        for term in self:
            term.classify_args( variables )
            
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

Volume = 'Volume'
Surface = 'Surface'
Edge = 'Edge'
Point = 'Point'
SurfaceExtra = 'SurfaceExtra'

##
# 21.07.2006, c
class Term( Struct ):
    name = ''
    arg_types = ()
    geometry = []

    @staticmethod
    def from_desc(constructor, desc, regions):
        try:
            region = regions[desc.region]
        except IndexError:
            raise KeyError('region "%s" does not exist!' % desc.region)
        obj = constructor(desc.name, desc.sign,
                          region=region, integral_name=desc.integral)

        arg_names = []
        arg_steps = {}
        arg_derivatives = {}
        arg_traces = {}
        for arg in desc.args:
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

            arg_names.append( name )
            arg_steps[name] = step
            arg_derivatives[name] = derivative
            arg_traces[name] = trace

        obj.arg_names = arg_names
        obj.arg_steps = arg_steps
        obj.arg_derivatives = arg_derivatives
        obj.arg_traces = arg_traces

        return obj

    def __init__(self, name, sign, region=None, integral_name=None,
                 dof_conn_type='volume', function=None):
        self.name = name
        self.sign = sign
        self.ats = list(self.arg_types)

        self.char_fun = CharacteristicFunction(region)
        self.region = region
        self.integral_name = integral_name
        self.dof_conn_type = dof_conn_type
        self.function = function
        self.step = 0
        self.dt = 1.0
        self.has_integral = True
        self.has_region = True
        self.has_geometry = True
        
        self.itype = itype = None
        aux = re.compile('([a-z]+)_.*').match(name)
        if aux:
            itype = aux.group(1)
        self.raw_itype = itype
    
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """Subclasses either implement __call__ or plug in a proper _call()."""
        return self._call( diff_var, chunk_size, **kwargs )

    def _call( self, diff_var = None, chunk_size = None, **kwargs ):
        msg = 'base class method "_call" called for %s' % self.__class__.__name__
        raise RuntimeError( msg )
    
    ##
    # 16.11.2005, c
    def get_arg_names( self ):
        return self.__arg_names

    def set_arg_names( self, val ):
        if isinstance( self.__class__.arg_types[0], tuple ):
            n_arg = len( self.__class__.arg_types[0] )
        else:
            n_arg = len( self.__class__.arg_types )

        if len( val ) != n_arg:
            raise ValueError( 'equal shapes: %s, %s' \
                              % (val, self.__class__.arg_types) )
        self.__arg_names = val
    arg_names = property( get_arg_names, set_arg_names )

    ##
    # 24.07.2006, c
    def classify_args(self, variables):
        """state variable can be in place of parameter variable and vice
        versa."""
        self.names = Struct( name = 'arg_names',
                             material = [], variable = [], user = [],
                             state = [], virtual = [], parameter = [],
                             material_split = [] )

        msg = "variable '%s' requested by term '%s' does not exist!"

        if isinstance( self.arg_types[0], tuple ):
            assert_( len( self.modes ) == len( self.arg_types )\
                     == len( self.__class__.geometry ) )
            # Find matching call signature.
            matched = []
            for it, arg_types in enumerate( self.arg_types ):
                failed = False
                for ii, arg_type in enumerate( arg_types ):
                    name = self.__arg_names[ii]
                    if _match_var( arg_type ):
                        names = self.names.variable
                        try:
                            var = variables[name]
                        except KeyError:
                            raise KeyError( msg )

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
                self.geometry = self.__class__.geometry[i_match]
                self.mode = self.modes[i_match]
            elif len( matched ) == 0:
                msg = 'cannot match arguments! (%s)' % self.__arg_names
                raise ValueError( msg )
            else:
                msg = 'ambiguous arguments! (%s)' % self.__arg_names
                raise ValueError( msg )
        else:
            arg_types = self.arg_types
            self.mode = None

        # Set actual argument types.
        self.ats = list( arg_types )

        for ii, arg_type in enumerate( arg_types ):
            name = self.__arg_names[ii]
            if _match_var( arg_type ):
                names = self.names.variable
                try:
                    var = variables[name]
                except KeyError:
                    raise KeyError( msg )

                if _match_state( arg_type ) and \
                       var.is_state_or_parameter():
                    self.names.state.append( name )
                elif (arg_type == 'virtual') and var.is_virtual():
                    self.names.virtual.append( name )
                elif _match_parameter( arg_type ) and \
                         var.is_state_or_parameter():
                    self.names.parameter.append( name )

            elif _match_material( arg_type ):
                names = self.names.material
                match = _match_material_root( name )
                if match:
                    self.names.material_split.append( (match.group( 1 ),
                                                      match.group( 2 )) )
                else:
                    self.names.material_split.append( (name, None) )
            else:
                names = self.names.user
            names.append( name )

        self.n_virtual = len( self.names.virtual )
        if self.n_virtual > 1:
            raise ValueError( 'at most one virtial variable is allowed! (%d)'\
                              % self.n_virtual )

        self.set_arg_types()

        if (self.raw_itype == 'dw') and (self.mode == 'eval'):
            self.itype = 'd'
        else:
            self.itype = self.raw_itype

    def set_arg_types( self ):
        pass

    def check_args(self, variables, materials, user=None):
        """Common checking to all terms."""
        check_names(self.get_variable_names(), variables.names,
                    'variable(s) "%s" not found!')
        check_names(self.get_material_names(), materials.names,
                    'material(s) "%s" not found!')

        if isinstance(user, dict):
            check_names(self.get_user_names(), user.keys(),
                        'user data "%s" not found!')

        elif user is not None:
            raise ValueError('user data must be a dict or None!')

        igs = self.char_fun.igs
        vns = self.get_variable_names()
        for name in vns:
            field = variables[name].get_field()
            if field is None:
                continue

            if self.arg_traces[name]:
                if not nm.all(nm.setmember1d(self.region.all_vertices,
                                             field.region.all_vertices)):
                    msg = ('%s: incompatible regions: (self, trace of field %s)'
                           + '(%s in %s)') %\
                           (self.name, field.name,
                            self.region.all_vertices, field.region.all_vertices)
                    raise ValueError(msg)
            else:
                if not set( igs ).issubset( set( field.aps.igs ) ):
                    msg = ('%s: incompatible regions: (self, field)'
                           + ' (%s(%s) in %s(%s)') %\
                             (self.name, igs, name, field.igs(), field.name)
                    raise ValueError(msg)

        mns = self.get_material_names()
        for name in mns:
            mat = materials[name]

            if not set( igs ).issubset( set( mat.igs ) ):
                msg= ('%s: incompatible regions: (self, material)'
                      + ' (%s(%s) in %s(%s)') %\
                      (self.name, igs, name, mat.igs, mat.name)
                raise ValueError(msg)

    ##
    # 24.07.2006, c
    def get_variable_names( self ):
        return self.names.variable

    ##
    # 24.07.2006, c
    # 02.08.2006
    def get_material_names( self ):
        return [aux[0] for aux in self.names.material_split]

    ##
    # 24.07.2006, c
    def get_user_names( self ):
        return self.names.user

    ##
    # c: 24.07.2006, r: 04.07.2008
    def get_virtual_name( self, variables = None ):
        if not self.names.virtual:
            return None
        
        name = self.names.virtual[0]
        if variables is None:
            return name
        else:
            if variables[name].kind == 'test':
                return name
            else:
                msg = 'variable %s is not virtual!' % name
                raise ValueError( msg )

    ##
    # c: 24.07.2006, r: 04.07.2008
    def get_state_names( self, variables = None ):
        """If variables are given, return only true unknowns whose data are of
        the current time step (0)."""
        if variables is None:
            return copy( self.names.state )
        else:
            return [st for st in self.names.state
                    if (variables[st].kind == 'unknown') and
                    (self.arg_steps[st] == 0)]
            
    ##
    # 26.07.2007, c
    def get_parameter_names( self ):
        return copy( self.names.parameter )

    def get_conn_key(self):
        """The key to be used in DOF connectivity information."""
        key = (self.name,) + tuple(self.arg_names)
        key += (self.integral_name, self.region.name)

        return key

    def get_args( self, arg_types = None, **kwargs ):
        """Extract arguments from **kwargs by type as specified in arg_types
        (or self.ats)."""
        ats = self.ats
        if arg_types is None:
            arg_types = ats
        args = []

        iname, region_name, ig = self.get_current_group()
        for at in arg_types:
            ii = ats.index( at )
            name = self.arg_names[ii]
##             print at, ii, name
            if at[:8] == 'material':
##                 print self.names.material
                im = self.names.material.index( name )
                split = self.names.material_split[im]
                mat = kwargs[split[0]]
                mat_data = mat.get_data((region_name, self.integral_name),
                                        ig, split[1])
                args.append(mat_data)
            else:
                args.append( kwargs[name] )

        return args

    def get_kwargs( self, keys, **kwargs ):
        """Extract arguments from **kwargs listed in keys (default is
        None)."""
        return [kwargs.get( name ) for name in keys]

    ##
    # 24.07.2006, c
    def get_arg_name( self, arg_type, full = False ):
        try:
            ii = self.ats.index( arg_type )
        except ValueError:
            return None

        name = self.arg_names[ii]
        if full:
            # Include derivatives.
            if self.arg_derivatives[name]:
                name = 'd%s/%s' % (name, self.arg_derivatives[name] )

        return name

    ##
    # c: 29.11.2007, r: 10.04.2008
    def describe_geometry( self, geometries, variables, integrals ):
        """Takes reference to the used integral."""
        if not self.has_geometry: return

        try:
            integral = integrals[self.integral_name]
        except IndexError:
            msg = 'integral %s is not defined!' % self.integral_name
            raise ValueError(msg)
            
        integral.create_qp()
        self.integral = integral
        
        tgs = self.get_geometry()
        for var_name in self.get_variable_names():
            print '>>>>>', self.name, var_name

            variable = variables[var_name]
            if not variable.has_field: continue

            field = variable.field

            is_trace = self.arg_traces[variable.name]
            if not is_trace:
                assert_( field.region.contains( self.region ) )

##             print field.name, field.region_name
##             print field.bases

            if tgs.has_key( var_name ):
##                 print ':', tgs[var_name]

                if is_trace:
                    aux = variables.get_mirror_region(self.region,
                                                      return_ig_map=True)
                    region, ig_map = aux

                else:
                    region, ig_map = self.region, None
                    
                field.aps.describe_geometry(field, geometries,
                                            tgs[var_name], region, self.region,
                                            integral, ig_map=ig_map)

    def get_region(self):
        return self.region

    def get_geometry( self ):
        geom = self.geometry
        out = {}
        if geom:
            for (gtype, arg_type) in geom:
                arg_name = self.get_arg_name( arg_type )
                out[arg_name] = gtype
        return out

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
    # c: 27.02.2007, r: 15.04.2008
    def get_cache( self, base_name, ii ):
        args = self.use_caches[base_name][ii]
        ans = [self.get_arg_name( arg, full = True ) for arg in args
               if not type( arg ) == dict]
##         print args, ans
##         pause()
        cname = '_'.join( [base_name] + ans )
        return self.caches[cname]

    ##
    # 02.03.2007, c
    def igs( self ):
        return self.char_fun.igs

    def needs_local_chunk(self):
        """Returns a tuple of booleans telling whether the term requires local
        element numbers in an assembling chunk to pass to an element
        contribution function, and to an assembling function."""
        ret = [False, False]
        if self.dof_conn_type == 'surface':
            ret[1] = True
            
        return tuple(ret)

    ##
    # c: 05.12.2007, r: 15.01.2008
    def iter_groups( self ):
        if self.dof_conn_type == 'point':
            igs = self.igs()[0:1]
        else:
            igs = self.igs()

        for ig in igs:
            self.set_current_group( ig )
            yield ig

    def time_update( self, ts ):
        self.step = ts.step
        self.dt = ts.dt

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

    def get_approximation(self, variable, kind = 'Volume' ):
        is_trace = self.arg_traces[variable.name]
        key = self.get_current_group()
        out = variable.get_approximation(key, kind=kind, is_trace=is_trace)
        return out

##     def get_quadrature_orders(self):
##         """Curently, it just takes the order of the term integral."""
##         return self.integral.order
