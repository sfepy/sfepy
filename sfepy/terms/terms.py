from sfepy.base.base import *
try:
    import extmods.terms as terms
except:
    output( 'warning: sfepy extension modules are not compiled!' )
    output( 'type "make"' )
from sfepy.base.la import split_range
#from sfepy.base.ioutils import read_cache_data, write_cache_data

_match_var = re.compile( '^virtual$|^state(_[a-zA-Z0-9]+)?$'\
                        + '|^parameter(_[a-zA-Z0-9]+)?$' ).match
_match_state = re.compile( '^state(_[a-zA-Z0-9]+)?$' ).match
_match_parameter = re.compile( '^parameter(_[a-zA-Z0-9]+)?$' ).match
_match_material = re.compile( '^material(_[a-zA-Z0-9]+)?$' ).match
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

    def __call__( self, chunk_size, shape_in, zero = False, set_shape = True,
                  dtype = nm.float64 ):
        els = self.region.cells[self.ig]
        for out, chunk in vector_chunk_generator( els.shape[0], chunk_size,
                                                  shape_in, zero, set_shape,
                                                  dtype ):
            self.local_chunk = chunk
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

    ##
    # 24.07.2006, c
    # 25.07.2006
    # 11.08.2006
    # 24.08.2006
    # 11.10.2006
    # 29.11.2006
    def __init__( self, region, name, sign, function = None ):
        self.char_fun  = CharacteristicFunction( region )
        self.region = region
        self.name = name
        self.sign = sign
        self.dof_conn_type = 'volume'
        self.function = function
        
        itype = None
        aux = re.compile( '([a-z]+)_.*' ).match( name )
        if aux:
            itype = aux.group( 1 )
        self.itype = itype
        self.ats = list( self.arg_types )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """Subclasses either implement __call__ or plug in a proper _call()."""
        return self._call( diff_var, chunk_size, **kwargs )

    def _call( self, diff_var = None, chunk_size = None, **kwargs ):
        output( 'base class method called for %s' % self.__class__.__name__ )
        raise RuntimeError
    
    ##
    # c: 21.03.2008, r: 21.03.2008
    def get_shape( self, diff_var, apr, apc = None ):
        output( 'base class method called for %s' % self.__class__.__name__ )
        raise RuntimeError

    ##
    # c: 21.03.2008, r: 21.03.2008
    def build_c_fun_args( self, *args, **kwargs ):
        output( 'base class method called for %s' % self.__class__.__name__ )
        raise RuntimeError

    ##
    # 16.11.2005, c
    def get_arg_names( self ):
        return self.__arg_names

    def set_arg_names( self, val ):
        if len( val ) != len( self.__class__.arg_types ):
            raise ValueError, 'equal shapes: %s, %s' \
                  % (val, self.__class__.arg_types)
        self.__arg_names = val
    arg_names = property( get_arg_names, set_arg_names )

    ##
    # 24.07.2006, c
    def classify_args( self, variables ):
        """state variable can be in place of parameter variable and vice
        versa."""
        self.names = Struct( name = 'arg_names',
                             material = [], variable = [], user = [],
                             state = [], virtual = [], parameter = [],
                             material_split = [] )

        msg = "variable '%s' requested by term '%s' does not exist!"
        arg_types = []
        for ii, arg_type in enumerate( self.arg_types ):
            name = self.__arg_names[ii]
            for single_type in arg_type.split( '|' ):
                if _match_var( single_type ):
                    names = self.names.variable
                    try:
                        var = variables[name]
                    except KeyError:
                        raise KeyError( msg )
##                     print arg_types, single_type, name, var.flags

                    if _match_state( single_type ) and \
                           var.is_state_or_parameter():
                        self.names.state.append( name )
                    elif (single_type == 'virtual') and var.is_virtual():
                        self.names.virtual.append( name )
                    elif _match_parameter( single_type ) and \
                             var.is_state_or_parameter():
                        self.names.parameter.append( name )
                    else:
                        continue
                    arg_types.append( single_type )

                elif _match_material( arg_type ):
                    names = self.names.material
                    match = _match_material_root( name )
                    if match:
                        self.names.material_split.append( (match.group( 1 ),
                                                          match.group( 2 )) )
                    else:
                        self.names.material_split.append( (name, None) )
                    arg_types.append( single_type )
                else:
                    names = self.names.user
                    arg_types.append( single_type )
                break
            names.append( name )

        self.n_virtual = len( self.names.virtual )
        if self.n_virtual > 1:
            raise ValueError( 'at most one virtial variable is allowed! (%d)'\
                              % self.n_virtual )

        self.ats = arg_types

        self.set_arg_types()

    def set_arg_types( self ):
        pass
        
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
        name = self.names.virtual[0]
        if variables is None:
            return name
        else:
            if variables[name].kind == 'test':
                return name
            else:
                output( 'variable %s is not virtual!' % name )
                raise ValueError

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

    ##
    # c: 24.07.2006, r: 15.01.2008
    def get_args( self, arg_types = None, **kwargs ):
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
                args.append( mat.get_data( region_name, ig, split[1] ) )
            else:
                args.append( kwargs[name] )
        return args

    ##
    # 24.07.2006, c
    def get_arg_name( self, arg_type ):
        ii = self.ats.index( arg_type )
        return self.arg_names[ii]

    ##
    # c: 29.11.2007, r: 10.04.2008
    def describe_geometry( self, geometries, variables, integrals ):

        try:
            integral = integrals[self.integral_name]
        except ValueError:
            output( 'integral %s is not defined!' % self.integral_name )
            raise
            
        integral.create_qp()
        tgs = self.get_geometry()
        for var_name in self.get_variable_names():
##             print '>>>>>', self.name, var_name

            variable = variables[var_name]
            field = variable.field
            assert field.region.contains( self.region )
            
##             print field.name, field.region_name
##             print field.bases

            if tgs.has_key( var_name ):
##                 print tgs[var_name]
                field.aps.describe_geometry( field, geometries,
                                             tgs[var_name], integral )

##             print field.aps.aps_per_group
##             pause()

    ##
    # 23.08.2006, c
    # 24.08.2006
    def get_geometry( self ):
        geom = self.__class__.geometry
        if geom:
            out = {}
            for (gtype, arg_type) in geom:
                arg_name = self.get_arg_name( arg_type )
                out[arg_name] = Struct( gtype = gtype,
                                       region = self.region )
            return out
        else:
            return None

    ##
    # c: 28.08.2006, r: 15.01.2008
    def get_current_group( self ):
        return (self.integral_name, self.region.name, self.char_fun.ig)

    ##
    # 11.10.2006, c
    def get_dof_conn_type( self ):
        return self.dof_conn_type, self.region.name

    ##
    # c: 16.02.2007, r: 15.01.2008
    def set_current_group( self, ig ):
        self.char_fun.set_current_group( ig )

    ##
    # c: 27.02.2007, r: 15.04.2008
    def get_cache( self, base_name, ii ):
        args = self.use_caches[base_name][ii]
        ans = [self.get_arg_name( arg ) for arg in args if not type( arg ) == dict]
##         print args, ans
##         pause()
        cname = '_'.join( [base_name] + ans )
        return self.caches[cname]

    ##
    # 02.03.2007, c
    def igs( self ):
        return self.char_fun.igs

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
