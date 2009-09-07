##
# 16.02.2005, c
import numpy as nm
import scipy as sc
import scipy.linalg as nla
import scipy.sparse as sp
import pdb

import glob, re, time, sys, os
from copy import copy, deepcopy
from types import MethodType, UnboundMethodType
from getch import getch

import atexit

real_types = [nm.float64]
complex_types = [nm.complex128]

nm.set_printoptions( threshold = 100 )

sfepy_config_dir = os.path.expanduser('~/.sfepy')
if not os.path.exists(sfepy_config_dir):
    os.makedirs(sfepy_config_dir)

##
# 22.09.2005, c
# 24.10.2005
if sys.version[:5] < '2.4.0':
    def sorted( sequence ):
        tmp = copy( sequence )
        tmp.sort()
        return tmp

##
# 26.05.2006, c
def debug():
    pdb.set_trace()

##
# c: 20.06.2008, r: 20.06.2008
def import_file( filename ):
    path = os.path.dirname( filename )
    if not path in sys.path:
        sys.path.append( path )
    name = os.path.splitext( os.path.basename( filename ) )[0]
    mod = __import__( name )
    return mod


def assert_( condition ):
    if not condition:
        raise ValueError( 'assertion failed!' )

##
# c: 06.04.2005, r: 05.05.2008
def pause( msg = None ):
    """
    Prints the line number and waits for a keypress.

    If you press:
    "q" ............. it will call sys.exit()
    any other key ... it will continue execution of the program

    This is useful for debugging.
    """
    f = sys._getframe(1)
    ff = f.f_code
    if (msg):
        print '%s, %d: %s(), %d: %s' % (ff.co_filename, ff.co_firstlineno,
                                        ff.co_name, f.f_lineno, msg)
    else:
        print '%s, %d: %s(), %d' % (ff.co_filename, ff.co_firstlineno,
                                    ff.co_name, f.f_lineno)
    spause()

##
# Silent pause.
# 18.02.2005, c
# 12.02.2007
def spause( msg = None ):
    """
    Waits for a keypress.

    If you press:
    "q" ............. it will call sys.exit()
    any other key ... it will continue execution of the program

    This is useful for debugging. This function is called from pause().
    """
    if (msg):
        print msg
    sys.stdout.flush()
    ch = getch()
    if ch == 'q':
        sys.exit()

##
# 02.01.2005
class Struct( object ):
    # 03.10.2005, c
    # 26.10.2005
    def __init__( self, **kwargs ):
        if kwargs:
            self.__dict__.update( kwargs )
        
    # 08.03.2005
    def __str__( self ):
        """Print instance class, name and items in alphabetical order."""
        ss = "%s" % self.__class__.__name__
        if hasattr( self, 'name' ):
            ss += ":%s" % self.name
        ss += '\n'

        keys, vals = self.__dict__.keys(), self.__dict__.values()
        order = nm.argsort(keys)
        for ii in order:
            key, val = keys[ii], vals[ii]

            if issubclass( val.__class__, Struct ):
                ss += "  %s:\n    %s" % (key, val.__class__.__name__)
                if hasattr( val, 'name' ):
                    ss += ":%s" % val.name
                ss += '\n'
            else:
                aux = "\n" + str( val )
                aux = aux.replace( "\n", "\n    " );
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return( ss.rstrip() )

    def __repr__( self ):
        ss = "%s" % self.__class__.__name__
        if hasattr( self, 'name' ):
            ss += ":%s" % self.name
        return ss

    ##
    # 28.08.2007, c
    def __add__( self, other ):
        """Merge Structs. Attributes of new are those of self unless an
        attribute and its counterpart in other are both Structs - these are
        merged then."""
        new = copy( self )
        for key, val in other.__dict__.iteritems():
            if hasattr( new, key ):
                sval = getattr( self, key )
                if issubclass( sval.__class__, Struct ) and \
                       issubclass( val.__class__, Struct ):
                    setattr( new, key, sval + val )
                else:
                    setattr( new, key, sval )
            else:
                setattr( new, key, val )
        return new

    ##
    # 28.08.2007, c
    def __iadd__( self, other ):
        """Merge Structs in place. Attributes of self are left unchanged
        unless an attribute and its counterpart in other are both Structs -
        these are merged then."""
        for key, val in other.__dict__.iteritems():
            if hasattr( self, key ):
                sval = getattr( self, key )
                if issubclass( sval.__class__, Struct ) and \
                       issubclass( val.__class__, Struct ):
                    setattr( self, key, sval + val )
            else:
                setattr( self, key, val )
        return self

    # 08.03.2005, c
    def str_all( self ):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.iteritems():
            if issubclass( self.__dict__[key].__class__, Struct ):
                ss += "  %s:\n" % key
                aux = "\n" + self.__dict__[key].str_all()
                aux = aux.replace( "\n", "\n    " );
                ss += aux[1:] + "\n"
            else:
                aux = "\n" + str( val )
                aux = aux.replace( "\n", "\n    " );
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return( ss.rstrip() )

    ##
    # 09.07.2007, c
    def to_dict( self ):
        return copy( self.__dict__ )

    def get_default_attr( self, key, default = None, msg_if_none = None ):
        """Behaves like dict.get() if msg_if_none is None."""
        return get_default_attr( self, key, default, msg_if_none )

    def set_default_attr( self, key, default = None ):
        """Behaves like dict.setdefault()."""
        return self.__dict__.setdefault( key, default )

    def copy(self, deep=False, name=None):
        """Make a (deep) copy of self.

        Parameters
        ----------
        deep : bool
            Make a deep copy.
        name : str
            Name of the copy, with default self.name + '_copy'.
        """
        if deep:
            other = deepcopy(self)
        else:
            other = copy(self)

        if hasattr(self, 'name'):
            other.name = get_default(name, self.name + '_copy')

        return other
#
# 12.07.2007, c
class IndexedStruct( Struct ):

    ##
    # 12.07.2007, c
    def __getitem__( self, key ):
        return getattr( self, key )

    ##
    # 12.07.2007, c
    def __setitem__( self, key, val ):
        setattr( self, key, val )

##
# 14.07.2006, c
class Container( Struct ):

    def __init__( self, objs = None, **kwargs ):
        Struct.__init__( self, **kwargs )

        if objs is not None:
            self._objs = objs
            self.update()
        else:
            self._objs = []
            self.names = []

    def update( self, objs = None ):
        if objs is not None:
            self._objs = objs

        self.names = [obj.name for obj in self._objs]

    ##
    # 14.11.2005, c
    def __getitem__( self, ii ):
        if isinstance( ii, int ):
            return self._objs[ii]
        elif isinstance( ii, str ):
            return self._objs[self.names.index( ii )]
        else:
            raise IndexError

    def __iter__( self ):
        return self._objs.__iter__()
    ##
    # 18.07.2006, c
    def __len__( self ):
        return len( self._objs )

    def append( self, obj ):
        self._objs.append( obj )
        self.names.append( obj.name )

    def remove_name( self, name ):
        ii = self.names.index[name]
        del self.names[ii]
        del self._objs[ii]

    ##
    # dict-like methods.
    def itervalues( self ):
        return self._objs.__iter__()

    def iterkeys( self ):
        return self.get_names().__iter__()

    def iteritems( self ):
        for obj in self._objs:
            yield obj.name, obj

    ##
    # 20.09.2006, c
    def has_key( self, ii ):
        if isinstance( ii, int ):
            if (ii < len( self )) and (ii >= (-len( self ))):
                return True
            else:
                return False
        elif isinstance( ii, str ):
            try:
                self.names.index( ii )
                return True
            except:
                return False
        else:
            raise IndexError, 'unsupported index type: %s' % key
        
    ##
    # 12.06.2007, c
    def print_names( self ):
        print [obj.name for obj in self._objs]

    def get_names( self ):
        return [obj.name for obj in self._objs]

        
##
# 30.11.2004, c
# 01.12.2004
# 01.12.2004
class OneTypeList( list ):
    def __init__( self, item_class ):
        self.item_class = item_class
        pass
    
    def __setitem__( self, key, value ):
        if (type( value ) in (list, tuple)):
            for ii, val in enumerate( value ):
                if (val.__class__ != self.item_class):
                    raise TypeError
        else:
            if (value.__class__ != self.item_class):
                raise TypeError
        list.__setitem__( self, key, value )

    ##
    # 21.11.2005, c
    def __getitem__( self, ii ):
        if isinstance( ii, int ):
            return list.__getitem__( self, ii )
        elif isinstance( ii, str ):
            ir = self.find( ii, ret_indx = True )
            if ir:
                return list.__getitem__( self, ir[0] )
            else:
                raise IndexError, ii
        else:
            raise IndexError, ii
    

    def __str__( self ):
        ss = "[\n"
        for ii in self:
            aux = "\n" + ii.__str__()
            aux = aux.replace( "\n", "\n  " );
            ss += aux[1:] + "\n"
        ss += "]"
        return( ss )
    
    def find( self, name, ret_indx = False ):
        for ii, item in enumerate( self ):
            if item.name == name:
                if ret_indx:
                    return ii, item
                else:
                    return item
        return None

    ##
    # 12.06.2007, c
    def print_names( self ):
        print [ii.name for ii in self]

    def get_names( self ):
        return [ii.name for ii in self]

class Output( Struct ):
    """Factory class providing output (print) functions.

    Example:

    >>> output = Output( 'sfepy:' )
    >>> output( 1, 2, 3, 'hello' )
    >>> output.prefix = 'my_cool_app:'
    >>> output( 1, 2, 3, 'hello' )
    """

    def __init__(self, prefix, filename=None, combined=False, **kwargs):
        Struct.__init__(self, **kwargs)

        self.prefix = prefix

        self.set_output(filename, combined)
        
    def __call__(self, *argc, **argv):
        self.output_function(*argc, **argv)

    def set_output(self, filename=None, combined=False, append=False):
        """Set the output function - all SfePy printing is accomplished by
        it. If filename is None, output is to screen only, otherwise it is to
        the specified file, moreover, if combined is True, both the ways are
        used.

        Arguments:
                filename - print into this file
                combined - print both on screen and into a file
                append - append to an existing file instead of overwriting it
        """
        self.level = 0
        def output_screen( *argc, **argv ):
            format = '%s' + ' %s' * (len( argc ) - 1)
            msg =  format % argc

            if msg.startswith( '...' ):
                self.level -= 1

            print self._prefix + ('  ' * self.level) + msg

            if msg.endswith( '...' ):
                self.level += 1

        def output_file( *argc, **argv ):
            format = '%s' + ' %s' * (len( argc ) - 1)
            msg =  format % argc

            if msg.startswith( '...' ):
                self.level -= 1

            fd = open( filename, 'a' )
            print >>fd, self._prefix + ('  ' * self.level) + msg
            fd.close()

            if msg.endswith( '...' ):
                self.level += 1

        def output_combined( *argc, **argv ):
            output_screen( *argc, **argv )
            output_file( *argc, **argv )
    
        if filename is None:
            self.output_function = output_screen

        else:
            if not append:
                fd = open( filename, 'w' )
                fd.close()

            if combined:
                self.output_function = output_combined
            else:
                self.output_function = output_file

    def get_output_function(self):
        return self.output_function

    def set_output_prefix( self, prefix ):
        assert_( isinstance( prefix, str ) )
        if len( prefix ) > 0:
            prefix += ' '
        self._prefix = prefix
        
    def get_output_prefix( self ):
        return self._prefix[:-1]
    prefix = property( get_output_prefix, set_output_prefix )
    
output = Output( 'sfepy:' )

def print_structs(objs):
    """Print Struct instances in a container, works recursively. Debugging
    utility function."""
    if isinstance(objs, dict):
        for key, vals in objs.iteritems():
            print key
            print_structs(vals)
    elif isinstance(objs, list):
        for vals in objs:
            print_structs(vals)
    else:
        print objs

def iter_dict_of_lists(dol, return_keys=False):
    for key, vals in dol.iteritems():
        for ii, val in enumerate(vals):
            if return_keys:
                yield key, ii, val
            else:
                yield val

##
# 19.07.2005, c
# 26.05.2006
# 17.10.2007
def dict_to_struct( *args, **kwargs ):

    try:
        level = kwargs['level']
    except:
        level = 0
        
    try:
        flag = kwargs['flag']
    except:
        print 'mask undefined!'
        raise

    # For level 0 only...
    try:
        constructor = kwargs['constructor']
    except:
        constructor = Struct

    out = []
    for arg in args:
        if type( arg ) == dict:
            if flag[level]:
                aux = constructor()
            else:
                aux = {}
                
            for key, val in arg.iteritems():
                if type( val ) == dict:
                    try:
                        flag[level+1]
                    except:
                        flag = flag + (0,)
                    val2 = dict_to_struct( val, level = level + 1, flag = flag )
                    if flag[level]:
                        aux.__dict__[key] = val2
                    else:
                        aux[key] = val2
                else:
                    if flag[level]:
                        aux.__dict__[key] = val
                    else:
                        aux[key] = val
            out.append( aux )
        else:
            out.append( arg )

    if len( out ) == 1:
        out = out[0]

    return out

##
# 23.01.2006, c
def is_sequence( var ):
    if issubclass( var.__class__, tuple ) or issubclass( var.__class__, list ):
        return True
    else:
        return False

##
# 17.10.2007, c
def is_derived_class( cls, parent ):
    return issubclass( cls, parent ) and (cls is not parent)

##
# 23.10.2007, c
def insert_static_method( cls, function ):
    setattr( cls, function.__name__, staticmethod( function ) )

##
# 23.10.2007, c
def insert_method( instance, function ):
    setattr( instance, function.__name__,
             UnboundMethodType( function, instance, instance.__class__ ) )

def use_method_with_name( instance, method, new_name ):
    setattr( instance, new_name, method )

def insert_as_static_method( cls, name, function ):
    setattr( cls, name, staticmethod( function ) )

##
# 09.08.2006, c
def invert_dict( d, is_val_tuple = False ):
    di = {}
    for key, val in d.iteritems():
        if is_val_tuple:
            for v in val:
                di[v] = key
        else:
            di[val] = key
    return di

##
# 24.08.2006, c
# 05.09.2006
def dict_from_keys_init( keys, seq_class = None ):

    if seq_class is None:
        return {}.fromkeys( keys )
    
    out = {}
    for key in keys:
        out[key] = seq_class()
    return out

##
# 16.10.2006, c
def dict_extend( d1, d2 ):
    for key, val in d1.iteritems():
        val.extend( d2[key] )

def set_defaults( dict_, defaults ):
    for key, val in defaults.iteritems():
        dict_.setdefault( key, val )

##
# c: 12.03.2007, r: 04.04.2008
def get_default( arg, default, msg_if_none = None ):
    
    if arg is None:
        out = default
    else:
        out = arg

    if (out is None) and (msg_if_none is not None):
        raise ValueError( msg_if_none )

    return out

##
# c: 28.04.2008, r: 28.04.2008
def get_default_attr( obj, attr, default, msg_if_none = None ):
    if hasattr( obj, attr ):
        out = getattr( obj, attr )
    else:
        out = default

    if (out is None) and (msg_if_none is not None):
        raise ValueError( msg_if_none )

    return out

##
# c: 27.02.2008, r: 27.02.2008
def select_by_names( objs_all, names, replace = None, simple = True ):
    objs = {}
    for key, val in objs_all.iteritems():
        if val.name in names:
            if replace is None:
                objs[key] = val
            else:
                new_val = copy( val )
                old_attr = getattr( val, replace[0] )
                if simple:
                    new_attr = old_attr % replace[1]
                    setattr( new_val, replace[0], new_attr )
                else:
                    new_attr = replace[1].get( val.name, old_attr )
                    setattr( new_val, replace[0], new_attr )
                objs[key] = new_val
    return objs

def ordered_iteritems(adict):
    keys = adict.keys()
    order = nm.argsort(keys)
    for ii in order:
        key = keys[ii]
        yield key, adict[key]
