##
# 16.02.2005, c
import numpy as nm
import scipy as sc
import scipy.linalg as nla
import scipy.sparse as sp
import pdb

import glob, re, time, sys, os
from copy import *
from sets import Set
from types import UnboundMethodType

from getch import getch

import atexit

nm.set_printoptions( threshold = 100 )

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
# c: 02.04.2008, r: 02.04.2008
prefix = 'sfe:'
def setOutputPrefix( prefix  = 'sfe:' ):
    globals()['prefix'] = prefix

##
# c: 05.06.2006, r: 02.04.2008
level = 0
def output( *argc, **argv ):
    global level
    format = ' %s' * len( argc )
    msg =  format % argc

    if msg.startswith( ' ...' ):
        level -= 1

    print prefix + ('  ' * level) + msg

    if msg.endswith( '...' ):
        level += 1

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
        ss = "%s" % self.__class__.__name__
        if hasattr( self, 'name' ):
            ss += ":%s" % self.name
        ss += '\n'
        for key, val in self.__dict__.iteritems():
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
    def strAll( self ):
        ss = "%s\n" % self.__class__
        for key, val in self.__dict__.iteritems():
            if issubclass( self.__dict__[key].__class__, Struct ):
                ss += "  %s:\n" % key
                aux = "\n" + self.__dict__[key].strAll()
                aux = aux.replace( "\n", "\n    " );
                ss += aux[1:] + "\n"
            else:
                aux = "\n" + str( val )
                aux = aux.replace( "\n", "\n    " );
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return( ss.rstrip() )

    ##
    # 09.07.2007, c
    def toDict( self ):
        return copy( self.__dict__ )

##
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

        self._objs = objs
        if objs is not None:
            self.update()
        else:
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
    def printNames( self ):
        print [obj.name for obj in self._objs]

##
# 30.11.2004, c
# 01.12.2004
# 01.12.2004
class OneTypeList( list ):
    def __init__( self, itemClass ):
        self.itemClass = itemClass
        pass
    
    def __setitem__( self, key, value ):
        if (type( value ) in (list, tuple)):
            for ii, val in enumerate( value ):
                if (val.__class__ != self.itemClass):
                    raise TypeError
        else:
            if (value.__class__ != self.itemClass):
                raise TypeError
        list.__setitem__( self, key, value )

    ##
    # 21.11.2005, c
    def __getitem__( self, ii ):
        if isinstance( ii, int ):
            return list.__getitem__( self, ii )
        elif isinstance( ii, str ):
            ir = self.find( ii, retIndx = True )
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
    
    def find( self, name, retIndx = False ):
        for ii, item in enumerate( self ):
            if item.name == name:
                if retIndx:
                    return ii, item
                else:
                    return item
        return None

    ##
    # 12.06.2007, c
    def printNames( self ):
        print [ii.name for ii in self]

##
# 08.03.2005, c
# 06.10.2005
# 11.10.2005
# 17.10.2005
## def printMemStats():
##     import ccore
##     modules = [ii.__name__ for ii in vars( ccore ).values()
##                if (type( ii ) == type( ccore ))]
##     modules = sorted( [ii for ii in modules if '_' not in ii] )

##     fd = open( 'memory.log', 'w' )
##     for module in modules:
##         fd.write( '*** ' + module + ':\n' )
##         try:
##             sys.modules[module].mem_print( fd, 0 )
##         except:
##             pass

##     fd.close()
## atexit.register( printMemStats  )

##
# 19.07.2005, c
# 26.05.2006
# 17.10.2007
def dictToStruct( *args, **kwargs ):

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
                    val2 = dictToStruct( val, level = level + 1, flag = flag )
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
def isSequence( var ):
    if issubclass( var.__class__, tuple ) or issubclass( var.__class__, list ):
        return True
    else:
        return False

##
# 17.10.2007, c
def isDerivedClass( cls, parent ):
    return issubclass( cls, parent ) and (cls is not parent)

##
# 23.10.2007, c
def insertStaticMethod( cls, function ):
    setattr( cls, function.__name__, staticmethod( function ) )

##
# 23.10.2007, c
def insertMethod( instance, function ):
    setattr( instance, function.__name__,
             UnboundMethodType( function, instance, instance.__class__ ) )

##
# 09.08.2006, c
def invertDict( d, isValTuple = False ):
    di = {}
    for key, val in d.iteritems():
        if isValTuple:
            for v in val:
                di[v] = key
        else:
            di[val] = key
    return di

##
# 24.08.2006, c
# 05.09.2006
def dictFromKeysInit( keys, seqClass = None ):

    if seqClass is None:
        return {}.fromkeys( keys )
    
    out = {}
    for key in keys:
        out[key] = seqClass()
    return out

##
# 16.10.2006, c
def dictExtend( d1, d2 ):
    for key, val in d1.iteritems():
        val.extend( d2[key] )

##
# c: 12.03.2007, r: 04.04.2008
def getDefault( arg, default, msgIfNone = None ):
    
    if arg is None:
        out = default
    else:
        out = arg

    if (out is None) and (msgIfNone is not None):
        output( msgIfNone )
        raise ValueError

    return out

##
# c: 28.04.2008, r: 28.04.2008
def getDefaultAttr( obj, attr, default ):
    if hasattr( obj, attr ):
        out = getattr( obj, attr )
    else:
        out = default
    return out

##
# c: 27.02.2008, r: 27.02.2008
def selectByNames( objs_all, names, replace = None, simple = True ):
    objs = {}
    for key, val in objs_all.iteritems():
        if val.name in names:
            if replace is None:
                objs[key] = val
            else:
                newVal = copy( val )
                oldAttr = getattr( val, replace[0] )
                if simple:
                    newAttr = oldAttr % replace[1]
                    setattr( newVal, replace[0], newAttr )
                else:
                    newAttr = replace[1].get( val.name, oldAttr )
                    setattr( newVal, replace[0], newAttr )
                objs[key] = newVal
    return objs
