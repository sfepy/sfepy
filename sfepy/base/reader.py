from base import Struct
import os.path as op

##
# 16.06.2005, c
class Reader( Struct ):
    """
    Reads and executes a Python file as a script with execfile(), storing its
    locals. Then sets the __dict__ of a new instance of objClass to the stored
    locals.

    Example:

    >>> class A:
    >>>    pass

    >>> read = Reader( '.' )
    >>> instanceOfA = read( A, 'file.py' )

    It is equivalent to:

    >>> mod = __import__( 'file' )
    >>> instanceOfA = A()
    >>> instanceOfA.__dict__.update( mod.__dict__ )

    The first way does not create the 'file.pyc'...
    """
    ##
    # 16.06.2005, c
    def __init__( self, directory ):
        self.directory = directory

    ##
    # 16.06.2005, c
    # 17.10.2005
    # 09.02.2006
    def __call__( self, objClass, name ):
        fileName = op.join( self.directory, name + '.py' )

        aux = {}
        execfile( fileName, {}, aux )

        obj = objClass()
        for key, val in aux.iteritems():
            obj.__dict__[key] = val

        return obj
