from __future__ import absolute_import
from .base import Struct
import os.path as op
import six

##
# 16.06.2005, c
class Reader( Struct ):
    """
    Reads and executes a Python file as a script with execfile(), storing its
    locals. Then sets the __dict__ of a new instance of obj_class to the stored
    locals.

    Example:

    >>> class A:
    >>>    pass

    >>> read = Reader( '.' )
    >>> instance_of_a = read( A, 'file.py' )

    It is equivalent to:

    >>> mod = __import__( 'file' )
    >>> instance_of_a = A()
    >>> instance_of_a.__dict__.update( mod.__dict__ )

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
    def __call__( self, obj_class, name ):
        filename = op.join( self.directory, name + '.py' )

        aux = {}
        execfile( filename, {}, aux )

        obj = obj_class()
        for key, val in six.iteritems(aux):
            obj.__dict__[key] = val

        return obj
