from base import Struct
import os.path as op

##
# 16.06.2005, c
class Reader( Struct ):
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
