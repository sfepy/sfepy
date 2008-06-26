from sfepy.base.base import *
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.meshio import HDF5MeshIO

##
# 14.06.2007, c
class Histories( Container ):

    ##
    # c: 14.06.2007, r: 23.06.2008
    def fromFileHDF5( fileName, varNames ):
        """TODO: do not read entire file, provide data on demand."""
        io = HDF5MeshIO( fileName )
        ts = TimeStepper( *io.readTimeStepper() )
        ths = io.readVariablesTimeHistory( varNames, ts )

        steps = nm.arange( ts.nStep, dtype = nm.int32 )
        objs = OneTypeList( History )
        for name, th in ths.iteritems():
            hist = History( name = name,
                            steps = steps,
                            th = th )
            objs.append( hist )
            
        obj = Histories( objs,
                         name = ' '.join( varNames ) )
        return obj
    fromFileHDF5 = staticmethod( fromFileHDF5 )


##
# 14.06.2007, c
class History( Struct ):

    ##
    # 14.06.2007, c
    def fromSequence( seq, name ):
        obj = History( name = name,
                       steps = nm.arange( len( seq ), dtype = nm.int32 ),
                       th = list( seq ) )
        return obj
    fromSequence = staticmethod( fromSequence )

##     def __setitem__( self, key, value ):
##         pass

    ##
    # 14.06.2007, c
    def __getitem__( self, ii ):
        return self.th[ii]

    ##
    # 15.06.2007, c
    def __len__( self ):
        return len( self.th )
