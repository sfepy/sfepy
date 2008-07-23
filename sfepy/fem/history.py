from sfepy.base.base import *
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.meshio import HDF5MeshIO

##
# 14.06.2007, c
class Histories( Container ):

    def from_file_hdf5( file_name, var_names ):
        """TODO: do not read entire file, provide data on demand."""
        io = HDF5MeshIO( file_name )
        ts = TimeStepper( *io.read_time_stepper() )
        ths = io.read_variables_time_history( var_names, ts )

        steps = nm.arange( ts.n_step, dtype = nm.int32 )
        objs = OneTypeList( History )
        for name, th in ths.iteritems():
            hist = History( name,
                            steps = steps,
                            th = th )
            objs.append( hist )
            
        obj = Histories( objs,
                         name = ' '.join( var_names ) )
        return obj
    from_file_hdf5 = staticmethod( from_file_hdf5 )


##
# 14.06.2007, c
class History( Struct ):

    ##
    # 14.06.2007, c
    def from_sequence( seq, name ):
        obj = History( name = name,
                       steps = nm.arange( len( seq ), dtype = nm.int32 ),
                       th = list( seq ) )
        return obj
    from_sequence = staticmethod( from_sequence )

    def __init__( self, name, th = None, steps = None, times = None ):
        Struct.__init__( self, name = name,
                         th = get_default( th, [] ),
                         steps = get_default( steps, [] ),
                         times = get_default( times, [] ) )

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

    def append( self, item, step, time ):
        self.th.append( item )
        self.steps.append( step )
        self.times.append( time )
