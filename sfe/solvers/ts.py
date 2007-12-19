from sfe.base.base import *

##
# 17.07.2006, c
class TimeStepper( Struct ):
    
    ##
    # 17.07.2006, c
    # 04.08.2006
    # 07.09.2006
    # 19.09.2006
    def fromConf( conf ):
        return TimeStepper( **conf )
    fromConf = staticmethod( fromConf )

    ##
    # 19.09.2006, c
    def __init__( self, t0, t1, dt, nStep ):
        self.t0, self.t1, self.dt, self.nStep = t0, t1, dt, int( nStep )

        self.time = None
        self.step = None
        self.nt = None

        if not hasattr( self, 'nStep' ):
            self.nStep = round( nm.floor( ((self.t1 - self.t0) / self.dt)
                                          + 0.5 ) + 1.0 );

        if self.nStep > 1:
            self.times, self.dt = nm.linspace( self.t0, self.t1, self.nStep,
                                               endpoint = True, retstep = True )
        else:
            self.times = nm.array( (self.t0,), dtype = nm.float64 )
        
    ##
    # 17.07.2006, c
    # 01.08.2006
    # 07.09.2006
    # 19.09.2006
    # 17.07.2007
    def __iter__( self ):
        """ts.step, ts.time is consistent with step, time returned here
        ts.nt is normalized time in [0, 1]"""
        return self.iterFrom( 0 )

    ##
    # 17.07.2007, c
    def iterFrom( self, step ):
        self.step = step - 1

        for time in self.times[step:]:

            self.time = time
            self.step += 1
            self.normalizeTime()
    
            yield self.step, self.time

    ##
    # 19.09.2006, c
    def normalizeTime( self ):
        if self.nStep > 1:
            self.nt = float( self.step ) / (self.nStep - 1)
        else:
            self.nt = 0.0
        
    ##
    # 19.09.2006, c
    # 29.09.2006
    def setStep( self, step = -1, nt = 1.0 ):
        if step is None:
            step = int( round( nt * (self.nStep - 1) ) )
#            print step
        if step < 0:
            step = self.nStep + step

        self.step = step
        self.time = self.times[step]
        self.normalizeTime()

    ##
    # 14.06.2007, c
    def __eq__( self, other ):

        if type( other ) == type( self ):
            return (abs( self.t0 == other.t0 ) < 1e-15) and \
                   (abs( self.t1 == other.t1 ) < 1e-15) and \
                   (self.nStep == other.nStep)
        else:
            raise ValueError
