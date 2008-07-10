from sfepy.base.base import *
from sfepy.solvers.solvers import TimeSteppingSolver


##
# c: 09.07.2008, r: 10.07.2008
def getPrintInfo( nStep ):
    nDigit = int( nm.log10( nStep - 1 ) + 1 )
    format = '%%%dd of %%%dd' % (nDigit, nDigit)
    suffix = '.%%0%dd.vtk' % nDigit
    return nDigit, format, suffix

##
# 17.07.2006, c
class TimeStepper( Struct ):
    
    ##
    # c: 17.07.2006, r: 06.02.2008
    def fromConf( conf ):
        return TimeStepper( conf.t0, conf.t1, conf.dt, conf.nStep )
    fromConf = staticmethod( fromConf )

    ##
    # c: 19.09.2006, r: 10.07.2008
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

        self.nDigit, self.format, self.suffix = getPrintInfo( self.nStep )
        
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
    # c: 19.09.2006, r: 13.06.2008
    def setStep( self, step = -1, nt = 1.0 ):
        nm1 = self.nStep - 1
        if step is None:
            step = int( round( nt * nm1 ) )
#            print step
        if step < 0:
            step = self.nStep + step
        if (step >= self.nStep) or (step < 0):
            output( 'time step must be in [%d, %d]' % (-nm1, nm1)  )
            raise ValueError

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

##
# c: 06.02.2008, r: 06.02.2008
class SimpleTimeSteppingSolver( TimeSteppingSolver ):
    name = 'ts.simple'

    ##
    # c: 06.02.2008, r: 06.02.2008
    def __init__( self, conf, **kwargs ):
        TimeSteppingSolver.__init__( self, conf, **kwargs )

        self.ts = TimeStepper.fromConf( conf )
        nd = self.ts.nDigit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format

    ##
    # c: 06.02.2008, r: 06.02.2008
    def __call__( self, state0 = None, conf = None,
                  stepFun = None, stepArgs = None ):

        stepFun = getDefault( stepFun, self.stepFun )
        stepArgs = getDefault( stepArgs, self.stepArgs )

        for step, time in self.ts:
            output( self.format % (time, step + 1, self.ts.nStep) )

            state = stepFun( self.ts, state0, *stepArgs )
            state0 = state.copy()
            yield step, time, state
