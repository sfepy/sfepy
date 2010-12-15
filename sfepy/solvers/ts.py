from sfepy.base.base import *
from sfepy.solvers.solvers import TimeSteppingSolver


def get_print_info( n_step ):
    if n_step > 1:
        n_digit = int( nm.log10( n_step - 1 ) + 1 )

    else:
        n_digit = 1

    format = '%%%dd of %%%dd' % (n_digit, n_digit)
    suffix = '%%0%dd' % n_digit

    return n_digit, format, suffix

class TimeStepper( Struct ):
    
    def from_conf( conf ):
        return TimeStepper( conf.t0, conf.t1, conf.dt, conf.n_step )
    from_conf = staticmethod( from_conf )

    def __init__(self, t0, t1, dt, n_step, step=None, is_quasistatic=False):
        self.set_from_data(t0, t1, dt, n_step, step=step)
        self.is_quasistatic = is_quasistatic

    def set_from_data(self, t0, t1, dt, n_step, step=None):
        self.t0, self.t1, self.dt, self.n_step = t0, t1, dt, int(n_step)

        if not hasattr( self, 'n_step' ):
            self.n_step = round( nm.floor( ((self.t1 - self.t0) / self.dt)
                                          + 0.5 ) + 1.0 );

        if self.n_step > 1:
            self.times, self.dt = nm.linspace( self.t0, self.t1, self.n_step,
                                               endpoint = True, retstep = True )
        else:
            self.times = nm.array( (self.t0,), dtype = nm.float64 )
            self.dt = self.t1 - self.t0

        self.n_digit, self.format, self.suffix = get_print_info( self.n_step )

        self.set_step(step)

    def set_from_ts(self, ts, step=None):
        step = get_default(step, ts.step)
        self.set_from_data(ts.t0, ts.t1, ts.dt, ts.n_step, step=step)

    def __iter__( self ):
        """ts.step, ts.time is consistent with step, time returned here
        ts.nt is normalized time in [0, 1]"""
        return self.iter_from( 0 )

    def iter_from( self, step ):
        self.step = step - 1

        for time in self.times[step:]:

            self.time = time
            self.step += 1
            self.normalize_time()
    
            yield self.step, self.time

    def normalize_time( self ):
        if self.n_step > 1:
            self.nt = float( self.step ) / (self.n_step - 1)
        else:
            self.nt = 0.0
        
    def set_step( self, step = -1, nt = 1.0 ):
        nm1 = self.n_step - 1
        if step is None:
            step = int( round( nt * nm1 ) )
#            print step
        if step < 0:
            step = self.n_step + step
        if (step >= self.n_step) or (step < 0):
            output( 'time step must be in [%d, %d]' % (-nm1, nm1)  )
            raise ValueError

        self.step = step
        self.time = self.times[step]
        self.normalize_time()

    def __eq__( self, other ):

        if type( other ) == type( self ):
            return (abs( self.t0 == other.t0 ) < 1e-15) and \
                   (abs( self.t1 == other.t1 ) < 1e-15) and \
                   (self.n_step == other.n_step)
        else:
            raise ValueError

class SimpleTimeSteppingSolver( TimeSteppingSolver ):
    name = 'ts.simple'

    def __init__(self, conf, ts=None, **kwargs):
        TimeSteppingSolver.__init__( self, conf, **kwargs )

        self.ts = get_default(ts, TimeStepper.from_conf(conf))
        self.ts.is_quasistatic = conf.get_default_attr('quasistatic', False)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format

    def __call__( self, state0 = None, conf = None,
                  step_fun = None, step_args = None ):

        step_fun = get_default( step_fun, self.step_fun )
        step_args = get_default( step_args, self.step_args )

        for step, time in self.ts:
            output( self.format % (time, step + 1, self.ts.n_step) )

            state = step_fun( self.ts, state0, *step_args )
            state0 = state.copy(deep=True)
            yield self.ts, state
