import numpy as nm

from sfepy.base.base import output, get_default, Struct
from sfepy.solvers.solvers import make_get_conf, TimeSteppingSolver


def get_print_info( n_step ):
    if n_step > 1:
        n_digit = int( nm.log10( n_step - 1 ) + 1 )

    else:
        n_digit = 1

    format = '%%%dd of %%%dd' % (n_digit, n_digit)
    suffix = '%%0%dd' % n_digit

    return n_digit, format, suffix

class TimeStepper( Struct ):
    """
    Time stepper class.
    """

    @staticmethod
    def from_conf( conf ):
        return TimeStepper(conf.t0, conf.t1, conf.dt, conf.n_step,
                           conf.quasistatic)

    def __init__(self, t0, t1, dt=None, n_step=None, step=None,
                 is_quasistatic=False):
        self.set_from_data(t0, t1, dt=dt, n_step=n_step, step=step)
        self.is_quasistatic = is_quasistatic

    def _get_n_step(self, t0, t1, dt):
        n_step = int(round(nm.floor(((t1 - t0) / dt) + 0.5) + 1.0))
        return n_step

    def set_from_data(self, t0, t1, dt=None, n_step=None, step=None):
        self.t0, self.t1 = t0, t1

        dt = get_default(dt, t1 - t0)
        self.n_step = get_default(n_step,
                                  self._get_n_step(self.t0, self.t1, dt))

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

    def normalize_time(self):
        self.nt = (self.time - self.t0) / (self.t1 - self.t0)

    def set_step(self, step=0, nt=0.0):
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

class VariableTimeStepper(TimeStepper):
    """
    Time stepper class with a variable time step.
    """

    def set_from_data(self, t0, t1, dt=None, n_step=None, step=None):
        self.t0, self.t1 = t0, t1

        self.dtime = self.t1 - self.t0
        dt = get_default(dt, self.dtime)

        self.n_step0 = get_default(n_step,
                                   self._get_n_step(self.t0, self.t1, dt))

        if self.n_step0 > 1:
            self.dt = self.dtime / (self.n_step0 - 1)

        else:
            self.dt = self.dtime

        self.dt0 = self.dt

        self.n_digit, self.format, self.suffix = get_print_info(5)

        self.set_step(step)

    def set_step(self, step=0, nt=0.0):
        self.dts = []
        self.step = 0
        self.n_step = 1

        if step is None:
            self.times = [self.t0 + self.dt0 * nt]
            self.nt = nt

        else:
            nt = float(step) / (self.n_step0 - 1)
            self.times = [self.t0 + self.dt0 * nt]
            self.nt = nt

    def get_default_time_step(self):
        return self.dt0

    def set_time_step(self, dt):
        self.dt = dt

    def __iter__(self):
        """
        ts.step, ts.time is consistent with step, time returned here
        ts.nt is normalized time in [0, 1]
        """
        self.set_step(0)

        while self.nt < 1.0:
            self.time = self.times[self.step]

            yield self.step, self.time

            self.step += 1
            self.time += self.dt
            self.normalize_time()

            self.times.append(self.time)
            self.dts.append(self.dt)
            self.n_step = self.step + 1

class SimpleTimeSteppingSolver( TimeSteppingSolver ):
    name = 'ts.simple'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Process configuration options.
        """
        get = make_get_conf(conf, kwargs)
        common = TimeSteppingSolver.process_conf(conf)

        return Struct(t0=get('t0', 0.0),
                      t1=get('t1', 1.0),
                      dt=get('dt', None),
                      n_step=get('n_step', 10),
                      quasistatic=get('quasistatic', False)) + common

    def __init__(self, conf, **kwargs):
        TimeSteppingSolver.__init__(self, conf, **kwargs)

        self.ts = TimeStepper.from_conf(self.conf)

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
