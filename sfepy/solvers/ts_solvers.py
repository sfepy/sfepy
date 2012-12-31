from sfepy.base.base import output, get_default, Struct
from sfepy.solvers.solvers import make_get_conf, TimeSteppingSolver
from sfepy.solvers.ts import TimeStepper

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

class ExplicitTimeSteppingSolver(SimpleTimeSteppingSolver):
    name = 'ts.explicit'

    @staticmethod
    def process_conf(conf, kwargs):
        """
        Process configuration options.
        """
        get = make_get_conf(conf, kwargs)
        common = SimpleTimeSteppingSolver.process_conf(conf, kwargs)

        return Struct(mass=get('mass', None,
                               'missing "mass" in options!'),
                      lumped=get('lumped', False)) + common
