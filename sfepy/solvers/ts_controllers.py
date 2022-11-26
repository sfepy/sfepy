import numpy as nm

from sfepy.base.base import Struct
from sfepy.solvers.solvers import TimeStepController

class FixedTCS(TimeStepController):
    """
    Fixed (do-nothing) time step controller.
    """
    name = 'tsc.fixed'

class TimesSequenceTCS(TimeStepController):
    """
    Given times sequence time step controller.
    """
    name = 'tsc.time_sequence'

    _parameters = [
        ('times', 'iterable', range(1, 6), True,
         'A sequence of times to generate.'),
    ]

    def __init__(self, conf, **kwargs):
        TimeStepController.__init__(self, conf=conf, **kwargs)

        self.iter_times = iter(self.conf.times)

    def get_initial_dt(self, ts, vec, **kwargs):
        # This cannot be called repeatedly!
        self.t0 = ts.t0
        return next(self.iter_times) - self.t0

    def __call__(self, ts, vec0, vec1, **kwargs):
        self.t0 = ts.time
        try:
            t1 = next(self.iter_times)

        except StopIteration:
            t1 = ts.t1

        new_dt = t1 - self.t0
        if new_dt == 0.0:
            new_dt = 1.0

        status = Struct(u_err=None, v_err=None, emax=None, result='accept')

        return new_dt, status

class ElastodynamicsBasicTCS(TimeStepController):
    """
    Adaptive time step controller for elastodynamics by [1].

    [1] Grafenhorst, Matthias, Joachim Rang, and Stefan Hartmann.
    “Time-Adaptive Finite Element Simulations of Dynamical Problems for
    Temperature-Dependent Materials.” Journal of Mechanics of Materials and
    Structures 12, no. 1 (November 26, 2016): 57–91.
    https://doi.org/10.2140/jomms.2017.12.57.
    """
    name = 'tsc.ed_basic'

    _parameters = [
        ('eps_r', 'list of floats or float', None, True,
         'Relative tolerance(s).'),
        ('eps_a', 'list of floats or float', None, True,
         'Absolute tolerance(s).'),
        ('fmin', 'float', 0.3, False,
         'Minimum step size change on step rejection.'),
        ('fmax', 'float', 2.5, False,
         'Maximum step size change on step acceptance.'),
        ('fsafety', 'float', 0.8, False,
         'Step size change safety factor.'),
    ]

    def __call__(self, ts, vec0, vec1, unpack, **kwargs):
        conf = self.conf
        dt = ts.dt
        u_eps_a, u_eps_r, v_eps_a, v_eps_r, fmin, fmax, fsafety = (
            conf.eps_a[0], conf.eps_r[0],
            conf.eps_a[1], conf.eps_r[1],
            conf.fmin, conf.fmax, conf.fsafety,
        )
        u0, v0, a0 = unpack(vec0)
        u1, v1, a1 = unpack(vec1)

        # Backward Euler step.
        u1_be = u0 + dt * v1
        v1_be = v0 + dt * a1
        # Truncation error estimates.
        u_terr = u1 - u1_be
        v_terr = v1 - v1_be

        def get_err(terr, eps_a, eps_r):
            return nm.sqrt(
                1 / len(terr) *
                nm.sum((terr / (eps_r * nm.abs(terr) + eps_a))**2)
            )
        u_err = get_err(u_terr, u_eps_a, u_eps_r)
        v_err = get_err(v_terr, v_eps_a, v_eps_r)
        emax = max(u_err, v_err)
        status = Struct(u_err=u_err, v_err=v_err, emax=emax)
        if emax <= 1:
            # Step accepted.
            new_dt = dt * min(fmax, fsafety / nm.sqrt(emax))
            status.result = 'accept'

        else:
            new_dt = dt * max(fmin, fsafety / nm.sqrt(emax))
            status.result = 'reject'

        return new_dt, status
