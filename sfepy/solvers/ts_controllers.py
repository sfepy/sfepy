"""
Time step controllers.
"""
import numpy as nm

from sfepy.base.base import Struct
from sfepy.solvers.solvers import TimeStepController

class FixedTSC(TimeStepController):
    """
    Fixed (do-nothing) time step controller.
    """
    name = 'tsc.fixed'

class TimesSequenceTSC(TimeStepController):
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

def eval_scaled_norm(terr, eps_a, eps_r):
    return nm.sqrt(
        1 / len(terr) *
        nm.sum((terr / (eps_r * nm.abs(terr) + eps_a))**2)
    )

class ElastodynamicsBasicTSC(TimeStepController):
    """
    Adaptive time step I-controller for elastodynamics.

    The implementation is based on [1]. The default parameters correspond to
    the PID-Controller as implemented in ``tsc.ed_pid`` with P=D=0, I=1.

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
         'Minimum step size change factor on step rejection.'),
        ('fmax', 'float', 2.5, False,
         'Maximum step size change factor on step acceptance.'),
        ('fsafety', 'float', 0.8, False,
         'Step size change safety factor.'),
        ('error_order', 'float', 2, False,
         'The order of the solver error estimate.'),
        ('guess_dt0', 'bool', False, False,
         'Guess a good initial step size from initial conditions.'),
    ]

    @staticmethod
    def get_scaled_errors(dt, vec0, vec1, eps_as, eps_rs, unpack):
        u_eps_a, v_eps_a = eps_as
        u_eps_r, v_eps_r = eps_rs
        aux = unpack(vec0)
        u0, v0, a0 = aux if unpack.n_arg == 3 else (aux[0], aux[2], aux[3])
        aux = unpack(vec1)
        u1, v1, a1 = aux if unpack.n_arg == 3 else (aux[0], aux[2], aux[3])

        # Backward Euler step.
        u1_be = u0 + dt * v1
        v1_be = v0 + dt * a1
        # Truncation error estimates.
        u_terr = u1 - u1_be
        v_terr = v1 - v1_be

        u_err = eval_scaled_norm(u_terr, u_eps_a, u_eps_r)
        v_err = eval_scaled_norm(v_terr, v_eps_a, v_eps_r)

        return u_err, v_err

    def get_initial_dt(self, ts, vec, unpack, **kwargs):
        """
        Adapted from [1] for second order ODEs.

        [1] Hairer, Ernst, Gerhard Wanner, and Syvert P. Nørsett. Solving
        Ordinary Differential Equations I: Nonstiff Problems. Vol. 8. Springer
        Series in Computational Mathematics. Berlin, Heidelberg: Springer,
        1993. https://doi.org/10.1007/978-3-540-78862-1.
        """
        conf = self.conf
        if not conf.guess_dt0:
            return ts.dt

        error_order = conf.error_order
        eps_a = min(*conf.eps_a)
        eps_r = min(*conf.eps_r)

        aux = unpack(vec0)
        u0, v0, a0 = aux if unpack.n_arg == 3 else (aux[0], aux[2], aux[3])

        d0 = eval_scaled_norm(u0, eps_a, eps_r)
        d1 = eval_scaled_norm(v0, eps_a, eps_r)
        d2 = eval_scaled_norm(a0, eps_a, eps_r)
        md1d2 = max(d1, d2)

        h0 = 1e-6 if (d0 < 1e-5) or (d1 < 1e-5) else 0.01 * (d0 / d1)
        h1 = (max(1e-6, h0 * 1e-3) if md1d2 < 1e-15
              else (0.01 / md1d2) ** (1 / error_order))

        dt0 = min(100 * h0, h1)
        return dt0

    def __call__(self, ts, vec0, vec1, unpack, **kwargs):
        conf = self.conf
        dt = ts.dt
        fmin, fmax, fsafety, error_order = (
            conf.fmin, conf.fmax, conf.fsafety, conf.error_order,
        )
        aux = unpack(vec0)
        u0, v0, a0 = aux if unpack.n_arg == 3 else (aux[0], aux[2], aux[3])
        aux = unpack(vec1)
        u1, v1, a1 = aux if unpack.n_arg == 3 else (aux[0], aux[2], aux[3])

        u_err, v_err = self.get_scaled_errors(
            dt, vec0, vec1, conf.eps_a, conf.eps_r, unpack,
        )
        emax = max(u_err, v_err)
        status = Struct(u_err=u_err, v_err=v_err, emax=emax)
        if emax <= 1:
            # Step accepted.
            new_dt = dt * min(fmax, fsafety / emax**(1.0 / error_order))
            status.result = 'accept'

        else:
            new_dt = dt * max(fmin, fsafety / emax**(1.0 / error_order))
            status.result = 'reject'

        return new_dt, status

class ElastodynamicsPIDTSC(ElastodynamicsBasicTSC):
    """
    Adaptive time step PID controller for elastodynamics.

    The implementation is based on [1], [2] (PI Controller) and [3] (PID).
    The default parameters correspond to the I-Controller as implemented in
    ``tsc.ed_basic``.

    [1] Grafenhorst, Matthias, Joachim Rang, and Stefan Hartmann.
    “Time-Adaptive Finite Element Simulations of Dynamical Problems for
    Temperature-Dependent Materials.” Journal of Mechanics of Materials and
    Structures 12, no. 1 (November 26, 2016): 57–91.
    https://doi.org/10.2140/jomms.2017.12.57.
    [2] Hairer, Ernst, Syvert Paul Nørsett, and Gerhard Wanner. Solving
    Ordinary Differential Equations II: Stiff and Differential-Algebraic
    Problems. Springer Science & Business Media, 1993.
    [3] Söderlind, Gustaf. “Digital Filters in Adaptive Time-Stepping.” ACM
    Transactions on Mathematical Software 29, no. 1 (March 1, 2003): 1–26.
    https://doi.org/10.1145/641876.641877.
    """
    name = 'tsc.ed_pid'

    _parameters = ElastodynamicsBasicTSC._parameters + [
        ('pcoef', 'float', 0.0, False,
         'Proportional (P) coefficient of the step size control.'),
        ('icoef', 'float', 1.0, False,
         'Intregral (I) coefficient of the step size control.'),
        ('dcoef', 'float', 0.0, False,
         'Derivative (D) coefficient of the step size control.'),
    ]

    def __init__(self, conf, **kwargs):
        ElastodynamicsBasicTSC.__init__(self, conf=conf, **kwargs)

        self.emax0 = 1.0
        self.emax00 = 1.0

    def __call__(self, ts, vec0, vec1, unpack, **kwargs):
        conf = self.conf
        dt = ts.dt
        fmin, fmax, fsafety, pcoef, icoef, dcoef, error_order = (
            conf.fmin, conf.fmax, conf.fsafety,
            conf.pcoef, conf.icoef, conf.dcoef,
            conf.error_order,
        )
        b1 = -(pcoef + icoef + dcoef) / error_order
        b2 = (pcoef + 2 * dcoef) / error_order
        b3 = -dcoef / error_order
        u_err, v_err = self.get_scaled_errors(
            dt, vec0, vec1, conf.eps_a, conf.eps_r, unpack,
        )
        emax = max(u_err, v_err)

        c1 = 1 if b1 == 0 else emax**b1
        c2 = 1 if b2 == 0 else self.emax0**b2
        c3 = 1 if b3 == 0 else self.emax00**b3
        new_dt = dt * min(fmax, max(fmin, fsafety * c1 * c2 * c3))

        status = Struct(u_err=u_err, v_err=v_err, emax=emax)
        status.result = 'accept' if emax <= 1 else 'reject'

        self.emax00 = self.emax0
        self.emax0 = emax

        return new_dt, status

class ElastodynamicsLinearTSC(ElastodynamicsBasicTSC):
    """
    Adaptive time step controller for elastodynamics and linear problems.

    Simple heuristics around :class:`ElastodynamicsBasicTSC` that increases the
    step size only after a sufficient number of accepted iterations passed and
    the increase is large enough. In particular:

    - Let new_dt be the step size proposed by `tsc.ed_basic` and dt the current
      step size.
    - If the current step is rejected, the count attribute is reset to zero and
      fred * new_dt is returned.
    - If the current step is accepted:

      - If the count is lower than inc_wait, it is incremented and dt is
        returned.
      - Otherwise, if (new_dt / dt) >= min_finc (>= 1), the count
        is reset to zero and new_dt is returned.
      - Else, if (new_dt / dt) < min_finc, dt is returned.
    """
    name = 'tsc.ed_linear'

    _parameters = ElastodynamicsBasicTSC._parameters + [
        ('fred', 'float', 1.0, False,
         'Additional step size reduction factor w.r.t. `tsc.ed_basic`.'),
        ('inc_wait', 'int', 10, False,
         """The number of consecutive accepted steps to wait before increasing
            the step size."""),
        ('min_finc', 'float >= 1', 1.5, False,
         'Minimum step size increase factor.'),
    ]

    def __init__(self, conf, **kwargs):
        ElastodynamicsBasicTSC.__init__(self, conf=conf, **kwargs)

        if self.conf.min_finc < 1.0:
            raise ValueError(
                f'min_finc must be >= 1! (is {conf.min_finc})'
            )
        self.count = 0

    def __call__(self, ts, vec0, vec1, unpack, **kwargs):
        conf = self.conf
        dt = ts.dt

        new_dt, status = ElastodynamicsBasicTSC.__call__(
            self, ts, vec0, vec1, unpack, **kwargs
        )
        if status.result == 'reject':
            self.count = 0
            new_dt *= conf.fred

        else: # accept
            if self.count >= conf.inc_wait:
                if (new_dt / dt) >= conf.min_finc:
                    self.count = 0

                else:
                    new_dt = dt

            else:
                self.count += 1
                new_dt = dt

        return new_dt, status
