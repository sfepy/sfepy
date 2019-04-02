import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt
import numpy.linalg as nla

# sfepy imports
from sfepy.base.base import (get_default, output, assert_,
                             Struct, IndexedStruct)
from sfepy.solvers import SolverMeta, TimeSteppingSolver
from sfepy.solvers.ts import TimeStepper, VariableTimeStepper
from sfepy.solvers.ts_solvers import standard_ts_call
from sfepy.solvers.solvers import SolverMeta, NonlinearSolver
from sfepy.base.log import Log, get_logging_conf

class DGMultiStageTS(TimeSteppingSolver):
    """
    Explicit time stepping solver with multistage solve_step
    """
    __metaclass__ = SolverMeta
    name = "ts.multistaged"

    _parameters = [
        ('t0', 'float', 0.0, False,
         'The initial time.'),
        ('t1', 'float', 1.0, False,
         'The final time.'),
        ('dt', 'float', None, False,
         'The time step. Used if `n_step` is not given.'),
        ('n_step', 'int', 10, False,
         'The number of time steps. Has precedence over `dt`.'),
        ('quasistatic', 'bool', False, False,
         """If True, assume a quasistatic time-stepping. Then the non-linear
            solver is invoked also for the initial time."""),
        ('limiter', 'function', None, None,
         "Limiter for DG FEM"),
    ]

    def __init__(self, conf, nls=None, context=None, **kwargs):

        TimeSteppingSolver.__init__(self, conf, nls=nls, context=context,
                                    **kwargs)
        self.ts = TimeStepper.from_conf(self.conf)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format
        self.verbose = self.conf.verbose

        if hasattr(conf, "limiter"):
            if conf.limiter is not None:
                self.post_stage_hook = conf.limiter
        if hasattr(conf, "post_stage_hook"):
            if conf.post_stage_hook is not None:
                self.post_stage_hook = conf.post_stage_hook

        if "post_stage_hook" in kwargs.keys():
            self.post_stage_hook = kwargs["post_stage_hook"]
        else:
            self.post_stage_hook = lambda x: x

    def solve_step0(self, nls, vec0):
        res = nls.fun(vec0)
        err = nm.linalg.norm(res)
        output('initial residual: %e' % err, verbose=self.verbose)
        vec = vec0.copy()

        return vec

    def solve_step(self, ts, nls, vec, prestep_fun=None,  poststep_fun=None, status=None):
        raise NotImplementedError("Called abstract solver, call subclass.")

    def output_step_info(self, ts):
        output(self.format % (ts.time, ts.step + 1, ts.n_step),
               verbose=self.verbose)

    @standard_ts_call
    def __call__(self, vec0=None, nls=None, init_fun=None, prestep_fun=None,
                 poststep_fun=None, status=None):
        """
        Solve the time-dependent problem.
        """
        ts = self.ts
        nls = get_default(nls, self.nls)

        vec0 = init_fun(ts, vec0)

        self.output_step_info(ts)
        if ts.step == 0:
            prestep_fun(ts, vec0)

            vec = self.solve_step0(nls, vec0)

            poststep_fun(ts, vec)
            ts.advance()

        else:
            vec = vec0

        for step, time in ts.iter_from(ts.step):
            self.output_step_info(ts)

            prestep_fun(ts, vec)

            vect = self.solve_step(ts, nls, vec, prestep_fun, poststep_fun, status)

            poststep_fun(ts, vect)

            vec = vect

        return vec


class EulerStepSolver(DGMultiStageTS):
    """
    Updates solution using euler method
    # - unify structure of __call__ method, something like:
    #  1. prepare data
    #  2. call method computing all the stages
    #  3. provide stats for status, outputs
    #  4. return
    """
    name = 'ts.euler'
    __metaclass__ = SolverMeta

    def solve_step(self, ts, nls, vec_x0, status=None,
                   prestep_fun=None, poststep_fun=None):
        if ts is None:
            raise ValueError("Provide TimeStepper to explicit Euler solver")

        conf = nls.conf
        fun = nls.fun
        fun_grad = nls.fun_grad
        lin_solver = nls.lin_solver
        iter_hook = nls.iter_hook
        status = get_default(status, nls.status)

        ls_eps_a, ls_eps_r = lin_solver.get_tolerance()
        eps_a = get_default(ls_eps_a, 1.0)
        eps_r = get_default(ls_eps_r, 1.0)

        vec_x = vec_x0.copy()

        vec_r = fun(vec_x)

        mtx_a = fun_grad(vec_x)
        ls_status = {}

        vec_dx = lin_solver(vec_r, x0=vec_x,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status) - vec_x

        vec_e = mtx_a * vec_dx - vec_r
        lerr = nla.norm(vec_e)
        output('linear system sol error {}'.format(lerr))
        output('mtx max {}, min {}, trace {}'
               .format(mtx_a.max(), mtx_a.min(), nm.sum(mtx_a.diagonal())))

        vec_x = vec_x + ts.dt * vec_dx
        vec_x = self.post_stage_hook(vec_x)

        return vec_x


class TVDRK3StepSolver(DGMultiStageTS):
    """
    3rd order Runge-Kutta method
    """

    name = 'ts.tvd_runge_kutta_3'
    __metaclass__ = SolverMeta

    def solve_step(self, ts, nls, vec_x0, status=None,
                   prestep_fun=None, poststep_fun=None):
        if ts is None:
            raise ValueError("Provide TimeStepper to explicit Runge-Kutta solver")

        conf = nls.conf
        fun = nls.fun
        fun_grad = nls.fun_grad
        lin_solver = nls.lin_solver
        iter_hook = nls.iter_hook
        status = get_default(status, nls.status)

        ls_eps_a, ls_eps_r = lin_solver.get_tolerance()
        eps_a = get_default(ls_eps_a, 1.0)
        eps_r = get_default(ls_eps_r, 1.0)
        ls_status = {}

        # Add pre-stage hook?
        # ----1st stage----
        vec_x = vec_x0.copy()

        vec_r = fun(vec_x)
        mtx_a = fun_grad(vec_x)
        full_mtx_a = mtx_a.toarray()
        vec_dx = lin_solver(vec_r, x0=vec_x,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x1 = vec_x + ts.dt * (vec_dx - vec_x)

        from dg_field import get_unraveler, get_raveler
        unravel = get_unraveler(3, 99)
        un_vec_r = unravel(vec_r)
        un_vec_x1 = unravel(vec_x1)
        un_vec_x0 = unravel(vec_x0)

        vec_x1 = self.post_stage_hook(vec_x1)

        un_vec_x1_lim = unravel(vec_x1)

        # ----2nd stage----
        # ts.set_substep_time(1./2. * ts.dt)
        # prestep_fun(ts, vec_x1)
        vec_r = fun(vec_x1)
        mtx_a = fun_grad(vec_x1)
        vec_dx = lin_solver(vec_r, x0=vec_x1,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x2 = (3 * vec_x + vec_x1 + ts.dt * (vec_dx - vec_x1))/4

        un_vec_x2 = unravel(vec_x2)

        vec_x2 = self.post_stage_hook(vec_x2)

        un_vec_x2_lim = unravel(vec_x2)

        # ----3rd stage-----
        ts.set_substep_time(1./2. * ts.dt)
        prestep_fun(ts, vec_x2)
        vec_r = fun(vec_x2)
        mtx_a = fun_grad(vec_x2)
        vec_dx = lin_solver(vec_r, x0=vec_x2,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x3 = (vec_x + 2 * vec_x2 + 2*ts.dt * (vec_dx - vec_x2))/3

        un_vec_x3 = unravel(vec_x3)
        vec_x3 = self.post_stage_hook(vec_x3)
        un_vec_x3_lim = unravel(vec_x3)

        # vec_e = mtx_a * vec_dx - vec_r
        # lerr = nla.norm(vec_e)
        # output('linear system sol error {}'.format(lerr))
        # output('mtx max {}, min {}, trace {}'
        #        .format(mtx_a.max(), mtx_a.min(), nm.sum(mtx_a.diagonal())))
        # vec_x -= ts.dt * vec_dx

        return vec_x3


class RK4StepSolver(DGMultiStageTS):
    """
    Based on Hesthaven, J. S., & Warburton, T. (2008). Nodal Discontinuous Galerkin Methods.
    Journal of Physics A: Mathematical and Theoretical (Vol. 54). New York, NY: Springer New York.
    http://doi.org/10.1007/978-0-387-72067-8

    p. 63

    """
    name = 'ts.runge_kutta_4'
    __metaclass__ = SolverMeta

    def solve_step(self, ts, nls, vec_x0, status=None,
                   prestep_fun=None, poststep_fun=None):
        if ts is None:
            raise ValueError("Provide TimeStepper to explicit Runge-Kutta solver")

        from dg_field import get_unraveler, get_raveler
        unravel = get_unraveler(3, 99)

        conf = nls.conf
        fun = nls.fun
        fun_grad = nls.fun_grad
        lin_solver = nls.lin_solver
        iter_hook = nls.iter_hook
        status = get_default(status, self.status)

        ls_eps_a, ls_eps_r = lin_solver.get_tolerance()
        eps_a = get_default(ls_eps_a, 1.0)
        eps_r = get_default(ls_eps_r, 1.0)
        ls_status = {}

        # ----1st stage----
        vec_r = fun(vec_x0)
        mtx_a = fun_grad(vec_x0)
        vec_dx = lin_solver(vec_r, x0=vec_x0,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x1 = vec_dx - vec_x0

        # for debugging
        full_mtx_a = mtx_a.toarray()
        un_vec_r = unravel(vec_r)
        un_vec_x1 = unravel(vec_x1)
        un_vec_x0 = unravel(vec_x0)

        vec_x1 = self.post_stage_hook(vec_x1)

        un_vec_x1_lim = unravel(vec_x1)

        # ----2nd stage----
        vec_r = fun(vec_x0 + 1./2. * ts.dt * vec_x1)
        mtx_a = fun_grad(vec_x0 + 1./2. * ts.dt * vec_x1)
        vec_dx = lin_solver(vec_r, #x0=vec_x0 + 1./2. * ts.dt * vec_x1,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x2 = vec_dx - vec_x1

        un_vec_x2 = unravel(vec_x2)

        vec_x2 = self.post_stage_hook(vec_x2)

        un_vec_x2_lim = unravel(vec_x2)

        # ----3rd stage-----

        vec_r = fun(vec_x0 + 1./2. * ts.dt * vec_x2)
        mtx_a = fun_grad(vec_x0 + 1./2. * ts.dt * vec_x2)
        vec_dx = lin_solver(vec_r, #x0=vec_x0 + 1./2. * ts.dt * vec_x2,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x3 = vec_dx - vec_x2

        un_vec_x3 = unravel(vec_x3)

        vec_x3 = self.post_stage_hook(vec_x3)

        un_vec_x3_lim = unravel(vec_x3)

        # ----4th stage-----
        vec_r = fun(vec_x0 + ts.dt * vec_x3)
        mtx_a = fun_grad(vec_x0 + ts.dt * vec_x3)
        vec_dx = lin_solver(vec_r, #x0=vec_x0 + ts.dt * vec_x3,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x4 = vec_dx - vec_x3

        un_vec_x4 = unravel(vec_x4)

        vec_x4 = self.post_stage_hook(vec_x4)

        un_vec_x4_lim = unravel(vec_x4)

        vec_fin = vec_x0 + 1./6. * ts.dt * (vec_x1 + 2 * vec_x2 + 2 * vec_x3 + vec_x4)

        return vec_fin

