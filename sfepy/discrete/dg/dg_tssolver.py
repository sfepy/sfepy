import numpy as nm
from numpy import dot
import matplotlib.pyplot as plt
from numpy import newaxis as nax
import numpy.linalg as nla


# sfepy imports
from sfepy.base.base import (get_default, output, assert_,
                             Struct, IndexedStruct)
from sfepy.solvers import SolverMeta, TimeSteppingSolver
from sfepy.solvers.ts import TimeStepper, VariableTimeStepper
from sfepy.solvers.ts_solvers import standard_ts_call

class TSSolver:

    def __init__(self, eq, ic, bc, limiter, basis):
        self.equation = eq
        self.mesh = eq.mesh
        self.basis = basis
        self.limiter = limiter
        self.initial_cond = self.sampleIC(self.mesh, ic, self.intGauss2, self.basis)
        self.boundary_cond = bc

    def initialize(self, t0, tend, tsteps):
        dt = float(tend - t0) / tsteps
        dx = nm.max(self.mesh.coors[1:] - self.mesh.coors[:-1])
        dtdx = dt / dx
        maxa = abs(nm.max(self.equation.terms[1].a(self.mesh.coors)))
        print("Space divided into {0} cells, {1} steps, step size is {2}".format(self.mesh.n_el, len(self.mesh.coors),
                                                                                 dx))
        print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tsteps - 1, tsteps, dt))
        print("Courant number c = max(abs(u)) * dt/dx = {0}".format(maxa * dtdx))
        A = nm.zeros((2, self.mesh.n_el, self.mesh.n_el), dtype=nm.float64)
        b = nm.zeros((2, self.mesh.n_el, 1), dtype=nm.float64)
        u = nm.zeros((2, self.mesh.n_el + 2, tsteps, 1), dtype=nm.float64)

        # bc
        u[:, 0, 0] = self.boundary_cond["left"]
        u[:, -1, 0] = self.boundary_cond["right"]
        # ic
        u[:, 1:-1, 0] = self.initial_cond
        return A, b, dt, u

    # Move to problem class?
    def sampleIC(self, mesh, ic, quad, basis):
        sic = nm.zeros((2, self.mesh.n_el, 1), dtype=nm.float64)

        c = (mesh.coors[1:] + mesh.coors[:-1])/2  # center
        s = (mesh.coors[1:] - mesh.coors[:-1])/2  # scale
        sic[0, :] = quad(lambda t: ic(c + t*s))/2
        sic[1, :] = 3*quad(lambda t: t*ic(c + t*s))/2
        return sic

    # Will be taken care of in Integral class
    @staticmethod
    def intGauss2(f):

        x_1 = - nm.sqrt(1./3.)
        x_2 = nm.sqrt(1./3.)

        return f(x_1) + f(x_2)

    @staticmethod
    def intGauss3(f):
        x_0 = 0
        x_1 = - nm.sqrt(3./5.)
        x_2 = nm.sqrt(3./5.)

        w_0 = 8./9.
        w_1 = 5. / 9.

        return w_0 * f(x_0) + w_1 * f(x_1) + w_1 * f(x_2)

    # Separate class for limiters?
    @staticmethod
    def moment_limiter(u):
        """
        Krivodonova(2007): Limiters for high-order discontinuous Galerkin methods

        :param u: solution at time step n
        :return: limited solution
        """

        def minmod(a, b, c):
            seq = (nm.sign(a) == nm.sign(b)) & (nm.sign(b) == nm.sign(c))

            res = nm.zeros(nm.shape(a))
            res[seq] = nm.sign(a[seq]) * nm.minimum.reduce([nm.abs(b[seq]),
                                                            nm.abs(a[seq]),
                                                            nm.abs(c[seq])])

            return res

        idx = nm.arange(nm.shape(u[0, 1:-1])[0])
        nu = nm.copy(u)
        for l in range(1, 0, -1):
            tilu = minmod(nu[l, 1:-1][idx],
                          nu[l-1, 2:][idx] - nu[l-1, 1:-1][idx],
                          nu[l-1, 1:-1][idx] - nu[l-1, :-2][idx])
            idx = tilu != nu
            nu[l, 1:-1][idx] = tilu[idx]
        return nu

    def solve(self, t0, tend, tsteps=10):
        raise NotImplemented


class RK3Solver(TSSolver):
    """
    Runge-Kutta of order 3, with limiter
    """

    def solve(self, t0, tend, tsteps=10):

        A, b, dt, u = self.initialize(t0, tend, tsteps)

        # setup RK3 specific arrays
        u1 = nm.zeros((2, self.mesh.n_el + 2, 1), dtype=nm.float64)
        u2 = nm.zeros((2, self.mesh.n_el + 2, 1), dtype=nm.float64)

        for it in range(1, tsteps):
            # ----1st stage----
            # bcs
            u1[:, 0] = self.boundary_cond["left"]
            u1[:, -1] = self.boundary_cond["right"]

            # get RHS
            A[:] = 0
            b[:] = 0
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u[:, :, it-1])

            # get update u1
            # maybe use: for more general cases
            #                                     dot(nm.linalg.inv(A[0]), b[0])
            #                                     dot(nm.linalg.inv(A[1]), b[1])
            u1[0, 1:-1] = u[0, 1:-1, it-1] + dt * b[0] / nm.diag(A[0])[:, nax]
            u1[1, 1:-1] = u[1, 1:-1, it-1] + dt * b[1] / nm.diag(A[1])[:, nax]

            # limit
            u1 = self.limiter(u1)

            # ----2nd stage----
            # bcs
            u2[:, 0] = self.boundary_cond["left"]
            u2[:, -1] = self.boundary_cond["right"]

            # get RHS
            A[:] = 0
            b[:] = 0
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u1[:, :])

            # get update u2
            u2[0, 1:-1] = (3 * u[0, 1:-1, it - 1] + u1[0, 1:-1]
                           + dt * b[0] / nm.diag(A[0])[:, nax]) / 4
            u2[1, 1:-1] = (3 * u[1, 1:-1, it - 1] + u1[1, 1:-1]
                           + dt * b[1] / nm.diag(A[1])[:, nax]) / 4

            # limit
            u2 = self.limiter(u2)

            # ----3rd stage-----
            # get RHS
            A[:] = 0
            b[:] = 0
            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var=None, u=u2[:, :])

            # get update u3
            u[0, 1:-1, it] = (u[0, 1:-1, it - 1] + 2 * u2[0, 1:-1]
                              + 2*dt * b[0] / nm.diag(A[0])[:, nax]) / 3
            u[1, 1:-1, it] = (u[1, 1:-1, it - 1] + 2 * u2[1, 1:-1]
                              + 2*dt * b[1] / nm.diag(A[1])[:, nax]) / 3

            # limit
            u[:, :, it] = self.limiter(u[:, :, it])

        return u, dt

class EUSolver(TSSolver):
    """
    Euler method with limiter
    """

    def solve(self, t0, tend, tsteps=10):
        """

        :param t0:
        :param tend:
        :param tsteps:
        :return:
        """
        A, b, dt, u = self.initialize(t0, tend, tsteps)
        self.equation.terms[1].get_state_variables()[0].setup_dof_info()
        di = self.equation.terms[1].get_state_variables()[0].di
        self.equation.terms[1].get_state_variables()[0].setup_initial_conditions(self.ics)

        for it in range(1, tsteps):
            A[:] = 0
            b[:] = 0

            self.equation.evaluate(dw_mode="matrix", asm_obj=A, diff_var="u")
            self.equation.evaluate(dw_mode="vector", asm_obj=b, diff_var="u")

            u[0, 1:-1, it] = u[0, 1:-1, it - 1] + dt * b[0] / nm.diag(A[0])[:, nax]
            u[1, 1:-1, it] = u[1, 1:-1, it - 1] + dt * b[1] / nm.diag(A[1])[:, nax]

            u[:, :, it] = self.limiter(u[:, :, it])

        return u, dt

class DGTimeSteppingSolver(TimeSteppingSolver):
    """
    Explicit time stepping solver with a fixed time step.

    # TODO maybe inherit directly from SimpleTimeSteppingSolver and override solve_step method
    # TODO or create metaclass ExplicitTimeSteppingSolver
    """
    name = 'ts.euler'

    __metaclass__ = SolverMeta

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
    ]

    def __init__(self, conf, nls=None, context=None, **kwargs):

        TimeSteppingSolver.__init__(self, conf, nls=nls, context=context,
                                    **kwargs)
        self.ts = TimeStepper.from_conf(self.conf)

        nd = self.ts.n_digit
        format = '====== time %%e (step %%%dd of %%%dd) =====' % (nd, nd)

        self.format = format
        self.verbose = self.conf.verbose

    def solve_step0(self, nls, vec0):
        res = nls.fun(vec0)
        err = nm.linalg.norm(res)
        output('initial residual: %e' % err, verbose=self.verbose)
        vec = vec0.copy()

        return vec

    def solve_step(self, ts, nls, vec, prestep_fun=None):
        return nls(vec, ts=ts)

    def output_step_info(self, ts):
        output(self.format % (ts.time, ts.step + 1, ts.n_step),
               verbose=self.verbose)

    @standard_ts_call
    def __call__(self, vec0=None, nls=None, init_fun=None, prestep_fun=None,
                 poststep_fun=None, status=None, **kwargs):
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

            vect = self.solve_step(ts, nls, vec, prestep_fun)

            poststep_fun(ts, vect)

            vec = vect

        return vec


from sfepy.solvers.solvers import SolverMeta, NonlinearSolver
from sfepy.base.log import Log, get_logging_conf


class EulerStepSolver(NonlinearSolver):
    """
    Not actually nonlinear solver, solves only linear system.
    Updates solution using euler method
    # TODO create ExplicitStepSolver(?) class, inherit from it
    # - unify structure of __call__ method, something like:
    #  1. prepare data
    #  2. call method computing all the stages
    #  3. provide stats for status, outputs
    #  4. return
    """
    name = 'sls.euler'
    __metaclass__ = SolverMeta
    _parameters = []

    def __init__(self, conf, post_stage_hook=lambda x: x, **kwargs):
        NonlinearSolver.__init__(self, conf, **kwargs)

        conf = self.conf

        log = get_logging_conf(conf)
        conf.log = log = Struct(name='log_conf', **log)
        conf.is_any_log = (log.text is not None) or (log.plot is not None)

        if conf.is_any_log:
            self.log = Log([[r'$||r||$'], ['iteration']],
                           xlabels=['', 'all iterations'],
                           ylabels=[r'$||r||$', 'iteration'],
                           yscales=['log', 'linear'],
                           is_plot=conf.log.plot is not None,
                           log_filename=conf.log.text,
                           formats=[['%.8e'], ['%d']])

        else:
            self.log = None
        self.post_stage_hook = post_stage_hook

    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None, ts=None):
        if ts is None:
            raise ValueError("Provide TimeStepper to explicit Euler solver")

        conf = get_default(conf, self.conf)
        fun = get_default(fun, self.fun)
        fun_grad = get_default(fun_grad, self.fun_grad)
        lin_solver = get_default(lin_solver, self.lin_solver)
        iter_hook = get_default(iter_hook, self.iter_hook)
        status = get_default(status, self.status)

        ls_eps_a, ls_eps_r = lin_solver.get_tolerance()
        eps_a = get_default(ls_eps_a, 1.0)
        eps_r = get_default(ls_eps_r, 1.0)

        vec_x = vec_x0.copy()

        vec_r = fun(vec_x)

        mtx_a = fun_grad(vec_x)
        ls_status = {}
        vec_dx = lin_solver(vec_r, x0=vec_x,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_e = mtx_a * vec_dx - vec_r
        lerr = nla.norm(vec_e)
        output('linear system sol error {}'.format(lerr))
        output('mtx max {}, min {}, trace {}'
               .format(mtx_a.max(), mtx_a.min(), nm.sum(mtx_a.diagonal())))

        vec_x = vec_x - ts.dt * (vec_dx - vec_x)
        vec_x = self.post_stage_hook(vec_x)

        return vec_x


class RK3StepSolver(NonlinearSolver):
    """
    3rd order Runge-Kutta method

    # TODO create ExplicitStepSolver(?) class, inherit from it

    """

    name = 'sls.runge_kutta_3'
    __metaclass__ = SolverMeta
    _parameters = []

    def __init__(self, conf, post_stage_hook=lambda x: x, **kwargs):
        NonlinearSolver.__init__(self, conf, **kwargs)

        conf = self.conf

        log = get_logging_conf(conf)
        conf.log = log = Struct(name='log_conf', **log)
        conf.is_any_log = (log.text is not None) or (log.plot is not None)

        if conf.is_any_log:
            self.log = Log([[r'$||r||$'], ['iteration']],
                           xlabels=['', 'all iterations'],
                           ylabels=[r'$||r||$', 'iteration'],
                           yscales=['log', 'linear'],
                           is_plot=conf.log.plot is not None,
                           log_filename=conf.log.text,
                           formats=[['%.8e'], ['%d']])

        else:
            self.log = None

        self.post_stage_hook = post_stage_hook

    def __call__(self, vec_x0, conf=None, fun=None, fun_grad=None,
                 lin_solver=None, iter_hook=None, status=None, ts=None):
        if ts is None:
            raise ValueError("Provide TimeStepper to explicit Runge-Kutta solver")

        conf = get_default(conf, self.conf)
        fun = get_default(fun, self.fun)
        fun_grad = get_default(fun_grad, self.fun_grad)
        lin_solver = get_default(lin_solver, self.lin_solver)
        iter_hook = get_default(iter_hook, self.iter_hook)
        status = get_default(status, self.status)

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

        vec_x1 = vec_x - ts.dt * (vec_dx - vec_x)

        from dg_field import get_unraveler, get_raveler
        unravel = get_unraveler(3, 99)
        un_vec_r = unravel(vec_r)
        un_vec_x1 = unravel(vec_x1)
        un_vec_x0 = unravel(vec_x0)

        vec_x1 = self.post_stage_hook(vec_x1)

        un_vec_x1_lim = unravel(vec_x1)

        # ----2nd stage----
        vec_r = fun(vec_x1)
        mtx_a = fun_grad(vec_x1)
        vec_dx = lin_solver(vec_r, x0=vec_x1,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x2 = (3 * vec_x + vec_x1 - ts.dt * (vec_dx - vec_x1))/4
        un_vec_x2 = unravel(vec_x2)
        vec_x2 = self.post_stage_hook(vec_x2)
        un_vec_x2_lim = unravel(vec_x2)


        # ----3rd stage-----
        vec_r = fun(vec_x2)
        mtx_a = fun_grad(vec_x2)
        vec_dx = lin_solver(vec_r, x0=vec_x2,
                            eps_a=eps_a, eps_r=eps_r, mtx=mtx_a,
                            status=ls_status)

        vec_x3 = (vec_x + 2 * vec_x2 - 2*ts.dt * (vec_dx - vec_x2))/3

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


