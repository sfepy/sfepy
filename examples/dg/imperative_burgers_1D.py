"""
Burgers equation in 2D solved using discontinous galerkin method
"""

from os.path import join as pjoin

import numpy as nm

from examples.dg.example_dg_common import clear_folder, get_gen_1D_mesh_hook
from examples.dg.inits_consts import ghump
from script.dg_plot_1D import load_and_plot_fun

# sfepy imports
from sfepy.base.base import IndexedStruct
from sfepy.base.base import Struct, configure_output, output
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions
from sfepy.discrete.dg.fields import DGField
from sfepy.discrete.dg.limiters import MomentLimiter1D
from sfepy.discrete.fem import FEDomain
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_dg_solvers import TVDRK3StepSolver
from sfepy.terms.terms_dg import NonlinearHyperbolicDGFluxTerm, \
    NonlinearScalarDotGradTerm
from sfepy.terms.terms_dot import DotProductVolumeTerm

def main():
    # vvvvvvvvvvvvvvvv #
    approx_order = 1
    # ^^^^^^^^^^^^^^^^ #

    # Setup output names
    outputs_folder = "../outputs"

    domain_name = "domain_1D"
    problem_name = "iburgers_1D"
    output_folder = pjoin(outputs_folder, problem_name, str(approx_order))
    output_format = "vtk"
    save_timestn = 100
    clear_folder(pjoin(output_folder, "*." + output_format))
    configure_output({'output_screen': True,
                      'output_log_name':
                       pjoin(output_folder,
                             f"last_run_{problem_name}_{approx_order}.txt")})

    # ------------
    # | Get mesh |
    # ------------
    X1 = 0.
    XN = 1.
    n_nod = 100
    n_el = n_nod - 1
    mesh = get_gen_1D_mesh_hook(X1, XN, n_nod).read(None)

    # -----------------------------
    # | Create problem components |
    # -----------------------------

    integral = Integral('i', order=approx_order * 2)
    domain = FEDomain(domain_name, mesh)
    omega = domain.create_region('Omega', 'all')
    left = domain.create_region('Gamma1',
                                'vertices in x == %.10f' % X1,
                                'vertex')
    right = domain.create_region('Gamma2',
                                 'vertices in x == %.10f' % XN,
                                 'vertex')
    field = DGField('dgfu', nm.float64, 'scalar', omega,
                    approx_order=approx_order)

    u = FieldVariable('u', 'unknown', field, history=1)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u",
                                 integral, omega, u=u, v=v)

    velo = nm.array(1.0)


    def adv_fun(u):
        vu = velo.T * u[..., None]
        return vu


    def adv_fun_d(u):
        v1 = velo.T * nm.ones(u.shape + (1,))
        return v1


    burg_velo = velo.T / nm.linalg.norm(velo)


    def burg_fun(u):
        vu = burg_velo * u[..., None] ** 2
        return vu


    def burg_fun_d(u):
        v1 = 2 * burg_velo * u[..., None]
        return v1


    StiffT = NonlinearScalarDotGradTerm("burgers_stiff(f, df, u, v)",
                                        "fun , fun_d, u[-1], v",
                                        integral, omega,
                                        u=u, v=v,
                                        fun=burg_fun, fun_d=burg_fun_d)

    alpha = Material('alpha', val=[.0])
    # FluxT = AdvectDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
    #                          integral, omega, u=u, v=v, a=a, alpha=alpha)

    FluxT = NonlinearHyperbolicDGFluxTerm("burgers_lf_flux(f, df, u, v)",
                                          "fun , fun_d, v, u[-1]",
                                          integral, omega,
                                          u=u, v=v,
                                          fun=burg_fun, fun_d=burg_fun_d)

    eq = Equation('balance', MassT - StiffT + FluxT)
    eqs = Equations([eq])

    # ------------------------------
    # | Create boundary conditions |
    # ------------------------------
    left_fix_u = EssentialBC('left_fix_u', left, {'u.all': 1.0})
    right_fix_u = EssentialBC('right_fix_u', right, {'u.all': 0.0})


    # ----------------------------
    # | Create initial condition |
    # ----------------------------
    def ic_wrap(x, ic=None):
        return ghump(x - .3)


    ic_fun = Function('ic_fun', ic_wrap)
    ics = InitialCondition('ic', omega, {'u.0': ic_fun})

    # ------------------
    # | Create problem |
    # ------------------
    pb = Problem(problem_name,
                 equations=eqs,
                 conf=Struct(options={"save_times": save_timestn},
                             ics={}, ebcs={}, epbcs={}, lcbcs={},
                             materials={}),
                 active_only=False)
    pb.setup_output(output_dir=output_folder, output_format=output_format)
    pb.set_ics(Conditions([ics]))

    # ------------------
    # | Create limiter |
    # ------------------
    limiter = MomentLimiter1D

    # ---------------------------
    # | Set time discretization |
    # ---------------------------
    CFL = .2
    max_velo = nm.max(nm.abs(velo))
    t0 = 0
    t1 = .2
    dx = nm.min(mesh.cmesh.get_volumes(1))
    dt = dx / max_velo * CFL / (2 * approx_order + 1)
    tn = int(nm.ceil((t1 - t0) / dt))
    dtdx = dt / dx

    # ------------------
    # | Create solver |
    # ------------------
    ls = ScipyDirect({})
    nls_status = IndexedStruct()
    nls = Newton({'is_linear': True}, lin_solver=ls, status=nls_status)

    tss_conf = {'t0'     : t0,
                't1'     : t1,
                'n_step' : tn,
                'limiters': {"dgfu": limiter}}

    tss = TVDRK3StepSolver(tss_conf,
                           nls=nls, context=pb, verbose=True)

    # ---------
    # | Solve |
    # ---------
    pb.set_solver(tss)
    state_end = pb.solve()

    output("Solved equation \n\n\t\t u_t - div(f(u))) = 0\n")
    output(f"With IC: {ic_fun.name}")
    # output("and EBCs: {}".format(pb.ebcs.names))
    # output("and EPBCS: {}".format(pb.epbcs.names))
    output("-------------------------------------")
    output(f"Approximation order is {approx_order}")
    output(f"Space divided into {mesh.n_el} cells, " +
           f"{len(mesh.coors)} steps, step size is {dx}")
    output(f"Time divided into {tn - 1} nodes, {tn} steps, step size is {dt}")
    output(f"CFL coefficient was {CFL} and " +
          f"order correction {1 / (2 * approx_order + 1)}")
    output(f"Courant number c = max(abs(u)) * dt/dx = {max_velo * dtdx}")
    output("------------------------------------------")
    output(f"Time stepping solver is {tss.name}")
    output(f"Limiter used: {limiter.name}")
    output("======================================")

    # ----------
    # | Plot 1D|
    # ----------
    load_and_plot_fun(output_folder, domain_name,
                      t0, t1, min(tn, save_timestn),
                      ic_fun)


if __name__ == '__main__':
    main()
