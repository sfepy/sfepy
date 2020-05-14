import numpy as nm
import matplotlib.pyplot as plt
from os.path import join as pjoin

# sfepy imports
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct, configure_output, output
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver

from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm
from sfepy.base.ioutils import ensure_path

# DG imports
from sfepy.terms.terms_dg import AdvectionDGFluxTerm
from sfepy.solvers.ts_dg_solvers import TVDRK3StepSolver, RK4StepSolver, EulerStepSolver
from sfepy.discrete.dg.fields import DGField
from sfepy.discrete.dg.limiters import IdentityLimiter, MomentLimiter1D

from sfepy.discrete.dg.utils.inits_consts import \
    left_par_q, gsmooth, const_u, ghump, superic
from sfepy.discrete.dg.utils.visualizer import load_1D_vtks, plot1D_DG_sol
from run_dg_utils import clear_folder

# vvvvvvvvvvvvvvvv#
approx_order = 1
# ^^^^^^^^^^^^^^^^#
# Setup  names
outputs_folder = "../outputs"

domain_name = "domain_1D"
problem_name = "iadv_1D"
output_folder = pjoin(outputs_folder, problem_name, str(approx_order))
output_format = "vtk"
mesh_output_folder = "../mesh"
save_timestn = 100
clear_folder(pjoin(output_folder, "*" + output_format))
configure_output({'output_screen': True,
                  'output_log_name':
                      pjoin(output_folder, f"last_run_{problem_name}_{approx_order}.txt")})

# ------------
# | Get mesh |
# ------------
X1 = 0.
XN = 1.
n_nod = 100
n_el = n_nod - 1
coors = nm.linspace(X1, XN, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('uniform_1D{}'.format(n_nod), coors, None,
                      [conn], [mat_ids], descs)

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

velo = 1.0
a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, u, v)", "a.val, u[-1], v",
                                  integral, omega,
                                  u=u, v=v, a=a)

alpha = Material('alpha', val=[.0])
FluxT = AdvectionDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
                            integral, omega, u=u, v=v, a=a, alpha=alpha)

eq = Equation('balance', MassT - StiffT + FluxT)
eqs = Equations([eq])

# ------------------------------
# | Create bounrady conditions |
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
pb = Problem(problem_name, equations=eqs, conf=Struct(options={"save_times": save_timestn}, ics={},
                                                      ebcs={}, epbcs={}, lcbcs={}, materials={}),
             active_only=False)
pb.setup_output(output_dir=output_folder, output_format=output_format)
pb.set_ics(Conditions([ics]))

# ------------------
# | Create limiter |
# ------------------
limiter = IdentityLimiter

# ---------------------------
# | Set time discretization |
# ---------------------------
CFL = .1
max_velo = nm.max(nm.abs(velo))
t0 = 0
t1 = 1
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

tss = EulerStepSolver(tss_conf,
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
from sfepy.discrete.dg.utils.plot_1D_dg import load_and_plot_fun

load_and_plot_fun(output_folder, domain_name, t0, t1, min(tn, save_timestn), ic_fun)
