import numpy as nm
from numpy.linalg import norm
import matplotlib.pyplot as plt
from os.path import join as pjoin

# sfepy imports
from discrete.fem.periodic import match_x_line, match_y_line
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions, PeriodicBC,  DGPeriodicBC, DGEssentialBC
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import VTKMeshIO
from sfepy.base.ioutils import ensure_path

from sfepy.discrete.variables import DGFieldVariable

from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm

# local imports
from sfepy.discrete.dg.dg_terms import AdvectionDGFluxTerm
from sfepy.discrete.dg.dg_tssolver \
    import EulerStepSolver, TVDRK3StepSolver
from sfepy.discrete.dg.fields import DGField
from sfepy.discrete.dg.limiters import IdentityLimiter, MomentLimiter1D

from sfepy.discrete.dg.my_utils.inits_consts \
    import left_par_q, gsmooth, const_u, ghump, superic

from run_dg_utils import clear_folder

# vvvvvvvvvvvvvvvv#
approx_order = 2
CFL = .8
# ^^^^^^^^^^^^^^^^#
# Setup  names
domain_name = "iadv_2D_tens"
problem_name = "iadv_2D_tens"
output_folder = pjoin("output", problem_name, str(approx_order))
output_format = "msh"
mesh_output_folder = "output/mesh"
save_timestn = 100
clear_folder(pjoin(output_folder, "*." + output_format))

# ------------
# | Get mesh |
# -----------
# mesh = gen_block_mesh((1., 1.), (50, 50), (0.5, 0.5))

mesh_name = "tens_2D_mesh20"
mesh = Mesh.from_file("mesh/" + mesh_name + ".vtk")

angle = - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
velo = -nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
velo = nm.array([[1., 1.]]).T
max_velo = nm.max(nm.linalg.norm(velo))

# -----------------------------
# | Create problem components |
# -----------------------------

integral = Integral('i', order=approx_order * 2)
domain = FEDomain(domain_name, mesh)
omega = domain.create_region('Omega', 'all')

left = domain.create_region('left',
                            'vertices in x == 0',
                            'edge')

right = domain.create_region('right',
                             'vertices in x == 1',
                             'edge')

top = domain.create_region('top',
                           'vertices in y == 1',
                           'edge')

bottom = domain.create_region('bottom',
                              'vertices in y == 0',
                              'edge')

field = DGField('dgfu', nm.float64, 'scalar', omega,
                approx_order=approx_order)

u = DGFieldVariable('u', 'unknown', field, history=1)
v = DGFieldVariable('v', 'test', field, primary_var_name='u')

MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u", integral, omega, u=u, v=v)

a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, u, v)", "a.val, u[-1], v", integral, omega,
                                  u=u, v=v, a=a)

alpha = Material('alpha', val=[.0])
FluxT = AdvectionDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
                            integral, omega, u=u, v=v, a=a, alpha=alpha)

eq = Equation('balance', MassT + StiffT - FluxT)
eqs = Equations([eq])

# ------------------------------
# | Create bounrady conditions |
# ------------------------------
dirichlet_bc_u = DGEssentialBC('left_fix_u', left, {'u.all': 1.0})
periodic1_bc_u = DGPeriodicBC('top_bot', [top, bottom], {'u.all': 'u.all'}, match='match_x_line')
periodic2_bc_u = DGPeriodicBC('left_right', [right, left], {'u.all': 'u.all'}, match='match_y_line')


# right_fix_u = EssentialBC('right_fix_u', right, {'u.all' : 0.0})


# ----------------------------
# | Create initial condition |
# ----------------------------
def ic_wrap(x, ic=None):
    return gsmooth(x[..., 0:1]) * gsmooth(x[..., 1:])


ic_fun = Function('ic_fun', ic_wrap)
ics = InitialCondition('ic', omega, {'u.0': ic_fun})

# ------------------
# | Create problem |
# ------------------
pb = Problem(problem_name, equations=eqs, conf=Struct(options={"save_times": save_timestn}, ics={},
                                                      ebcs={}, epbcs={}, lcbcs={}, materials={},
                                                      ),
             active_only=False)
pb.setup_output(output_dir=output_folder, output_format=output_format)
pb.functions = {'match_x_line': Function("match_x_line", match_x_line),
                'match_y_line': Function("match_y_line", match_y_line)}
pb.set_ics(Conditions([ics]))
# pb.set_bcs(
#         # ebcs=Conditions([dirichlet_bc_u]),
#         epbcs=Conditions([
#             periodic1_bc_u,
#             periodic2_bc_u
#         ])
#     )

# ------------------
# | Create limiter |
# ------------------
limiter = IdentityLimiter

# ---------------------------
# | Set time discretization |
# ---------------------------
max_velo = nm.max(nm.linalg.norm(velo))
t0 = 0
t1 = 1
dx = nm.min(mesh.cmesh.get_volumes(2))
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
            "limiter": limiter}

tss = EulerStepSolver(tss_conf,
                      nls=nls, context=pb, verbose=True)

# ---------
# | Solve |
# ---------
print("Solving equation \n\n\t\t u_t - div(au)) = 0\n")
print("With IC: {}".format(ic_fun.name))
# print("and EBCs: {}".format(pb.ebcs.names))
# print("and EPBCS: {}".format(pb.epbcs.names))
print("-------------------------------------")
print("Approximation order is {}".format(approx_order))
print("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
print("CFL coefficient was {0} and order correction {1}".format(CFL, 1 / (2 * approx_order + 1)))
print("Courant number c = max(abs(u)) * dt/dx = {0}".format(max_velo * dtdx))
print("------------------------------------------")
print("Time stepping solver is {}".format(tss.name))
# print("Limiter used: {}".format(limiter.name))
print("======================================")

pb.set_solver(tss)
pb.solve()
