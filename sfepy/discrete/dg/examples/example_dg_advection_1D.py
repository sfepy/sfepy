import numpy as nm
import matplotlib.pyplot as plt

# sfepy imports
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver

from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm
from sfepy.discrete.fem.meshio import VTKMeshIO
from sfepy.base.ioutils import ensure_path

# local imports
from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm
from sfepy.discrete.dg.dg_tssolver import TVDRK3StepSolver, RK4StepSolver, EulerStepSolver
from sfepy.discrete.dg.dg_field import DGField
from sfepy.discrete.dg.dg_limiters import IdentityLimiter, Moment1DLimiter


from sfepy.discrete.dg.my_utils.inits_consts import \
    left_par_q, gsmooth, const_u, ghump, superic
from sfepy.discrete.dg.my_utils.visualizer import load_1D_vtks, plot1D_DG_sol
from sfepy.discrete.dg.my_utils.read_plot_1Ddata import clear_output_folder

# Setup output names
domain_name = "domain_1D"
output_folder = "output/adv_1D/"
output_folder_mesh = "output/mesh"
save_timestn = 100
clear_output_folder(output_folder)

#------------
#| Get mesh |
#------------
X1 = 0.
XN = 1.
n_nod = 100
n_el = n_nod - 1
coors = nm.linspace(X1, XN, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_1d', coors, None,
                      [conn], [mat_ids], descs)

outfile = "output/mesh/uniform_1D_mesh.vtk"
ensure_path(outfile)
meshio = VTKMeshIO(outfile)
meshio.write(outfile, mesh)


#-----------------------------
#| Create problem components |
#-----------------------------
#vvvvvvvvvvvvvvvv#
approx_order = 1
#^^^^^^^^^^^^^^^^#
integral = Integral('i', order=approx_order * 2)
domain = FEDomain('domain_1D', mesh)
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
# TODO u and v are flipped what now?!
StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, v, u)", "a.val, u[-1], v",
                                  integral, omega, u=u, v=v, a=a)

alpha = Material('alpha', val=[.0])
FluxT = AdvectDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
                         integral, omega, u=u, v=v, a=a, alpha=alpha)

eq = Equation('balance', MassT + StiffT - FluxT)
eqs = Equations([eq])


#------------------------------
#| Create bounrady conditions |
#------------------------------
left_fix_u = EssentialBC('left_fix_u', left, {'u.all' : 3.2})#, times="all")
right_fix_u = EssentialBC('right_fix_u', right, {'u.all' : 4.2})#, times="all")

#----------------------------
#| Create initial condition |
#----------------------------
def ic_wrap(x, ic=None):
    return ghump(x - .3)

ic_fun = Function('ic_fun', ic_wrap)
ics = InitialCondition('ic', omega, {'u.0': ic_fun})


#------------------
#| Create problem |
#------------------
pb = Problem('advection', equations=eqs, conf=Struct(options={"save_times": save_timestn}, ics={},
                                                     ebcs={}, epbcs={}, lcbcs={}, materials={}))
pb.setup_output(output_dir="output/adv_1D")  # , output_format="msh")
# pb.set_bcs(ebcs=Conditions([left_fix_u, right_fix_u]))
pb.set_ics(Conditions([ics]))

state0 = pb.get_initial_state()
pb.save_state("output/adv_1D/domain_1D_start.vtk", state=state0)


#------------------
#| Create limiter |
#------------------
limiter = Moment1DLimiter(field.n_el_nod, field.n_cell)


#---------------------------
#| Set time discretization |
#---------------------------
CFL = .4
max_velo = nm.max(nm.abs(velo))
t0 = 0
t1 = .2
dx = (XN - X1) / n_nod
dt = dx / nm.abs(velo) * CFL/(2*approx_order + 1)
tn = int(nm.ceil((t1 - t0) / dt))
dtdx = dt / dx
#------------------
#| Create solver |
#------------------
ls = ScipyDirect({})
nls_status = IndexedStruct()
nls = Newton({'is_linear': True}, lin_solver=ls, status=nls_status)

tss = EulerStepSolver({'t0': t0, 't1': t1, 'n_step': tn},
                         nls=nls, context=pb, verbose=True)
                        # ,post_stage_hook=limiter)
#
# tss = TVDRK3StepSolver({'t0': t0, 't1': t1, 'n_step': tn},
#                          nls=nls, context=pb, verbose=True)

# tss = RK4StepSolver({'t0': t0, 't1': t1, 'n_step': tn},
#                          nls=nls, context=pb, verbose=True, post_stage_hook=limiter)


#---------
#| Solve |
#---------
print("Solving equation \n\n\t\t u_t - div(au)) = 0\n")
print("With IC: {}".format(ic_fun.name))
# print("and EBCs: {}".format(pb.ebcs.names))
# print("and EPBCS: {}".format(pb.epbcs.names))
print("-------------------------------------")
print("Approximation order is {}".format(approx_order))
print("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
print("CFL coefficient was {0} and order correction {1}".format(CFL, 1/(2*approx_order + 1)))
print("Courant number c = max(abs(u)) * dt/dx = {0}".format(max_velo * dtdx))
print("------------------------------------------")
print("Time stepping solver is {}".format(tss.name))
print("Limiter used: {}".format(limiter.name))
print("======================================")

pb.set_solver(tss)
state_end = pb.solve()
pb.save_state("output/adv_1D/domain_1D_end.vtk", state=state_end)


#----------
#| Plot 1D|
#----------
# lmesh, u = load_1D_vtks("./output/adv_1D", "domain_1D", order=approx_order)
# plot1D_DG_sol(lmesh, t0, t1, u, tn=30, ic=ic_wrap,
#               delay=100, polar=False)
#
# from sfepy.discrete.dg.dg_field import get_unraveler, get_raveler
# from sfepy.discrete.dg.my_utils.visualizer import \
#     load_state_1D_vtk, plot_1D_legendre_dofs, reconstruct_legendre_dofs
# coors, u_end = load_state_1D_vtk("output/adv_1D/domain_1D_end.vtk", order=approx_order)
#
#
# u_start = get_unraveler(field.n_el_nod, field.n_cell)(state0.vec).swapaxes(0, 1)[..., None]
# # u_end = get_unraveler(field.n_el_nod, field.n_cell)(state_end.vec).swapaxes(0, 1)[..., None]
#
#
# plot_1D_legendre_dofs(coors, [u_start.swapaxes(0, 1)[:, :, 0], u_end.swapaxes(0, 1)[:, :, 0]])
#
# plt.figure("reconstructed")
# ww_s, xx = reconstruct_legendre_dofs(coors, None, u_end)
# ww_e, _ = reconstruct_legendre_dofs(coors, None, u_start)
#
# plt.plot(xx, ww_s[:, 0])
# plt.plot(xx, ww_e[:, 0])
# plt.show()
from sfepy.discrete.dg.my_utils.read_plot_1Ddata import load_and_plot_fun

load_and_plot_fun(output_folder, domain_name, t0, t1, min(tn, save_timestn), approx_order, ic_fun)
