import numpy as nm
from  numpy.linalg import norm
import matplotlib.pyplot as plt


# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions
from sfepy.terms.terms import Term
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import VTKMeshIO
from sfepy.base.ioutils import ensure_path
from sfepy.base.conf import ProblemConf
from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm


# local import
from sfepy.discrete.dg.dg_terms import AdvFluxDGTerm
# from dg_equation import Equation
from sfepy.discrete.dg.dg_tssolver \
    import EulerStepSolver, TVDRK3StepSolver
from sfepy.discrete.dg.dg_field import DGField

from sfepy.discrete.dg.my_utils.inits_consts \
    import left_par_q, gsmooth, const_u, ghump, superic

# Setup output names
domain_name = "domain_tensor_12D"
output_folder = "output/adv_tens_12D"
output_folder_mesh = "output/mesh"
save_timesn = 100


#------------
#| Get mesh |
#------------
mesh = gen_block_mesh((1., 1.), (20, 2), (.5, 0.5))
outfile = "output/mesh/tens_12D_mesh.vtk"
ensure_path(outfile)
meshio = VTKMeshIO(outfile)
meshio.write(outfile, mesh)


#-----------------------------
#| Create problem components |
#-----------------------------
#vvvvvvvvvvvvvvvv#
approx_order = 1
#^^^^^^^^^^^^^^^^#
integral = Integral('i', order=5)
domain = FEDomain(domain_name, mesh)
omega = domain.create_region('Omega', 'all')
dgfield = DGField('dgfu', nm.float64, 'scalar', omega,
                  approx_order=approx_order)

u = FieldVariable('u', 'unknown', dgfield, history=1)
v = FieldVariable('v', 'test', dgfield, primary_var_name='u')

MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u", integral, omega, u=u, v=v)

velo = nm.array([[1., 0.]]).T
a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, u, v)", "a.val, u[-1], v", integral, omega,
                                    u=u, v=v, a=a, mode="grad_virtual")

FluxT = AdvFluxDGTerm(integral, omega, u=u, v=v, a=a)

eq = Equation('balance', MassT + StiffT - FluxT)
eqs = Equations([eq])


#------------------------------
#| Create bounrady conditions |
#------------------------------
# TODO BCs


def ic_wrap(x, ic=None):
    return gsmooth(x[..., 0:1])
#-----------
#| Plot IC |
#-----------
# X = nm.arange(0, 1, 0.005)
# Y = nm.arange(0, 1, 0.005)
# X, Y = nm.meshgrid(X, Y)
# coors = nm.dstack((X[:, :, None], Y[:, :, None]))
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# from matplotlib import cm
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# Z = ic_wrap(coors)
# surf = ax.plot_surface(X, Y, Z[:, :, 0], cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False, alpha=.6)
#
# fig = plt.figure()
# ax = fig.gca()
# ax.contour(X, Y, Z[..., 0])
# plt.show()
ic_fun = Function('ic_fun', ic_wrap)
ics = InitialCondition('ic', omega, {'u.0': ic_fun})

#------------------
#| Create problem |
#------------------
pb = Problem('advection', equations=eqs, conf=Struct(options={"save_times": "all"}, ics={},
                                                     ebcs={}, epbcs={}, lcbcs={}, materials={}))
pb.setup_output(output_dir=output_folder, output_format="msh")
# pb.set_bcs(ebcs=Conditions([left_fix_u, right_fix_u]))
pb.set_ics(Conditions([ics]))

state0 = pb.get_initial_state()
pb.save_state("output/state0_tens_12D.msh", state=state0)


#---------------------------
#| Set time discretization |
#---------------------------
CFL = .4
max_velo = nm.max(nm.abs(velo))
t0 = 0
t1 = .2
dx = nm.min(mesh.cmesh.get_volumes(2))
dt = dx / norm(velo) * CFL/(2*approx_order + 1)
# time_steps_N = int((tf - t0) / dt) * 2
tn = int(nm.ceil((t1 - t0) / dt))
dtdx = dt / dx
#------------------
#| Create solver |
#------------------
ls = ScipyDirect({})
nls_status = IndexedStruct()
nls = Newton({}, lin_solver=ls, status=nls_status)

# tss = EulerStepSolver({'t0': t0, 't1': t1, 'n_step': tn},
#                         nls=nls, context=pb, verbose=True)

tss = TVDRK3StepSolver({'t0': t0, 't1': t1, 'n_step': tn},
                         nls=nls, context=pb, verbose=True)


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
# print("Limiter used: {}".format(limiter.name))
print("======================================")

pb.set_solver(tss)
state_end = pb.solve()
# pb.save_state(output_folder, state=state_end)
