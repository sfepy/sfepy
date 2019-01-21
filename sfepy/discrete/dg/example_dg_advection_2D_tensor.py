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
from sfepy.base.conf import ProblemConf
from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm

# local import
from dg_terms import AdvFluxDGTerm, AdvVolDGTerm, ScalarDotMGradScalarDGTerm
# from dg_equation import Equation
from dg_tssolver import EulerStepSolver, DGTimeSteppingSolver, RK3StepSolver
from dg_field import DGField

from my_utils.inits_consts import left_par_q, gsmooth, const_u, ghump, superic
from my_utils.visualizer import load_vtks, plot1D_DG_sol

mesh = gen_block_mesh((1., 1.), (100, 100), (0.5, 0.5))
outfile = "output/mesh/tensor_2D_mesh.vtk"
meshio = VTKMeshIO(outfile)
meshio.write(outfile, mesh)


domain = FEDomain('domain_tensor_2D', mesh)
omega = domain.create_region('Omega', 'all')
#vvvvvvvvvvvvvvvv#
approx_order = 2
#^^^^^^^^^^^^^^^^#
field = DGField('dgfu', nm.float64, 'scalar', omega,
                approx_order=approx_order)

u = FieldVariable('u', 'unknown', field, history=1)
v = FieldVariable('v', 'test', field, primary_var_name='u')

velo = nm.array([[1., 0.]]).T

t0 = 0
t1 = 1
max_velo = nm.max(nm.abs(velo))

dx = nm.min(mesh.cmesh.get_volumes(2))
dt = dx / norm(velo) * 1/2
# time_steps_N = int((tf - t0) / dt) * 2
tn = int(nm.ceil((t1 - t0) / dt))
dtdx = dt / dx
print("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
print("Courant number c = max(abs(u)) * dt/dx = {0}".format(max_velo * dtdx))

integral = Integral('i', order=5)
IntT = AdvVolDGTerm(integral, omega, u=u, v=v)

a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarDGTerm("adv_stiff(a.val, u, v)", "a.val, u, v", integral, omega,
                                    u=u, v=v, a=a, mode="grad_virtual")

FluxT = AdvFluxDGTerm(integral, omega, u=u, v=v, a=a)

eq = Equation('balance', IntT + StiffT + FluxT)
eqs = Equations([eq])


def ic_wrap(x, ic=None):
    return superic(x[..., 0:1])*gsmooth(x[..., 1:])


X = nm.arange(0, 1, 0.005)
Y = nm.arange(0, 1, 0.005)
X, Y = nm.meshgrid(X, Y)
coors = nm.dstack((X[:, :, None], Y[:, :, None]))

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

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

pb = Problem('advection', equations=eqs)
pb.setup_output(output_dir="./output/", output_format="h5")
# pb.set_bcs(ebcs=Conditions([left_fix_u, right_fix_u]))
pb.set_ics(Conditions([ics]))

state0 = pb.get_initial_state()
pb.save_state("output/state0_tensor_2D.msh", state=state0)

ls = ScipyDirect({})
nls_status = IndexedStruct()
nls = RK3StepSolver({}, lin_solver=ls, status=nls_status)

tss = DGTimeSteppingSolver({'t0': t0, 't1': t1, 'n_step': tn},
                                nls=nls, context=pb, verbose=True)
pb.set_solver(tss)
pb.solve()