import numpy as nm
from numpy.linalg import norm
import matplotlib.pyplot as plt
from os.path import join as pjoin

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
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import VTKMeshIO
from sfepy.base.ioutils import ensure_path

from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm

# local imports
from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm
from sfepy.discrete.dg.dg_tssolver \
    import EulerStepSolver, TVDRK3StepSolver
from sfepy.discrete.dg.dg_field import DGField
from sfepy.discrete.dg.dg_limiters import IdentityLimiter, Moment1DLimiter

from sfepy.discrete.dg.my_utils.inits_consts \
    import left_par_q, gsmooth, const_u, ghump, superic

from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder


#vvvvvvvvvvvvvvvv#
approx_order = 1
#^^^^^^^^^^^^^^^^#
# Setup  names
domain_name = "domain_2D"
problem_name = "iadv_2D_simp"
output_folder = pjoin("output", problem_name, str(approx_order))
output_format = "msh"
mesh_output_folder = "output/mesh"
save_timestn = 100

#------------
#| Get mesh |
#------------
# mesh = gen_block_mesh((1., 1.), (20, 20), (0.5, 0.5))
# mesh = triangulate(mesh)
mesh_name = "square_tri2"
mesh = Mesh.from_file("mesh/" + mesh_name + ".mesh")


#-----------------------------
#| Create problem components |
#-----------------------------

integral = Integral('i', order=approx_order * 2)
domain = FEDomain(domain_name, mesh)
omega = domain.create_region('Omega', 'all')
field = DGField('dgfu', nm.float64, 'scalar', omega,
                  approx_order=approx_order)

u = FieldVariable('u', 'unknown', field, history=1)
v = FieldVariable('v', 'test', field, primary_var_name='u')


MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u", integral, omega, u=u, v=v)

velo = nm.array([[1., 0.]]).T
a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, u, v)", "a.val, u[-1], v", integral, omega,
                                    u=u, v=v, a=a)

alpha = Material('alpha', val=[.0])
FluxT = AdvectDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
                         integral, omega, u=u, v=v, a=a, alpha=alpha)

eq = Equation('balance', MassT + StiffT - FluxT)
eqs = Equations([eq])


#------------------------------
#| Create bounrady conditions |
#------------------------------
# TODO BCs

#----------------------------
#| Create initial condition |
#----------------------------
def ic_wrap(x, ic=None):
    x = (x + 1.3) / 4
    return  gsmooth(x[..., 0:1]) * gsmooth(x[..., 1:])
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
pb = Problem(problem_name, equations=eqs, conf=Struct(options={"save_times": save_timestn}, ics={},
                                                     ebcs={}, epbcs={}, lcbcs={}, materials={}),
             active_only=False)
pb.setup_output(output_dir=output_folder, output_format=output_format)
pb.set_ics(Conditions([ics]))


#------------------
#| Create limiter |
#------------------
limiter = IdentityLimiter

#---------------------------
#| Set time discretization |
#---------------------------
CFL = .4
max_velo = nm.max(nm.linalg.norm(velo))
t0 = 0
t1 = .2
dx = nm.min(mesh.cmesh.get_volumes(2))
dt = dx / max_velo * CFL/(2*approx_order + 1)
tn = int(nm.ceil((t1 - t0) / dt))
dtdx = dt / dx

#------------------
#| Create solver |
#------------------
ls = ScipyDirect({})
nls_status = IndexedStruct()
nls = Newton({'is_linear': True}, lin_solver=ls, status=nls_status)

tss_conf = {'t0': t0,
            't1': t1,
            'n_step': tn,
            "limiter": limiter}

tss = EulerStepSolver(tss_conf,
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
