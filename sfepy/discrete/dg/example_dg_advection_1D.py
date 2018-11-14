import numpy as nm
import matplotlib.pyplot as plt

# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions
from sfepy.terms.terms import Term
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver


# local import
from dg_terms import AdvFluxDGTerm, AdvIntDGTerm
# from dg_equation import Equation
from dg_tssolver import TSSolver, RK3Solver, EUSolver
from dg_field import DGField

from my_utils.inits_consts import left_par_q, gsmooth, const_u, ghump, superic
from my_utils.visualizer import animate1d, sol_frame

X1 = 0.
XN1 = 1.
n_nod = 100
n_el = n_nod - 1
coors = nm.linspace(X1, XN1, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_1d', coors, None,
                      [conn], [mat_ids], descs)

ts = 0
te = 2
tn = 800

domain = FEDomain('domain', mesh)  # TODO DGDomain
omega = domain.create_region('Omega', 'all')
left = domain.create_region('Gamma1',
                              'vertices in x == %.10f' % X1,
                              'vertex')
right = domain.create_region('Gamma2',
                              'vertices in x == %.10f' % XN1,
                              'vertex')
field = DGField.from_args('fu', nm.float64, 'vector', omega,
                        approx_order=2)
u = FieldVariable('u', 'unknown', field, history=1)
v = FieldVariable('v', 'test', field, primary_var_name='u')
integral = Integral('i', order=2)

# TODO use sfepy volume term?

f = Material('f', val=[1.0])  # TODO how do materials really work?
IntT = Term.new("dw_volume_lvf(f.val, v)", integral, omega, v=v, f=f)

a = Material('a', val=[1.0])
FluxT = AdvFluxDGTerm(integral, omega, u=u, v=v, a=a)


eq = Equation('balance', IntT + FluxT)
eqs = Equations([eq])

left_fix_u = EssentialBC('left_fix_u', left, {'u.all' : 0.0})
right_fix_u = EssentialBC('right_fix_u', right, {'u.all' : 0.0})

ic_fun = Function('ic_fun', lambda x, ic: superic(x))  # why does IC function get ic object?
ics = InitialCondition('ic', omega, {'u.0': ic_fun})  # TODO how to initialize variable with IC?

pb = Problem('advection', equations=eqs)
pb.set_bcs(ebcs=Conditions([left_fix_u, right_fix_u]))
pb.set_ics(Conditions([ics]))

state0 = pb.get_initial_state()
# it kinda works up until now


geometry = Struct(n_vertex=2,
                  dim=1,
                  coors=coors.copy())

ic = superic
bc = {"right" : 0.0,
      "left" : 0.0}

# TODO use sfepy solver, which one? How do we get complete solution out?
# tss = EUSolver(eq, ic, bc, TSSolver.moment_limiter, LegendrePolySpace("legb", geometry, 1))

# from time_poisson_interactive.py
ls = ScipyDirect({})
nls_status = IndexedStruct()
nls = Newton({'is_linear' : True}, lin_solver=ls, status=nls_status)
tss = SimpleTimeSteppingSolver({'t0' : 0.0, 't1' : 1.0, 'n_step' : 100},
                               nls=nls, context=pb, verbose=True)
pb.set_solver(tss)

pb.time_update(tss.ts)
state0.apply_ebc()

tss_status = IndexedStruct()
pb.solve()
# tss(state0.get_vec(pb.active_only),
#     status=tss_status)

print(tss_status)
# u, dt = tss.solve(ts, te, tn)
# sic = tss.initial_cond


#--------
#| Plot |
#--------
# TODO move plotting to visualizer
plt.figure("Sampled Solution anim")
X = (mesh.coors[1:] + mesh.coors[:-1])/2
T = nm.linspace(ts, te, tn)
# sic = TSSolver.initial_cond

# Plot mesh
plt.vlines(mesh.coors[:, 0], ymin=0, ymax=.5, colors="grey")
plt.vlines((mesh.coors[0], mesh.coors[-1]), ymin=0, ymax=.5, colors="k")
plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")

# Plot IC and its sampling
# TODO get IC sampling, from where?
# c0 = plt.plot(X, sic[0, :, 0], label="IC-0", marker=".", ls="")[0].get_color()
# c1 = plt.plot(X, sic[1, :, 0], label="IC-1", marker=".", ls="")[0].get_color()
# # plt.plot(coors, .1*alones(n_nod), marker=".", ls="")
# plt.step(coors[1:], sic[0, :, 0], label="IC-0", color=c0)
# plt.step(coors[1:], sic[1, :, 0], label="IC-1", color=c1)
# plt.plot(coors[1:], sic[1, :], label="IC-1", color=c1)
xs = nm.linspace(X1, XN1, 500)[:, None]
plt.plot(xs, ic(xs), label="IC-ex")

# Animate sampled solution
anim = animate1d(u[:, :, :, 0].T, nm.append(coors, coors[-1]), T, ylims=[-1, 2], plott="step")
plt.xlim(coors[0]-.1, coors[-1]+.1)
plt.legend(loc="upper left")
plt.title("Sampled solution")

plt.figure("Reconstructed Solution anim")
# Plot mesh
plt.vlines(mesh.coors[:, 0], ymin=0, ymax=.5, colors="grey")
plt.vlines((mesh.coors[0], mesh.coors[-1]), ymin=0, ymax=.5, colors="k")
plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")

# Prepare reconstructed solution
# ww = nm.zeros((2*n_nod, tn, 1))
# ww[0, :] = u[0, 0, :] - u[1, 0, :]
# ww[-1, :] = u[0, -1, :] + u[1, -1, :]
# ww[:-2:2] = u[0, 1:-1, :] - u[1, 1:-1, :]
# ww[1:-1:2] = u[0, 1:-1, :] + u[1, 1:-1, :]
#
# # nodes for plotting reconstructed solution
# xx = nm.zeros((2*n_nod, 1))
# xx[0] = mesh.coors[0]
# xx[-1] = mesh.coors[-1]
# xx[2::2] = mesh.coors[1:]
# xx[1:-1:2] = mesh.coors[1:]
# # plt.vlines(xx, ymin=0, ymax=.3, colors="green")
#

# Plot discontinuously!
ww = nm.zeros((3*n_nod-1, tn, 1))
ww[0, :] = u[0, 0, :] - u[1, 0, :]  # left bc
ww[-1, :] = u[0, -1, :] + u[1, -1, :]  # right bc

ww[0:-2:3] = u[0, 1:-1, :] - u[1, 1:-1, :]  # left edges of elements
ww[1:-1:3] = u[0, 1:-1, :] + u[1, 1:-1, :]  # right edges of elements
ww[2::3, :] = nm.NaN  # NaNs ensure plotting of discontinuities at element borders

# nodes for plotting reconstructed solution
xx = nm.zeros((3*n_nod-1, 1))
xx[0] = mesh.coors[0]
xx[-1] = mesh.coors[-1]
# the ending ones are still a bit odd, but hey, it works!
xx[1:-1] = nm.repeat(mesh.coors[1:], 3)[:, None]
# plt.vlines(xx, ymin=0, ymax=.3, colors="green")

# plot reconstructed IC
plt.plot(xx, ww[:, 0], label="IC")

# Animate reconstructed
anim_disc = animate1d(ww[:, :, 0].T, xx, T, ylims=[-1, 2])
plt.xlim(coors[0]-.1, coors[-1]+.1)
plt.legend(loc="upper left")
plt.title("Reconstructed solution")

# sol_frame(u[:, :, :, 0].T, nm.append(coors, coors[-1]), T, t0=0., ylims=[-1, 1], plott="step")

plt.show()
