import numpy as nm
import matplotlib.pyplot as plt


# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct


# local import
from dg_terms import AdvFluxDGTerm, AdvIntDGTerm
from dg_equation import Equation
from dg_tssolver import TSSolver, RK3Solver
from dg_basis import LegendrePolySpace

from my_utils.inits_consts import left_par_q, gsmooth, const_u, ghump, superic
from my_utils.visualizer import animate1d

X1 = -3.
XN1 = 7.
n_nod = 100
n_el = n_nod - 1
coors = nm.linspace(X1, XN1, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_1d', coors, None,
                      [conn], [mat_ids], descs)

a = 1.0
ts = 0
te = 10
tn = 1000

IntT = AdvIntDGTerm(mesh)
FluxT = AdvFluxDGTerm(mesh, a)

eq = Equation((IntT, FluxT))

ic = gsmooth
bc = {"left" : 0,
      "right" : 0}

geometry = Struct(n_vertex=2,
                  dim=1,
                  coors=coors.copy())

tss = RK3Solver(eq, ic, bc, LegendrePolySpace("legb", geometry, 1))

u, dt = tss.solve(ts, te, tn)
sic = tss.initial_cond
X = (mesh.coors[1:] + mesh.coors[:-1])/2
T = nm.linspace(ts, te, tn)


plt.figure("Numeric Solution anim")
# Plot mesh
plt.vlines(mesh.coors[:, 0], ymin=0, ymax=.5, colors="grey")
plt.vlines((mesh.coors[0], mesh.coors[-1]), ymin=0, ymax=.5, colors="k")
plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")


# Plot IC
plt.plot(X, sic[0, :, 0], label="IC-0", marker=".", ls="")
plt.plot(X, sic[1, :, 0], label="IC-1", marker=".", ls="")
xs = nm.linspace(X1, XN1, 500)[:, None]
plt.plot(xs, gsmooth(xs), label="IC-ex")

# Animate solution
anim = animate1d(u[:, :, :, 0].T, nm.append(nm.append(coors[0], X), coors[-1]), T, ylims=[-1, 1])
plt.xlim(coors[0]-.1, coors[-1]+.1)
plt.legend(loc="upper left")
plt.title("Numeric solution")

plt.show()
