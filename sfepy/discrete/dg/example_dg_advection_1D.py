import numpy as nm
import matplotlib.pyplot as plt


# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

# local import
from dg_terms import AdvFluxDGTerm, AdvIntDGTerm
from dg_equation import Equation
from dg_tssolver import TSSolver, RK3Solver
from dg_basis import DGBasis

from my_utils.inits_consts import left_par_q
from my_utils.visualizer import animate1d

n_nod = 100
n_el = n_nod - 1
coors = nm.linspace(0.0, 1.0, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_1d', coors, None,
                      [conn], [mat_ids], descs)

a = 1
ts = 0
te = 1
tn = 300

IntT = AdvIntDGTerm(mesh)
FluxT = AdvFluxDGTerm(mesh, a)

eq = Equation((IntT, FluxT))

ic = left_par_q
bc = {"left" : 0,
      "right" : 0}

tss = RK3Solver(eq, ic, bc, DGBasis(1))

u, dt = tss.solve(ts, te, tn)
sic = tss.initial_cond
X = (mesh.coors[1:] + mesh.coors[:-1])/2
T = nm.linspace(ts, te, tn)


plt.figure("Numeric Solution anim")
plt.vlines(mesh.coors[:, 0], ymin=0, ymax=.5)
plt.vlines(X, ymin=0, ymax=.3)

plt.plot(X, sic[0, :, 0])
plt.plot(X, sic[1, :, 0])
anim = animate1d(u[:, 1:-1, :, 0].T, X, T, ylims=[-2, 2])
plt.legend(loc="upper left")
plt.title("Numeric solution")

plt.show()