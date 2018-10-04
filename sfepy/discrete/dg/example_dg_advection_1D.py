import numpy as nm
from my_utils.inits_consts import left_par_q


# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

# local import
from dg_terms import AdvFluxDGTerm, AdvIntDGTerm
from dg_equation import Equation
from dg_tssolver import TSSolver

n_nod = 101
coors = nm.linspace(0.0, 1.0, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_1d', coors, None,
                      [conn], [mat_ids], descs)

a = 1.2
IntT = AdvIntDGTerm(mesh)
FluxT = AdvFluxDGTerm(mesh, a)
eq = Equation((IntT, FluxT))
ic = left_par_q((mesh.coors[:-1] + mesh.coors[1:])/2)
bc = {"left" : 0,
      "right" : 0}

tss = TSSolver(eq, ic, bc)
tss.solve(0, 1, 10)
