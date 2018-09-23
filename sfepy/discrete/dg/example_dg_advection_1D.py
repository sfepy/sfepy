import numpy as nm


# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

# local import
from dg_terms import AdvFluxTerm, AdvIntTerm
from dg_equation import Equation

n_nod = 11
coors = nm.linspace(0.0, 1.0, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_1d', coors, None,
                      [conn], [mat_ids], descs)

a = 1.2
IntT = AdvIntTerm(mesh, a)
FluxT = AdvFluxTerm(mesh, a)
eq = Equation((IntT, FluxT))
A = eq.assemble()
print A
