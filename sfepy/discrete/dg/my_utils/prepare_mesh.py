"""
Script to quickly prepare base meshes for convergence testing
"""

import numpy as nm

from discrete.fem import Mesh

X1 = 0.
XN = 2*nm.pi
n_nod = 2
n_el = n_nod - 1
coors = nm.linspace(X1, XN, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('uniform_1D{}'.format(n_nod), coors, None,
                      [conn], [mat_ids], descs)
mesh.write("meshes/1d/mesh_tensr_1D_02pi_2.vtk")
