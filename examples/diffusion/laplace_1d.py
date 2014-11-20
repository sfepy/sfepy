r"""
Laplace equation in 1D.
"""
import numpy as nm
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

def mesh_hook(mesh, mode):
    """
    Generate the 1D mesh from data.
    """
    if mode == 'read':
        n_nod = 11

        coors = nm.linspace(0.0, 1.0, n_nod).reshape((n_nod, 1))
        conns = nm.array([[i, i+1] for i in range(n_nod - 1)], dtype=nm.int32)
        mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
        descs = ['1_2'] * n_nod

        mesh = Mesh.from_data('laplace_1d', coors, None,
                              [conns], [mat_ids], descs)
        return mesh

    elif mode == 'write':
        pass

filename_mesh = UserMeshIO(mesh_hook)

materials = {
    'coef' : ({'val' : 1.0},),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.99999)', 'facet'),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integrals = {
    'i' : 2,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega(coef.val, s, t) = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}
