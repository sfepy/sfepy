"""
"""
from __future__ import absolute_import

import numpy as nm

import sfepy.mechanics.matcoefs as mc
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh

# Material parameters.
E = 200e9
nu = 0.3
rho = 7800.0

lam, mu = mc.lame_from_youngpoisson(E, nu)
cl = nm.sqrt((lam + 2.0 * mu) / rho)
cs = nm.sqrt(mu / rho)

# Initial velocity.
v0 = 1

# Mesh dimensions.
d = 2.5e-3
L = 4 * d
dims = [L, 2 * d, 2 * d]

# Mesh resolution: increase to improve accuracy.
shape = [21, 6, 6]

# Element size.
H = L / (shape[0] - 1)

# Time-stepping parameters.
dt = H / cl
t1 = 0.9 * L / cl

def mesh_hook(mesh, mode):
    """
    Generate the block mesh.
    """
    if mode == 'read':
        mesh = gen_block_mesh(dims, shape, [0.5 * L, d, d], name='user_block',
                              verbose=False)
        return mesh

    elif mode == 'write':
        pass

filename_mesh = UserMeshIO(mesh_hook)

regions = {
    'Omega' : 'all',
    'Impact' : ('vertices in (x < 1e-12)', 'facet'),
    'Symmetry-y' : ('vertices in (y < 1e-12)', 'facet'),
    'Symmetry-z' : ('vertices in (z < 1e-12)', 'facet'),
}

# Iron.
materials = {
    'solid' : ({
            'D': mc.stiffness_from_youngpoisson(dim=3, young=E, poisson=nu),
            'rho': rho,
     },),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i' : 1,
}

variables = {
    'u' : ('unknown field', 'displacement', 0, 1),
    'du' : ('unknown field', 'displacement', 1, 1),
    'ddu' : ('unknown field', 'displacement', 2, 1),
    # 'du' : ('parameter field', 'displacement', {'ic' : 'get_ic_du'}, 1),
    # 'ddu' : ('parameter field', 'displacement', {'ic' : 'get_ic_ddu'}, 1),
    'v' : ('test field', 'displacement', 'u'),
    'dv' : ('test field', 'displacement', 'du'),
    'ddv' : ('test field', 'displacement', 'ddu'),
}

ebcs = {
    'Impact' : ('Impact', {'u.all' : 0.0, 'du.all' : 0.0, 'ddu.all' : 0.0}),
    'Symmtery-y' : ('Symmetry-y', {'u.1' : 0.0, 'du.1' : 0.0, 'ddu.1' : 0.0}),
    'Symmetry-z' : ('Symmetry-z', {'u.2' : 0.0, 'du.2' : 0.0, 'ddu.2' : 0.0}),
}

def get_ic(coor, ic, mode='u'):
    val = nm.zeros_like(coor)
    if mode == 'u':
        val[:, 0] = 0.0

    elif mode == 'du':
        val[:, 0] = -1.0

    return val

functions = {
    'get_ic_u' : (get_ic,),
    'get_ic_du' : (lambda coor, ic: get_ic(coor, None, mode='du'),),
}

ics = {
    'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
}

equations = {
    'balance_of_forces' :
    """dw_volume_dot.i.Omega(solid.rho, ddv, ddu)
     + dw_zero.i.Omega(dv, du)
     + dw_lin_elastic.i.Omega(solid.D, v, u) = 0""",
    # 'balance_of_forces' :
    # """dw_volume_dot.i.Omega(solid.rho, v, d2u/dt2)
    #  + dw_lin_elastic.i.Omega(solid.D, v, u) = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'ls-i' : ('ls.petsc', {
        'method' : 'cg',
        'precond' : 'icc',
        'i_max' : 150,
        'eps_a' : 1e-32,
        'eps_r' : 1e-8,
        'verbose' : 2,
    }),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-6,
    }),
    # 'tss' : ('ts.simple', {
    #     't0' : 0.0,
    #     't1' : 1.0,
    #     'dt' : 0.1,
    #     'n_step' : None,
    # }),
    # 'tsb' : ('ts.bathe', {
    #     't0' : 0.0,
    #     't1' : 1.0,
    #     'dt' : 0.1,
    #     'n_step' : None,
    # }),
    'tsn' : ('ts.newmark', {
        't0' : 0.0,
        't1' : t1,
        'dt' : dt,
        'n_step' : None,

        'u' : 'u',
        'v' : 'du',
        'a' : 'ddu',

        'verbose' : 1,
    }),
    # 'tsp' : ('ts.petsc', {
    #     'method' : 'rk',
    #     't0' : 0.0,
    #     't1' : 1.0,
    #     'dt' : 0.1,
    #     'n_step' : None,
    # }),
}

options = {
    'ts' : 'ts',
    'nls' : 'newton',
    'ls' : 'ls-i',

    'output_format' : 'h5',
    'output_dir' : 'output/ed',
}
