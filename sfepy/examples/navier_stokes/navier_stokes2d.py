# -*- coding: utf-8 -*-
r"""
Navier-Stokes equations for incompressible fluid flow in 2D.

Find :math:`\ul{u}`, :math:`p` such that:

.. math::
    \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}
    + \int_{\Omega} ((\ul{u} \cdot \nabla) \ul{u}) \cdot \ul{v}
    - \int_{\Omega} p\ \nabla \cdot \ul{v}
    = 0
    \;, \quad \forall \ul{v} \;,

    \int_{\Omega} q\ \nabla \cdot \ul{u}
    = 0
    \;, \quad \forall q \;.

The mesh is created by ``gen_block_mesh()`` function.

View the results using::

  sfepy-view user_block.vtk -2
"""
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh

# Mesh dimensions.
dims = [0.1, 0.1]

# Mesh resolution: increase to improve accuracy.
shape = [51, 51]

def mesh_hook(mesh, mode):
    """
    Generate the block mesh.
    """
    if mode == 'read':
        mesh = gen_block_mesh(dims, shape, [0, 0], name='user_block',
                              verbose=False)
        return mesh

    elif mode == 'write':
        pass

filename_mesh = UserMeshIO(mesh_hook)

regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < -0.0499)', 'facet'),
    'Right' : ('vertices in (x > 0.0499)', 'facet'),
    'Bottom' : ('vertices in (y < -0.0499)', 'facet'),
    'Top' : ('vertices in (y > 0.0499)', 'facet'),
    'Walls' : ('r.Left +v r.Right +v r.Bottom', 'facet'),
}

materials = {
    'fluid' : ({'viscosity' : 1.00e-2},),
}

fields = {
    'velocity': ('real', 'vector', 'Omega', 2),
    'pressure': ('real', 'scalar', 'Omega', 1),
}

variables = {
    'u' : ('unknown field', 'velocity', 0),
    'v' : ('test field', 'velocity', 'u'),
    'p' : ('unknown field', 'pressure', 1),
    'q' : ('test field', 'pressure', 'p'),
}

ebcs = {
    '1_Walls' : ('Walls', {'u.all' : 0.0}),
    '0_Driven' : ('Top', {'u.0' : 1.0, 'u.1' : 0.0}),
    'Pressure' : ('Bottom', {'p.0' : 0.0}),
}

integrals = {
    'i' : 4,
}

equations = {
    'balance' :
    """+ dw_div_grad.i.Omega(fluid.viscosity, v, u)
       + dw_convect.i.Omega(v, u)
       - dw_stokes.i.Omega(v, p) = 0""",

    'incompressibility' :
    """dw_stokes.i.Omega(u, q) = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 15,
        'eps_a'      : 1e-10,
        'eps_r'      : 1.0,
    }),
}
