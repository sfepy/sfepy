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
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/rectangle_fine_quad.mesh'

regions = {
    'Omega' : 'all',
    'Surface' : ('vertices of surface', 'facet'),
    'Right' : ('vertices in (x > 4.999)', 'facet'),
    'Bottom' : ('vertices in (y < -9.999)', 'facet'),
    'Top' : ('vertices in (y > 9.999)', 'facet'),
    'Walls' : ('r.Top +v r.Right +v r.Bottom', 'facet'),
    'Driven' : ('r.Surface -v r.Walls', 'facet'),
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
    'Walls' : ('Walls', {'u.all' : 0.0}),
    'Driven' : ('Driven', {'u.1' : 0.1, 'u.0' : 0.0}),
}

equations = {
    'balance' :
    """+ dw_div_grad.5.Omega(fluid.viscosity, v, u)
       + dw_convect.5.Omega(v, u)
       - dw_stokes.5.Omega(v, p) = 0""",

    'incompressibility' :
    """dw_stokes.5.Omega(u, q) = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 15,
        'eps_a'      : 1e-10,
        'eps_r'      : 1.0,
        'problem'   : 'nonlinear'
    }),
}
