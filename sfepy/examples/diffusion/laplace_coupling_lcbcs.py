r"""
Two Laplace equations with multiple linear combination constraints.

The two equations are coupled by a periodic-like boundary condition constraint
with a shift, given as a non-homogeneous linear combination boundary condition.

Find :math:`u` such that:

.. math::
    \int_{\Omega_1} \nabla v_1 \cdot \nabla u_1
    = 0
    \;, \quad \forall v_1 \;,

    \int_{\Omega_2} \nabla v_2 \cdot \nabla u_2
    = 0
    \;, \quad \forall v_2 \;,

    u_1 = 0 \mbox{ on } \Gamma_{bottom} \;,

    u_2 = 1 \mbox{ on } \Gamma_{top} \;,

    u_1(\ul{x}) = u_2(\ul{x}) + a(\ul{x}) \mbox{ for }
    \ul{x} \in \Gamma = \bar\Omega_1 \cap \bar\Omega_2

    u_1(\ul{x}) = u_1(\ul{y}) + b(\ul{y}) \mbox{ for }
    \ul{x} \in \Gamma_{left}, \ul{y} \in \Gamma_{right}, \ul{y} = P(\ul{x}) \;,

    u_1 = c_{11} \mbox{ in } \Omega_{m11} \subset \Omega_1 \;,

    u_1 = c_{12} \mbox{ in } \Omega_{m12} \subset \Omega_1 \;,

    u_2 = c_2 \mbox{ in } \Omega_{m2} \subset \Omega_2 \;,

where :math:`a(\ul{x})`, :math:`b(\ul{y})` are given functions (shifts),
:math:`P` is the periodic coordinate mapping and :math:`c_{11}`, :math:`c_{12}`
and :math:`c_2` are unknown constant values - the unknown DOFs in
:math:`\Omega_{m11}`, :math:`\Omega_{m12}` and :math:`\Omega_{m2}` are replaced
by the integral mean values.

View the results using::

  sfepy-view square_quad.vtk -f u1:wu1:p0 1:vw:p0 u2:wu2:p1 1:vw:p1
"""
import numpy as nm

import sfepy.discrete.fem.periodic as per
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/square_quad.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

def get_shift1(ts, coors, region):
    val = 0.1 * coors[:, 0]

    return val

def get_shift2(ts, coors, region):
    val = nm.empty_like(coors[:, 1])
    val.fill(0.3)

    return val

functions = {
    'get_shift1' : (get_shift1,),
    'get_shift2' : (get_shift2,),
    'match_y_line' : (per.match_y_line,),
    'match_x_line' : (per.match_x_line,),
}

fields = {
    'scalar1': ('real', 1, 'Omega1', 1),
    'scalar2': ('real', 1, 'Omega2', 1),
}

materials = {
}

variables = {
    'u1' : ('unknown field', 'scalar1', 0),
    'v1' : ('test field', 'scalar1', 'u1'),
    'u2' : ('unknown field', 'scalar2', 1),
    'v2' : ('test field', 'scalar2', 'u2'),
}

regions = {
    'Omega1' : 'cells of group 1',
    'Omega2' : 'cells of group 2',
    'Omega_m1' : 'r.Omega1 -v (r.Gamma +s vertices of surface)',
    'Omega_m11' : 'r.Omega_m1 *v vertices in (x < 0)',
    'Omega_m12' : 'r.Omega_m1 *v vertices in (x > 0)',
    'Omega_m2' : 'r.Omega2 -v (r.Gamma +s vertices of surface)',
    'Left' : ('vertices in (x < -0.499)', 'facet'),
    'Right' : ('vertices in (x > 0.499)', 'facet'),
    'Bottom' : ('vertices in ((y < -0.499))', 'facet'),
    'Top' : ('vertices in ((y > 0.499))', 'facet'),
    'Gamma' : ('r.Omega1 *v r.Omega2 -v vertices of surface', 'facet'),
    'Gamma1' : ('copy r.Gamma', 'facet', 'Omega1'),
    'Gamma2' : ('copy r.Gamma', 'facet', 'Omega2'),
}

ebcs = {
    'fix1' : ('Top', {'u2.all' : 1.0}),
    'fix2' : ('Bottom', {'u1.all' : 0.0}),
}

lcbcs = {
    'shifted1' : (('Gamma1', 'Gamma2'),
                  {'u1.all' : 'u2.all'},
                  'match_x_line', 'shifted_periodic',
                  'get_shift1'),
    'shifted2' : (('Left', 'Right'),
                  {'u1.all' : 'u1.all'},
                  'match_y_line', 'shifted_periodic',
                  'get_shift2'),
    'mean11' : ('Omega_m11', {'u1.all' : None}, None, 'integral_mean_value'),
    'mean12' : ('Omega_m12', {'u1.all' : None}, None, 'integral_mean_value'),
    'mean2' : ('Omega_m2', {'u2.all' : None}, None, 'integral_mean_value'),
}

integrals = {
    'i1' : 2,
}

equations = {
    'eq1' : """
        dw_laplace.i1.Omega1(v1, u1) = 0
    """,
    'eq2' : """
        dw_laplace.i1.Omega2(v2, u2) = 0
    """,
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}
