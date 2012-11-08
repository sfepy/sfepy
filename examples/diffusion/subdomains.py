r"""
Two Laplace equations, one in :math:`\Omega_1`, other in
:math:`\Omega_2`, connected on the interface region :math:`\Gamma_{12}`
using traces of variables.

Find :math:`p_1`, :math:`p_2` such that:

.. math::
    \int_{\Omega_1} c_{11} \nabla q_1 \cdot \nabla p_1
    = \int_{\Gamma_{12}} q_1 (p_1 - \tr(p_2) - c_{12})
    \;, \quad \forall q_1 \;,

    \int_{\Omega_2} c_{21} \nabla q_2 \cdot \nabla p_2
    = \int_{\Gamma_{12}} q_2 (p_2 - \tr(p_1) - c_{22})
    \;, \quad \forall q_2 \;.
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/special/two_rectangles.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

regions = {
    'Omega1' : ('elements of group 1', {}),
    'Omega2' : ('elements of group 2', {}),
    'Gamma1' : ('nodes in (y > 0.29)', {}),
    'Gamma2' : ('nodes in (y < 0.01)', {}),
    'Gamma12' : ('r.Omega1 *n r.Omega2', {}),
    'Gamma12_1' : ('copy r.Gamma12', {'forbid' : 'group 2'}),
    'Gamma12_2' : ('copy r.Gamma12', {'forbid' : 'group 1'}),
}

materials = {
    'one1' : ({'one' : 1.0},),
    'one2' : ({'one' : 1.0},),
    'jump1'  : ({'val' : 1.0},),
    'jump2'  : ({'val' : 1.0},),
}

fields = {
    'pressure1' : ('real', 1, 'Omega1', 1),
    'pressure2' : ('real', 1, 'Omega2', 1),
}

variables = {
    'p1'   : ('unknown field', 'pressure1', 0),
    'q1'   : ('test field',    'pressure1', 'p1'),
    'p2'   : ('unknown field', 'pressure2', 1),
    'q2'   : ('test field',    'pressure2', 'p2'),
}

ebcs = {
    'p1' : ('Gamma1', {'p1.all' : 0.0}),
    'p2' : ('Gamma2', {'p2.0' : 1.0}),
}

integrals = {
    'ivol' : ('v', 2),
    'isurf' : ('s', 2),
}

equations = {
    'eq_1' :
    """dw_laplace.ivol.Omega1( one1.one, q1, p1 )
     + dw_jump.isurf.Gamma12_1( jump1.val, q1, p1, tr(p2) )
     = 0""",
    'eq_2' :
    """dw_laplace.ivol.Omega2( one2.one, q2, p2 )
     + dw_jump.isurf.Gamma12_2( jump2.val, q2, p2, tr(p1) )
     = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {'i_max'      : 1,
                               'eps_a'      : 1e-1,
                               'eps_r'      : 1.0,
                               'macheps'   : 1e-16,
                               # Linear system error < (eps_a * lin_red).
                               'lin_red'    : 1e-1,
                               'ls_red'     : 0.1,
                               'ls_red_warp' : 0.001,
                               'ls_on'      : 1.1,
                               'ls_min'     : 1e-5,
                               'check'     : 0,
                               'delta'     : 1e-6,
                               'is_plot'    : False,
                               # 'nonlinear' or 'linear' (ignore i_max)
                               'problem'   : 'nonlinear',
                               } )
}
