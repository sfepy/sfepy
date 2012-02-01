r"""
Acoustic pressure distribution.

This example shows how to solve a problem in complex numbers, note the
'accoustic_pressure' field definition.

Find :math:`p` such that:

.. math::
    c^2 \int_{\Omega} \nabla q \cdot \nabla p
    - w^2 \int_{\Omega} q p
    - i w c \int_{\Gamma_{out}} q p
    = i w c^2 \rho v_n \int_{\Gamma_{in}} q
    \;, \quad \forall q \;.
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/special/two_rectangles.mesh'

v_n = 1.0 # m/s
w = 1000.0
c = 343.0 # m/s
rho = 1.55 # kg/m^3

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

materials = {
    'one' : ({'one' : 1.0},),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_in' : ('nodes in (x < 0.01)', {}),
    'Gamma_out' : ('nodes in (x > 0.99)', {}),
}

fields = {
    'accoustic_pressure' : ('complex', 1, 'Omega', 1),
}

variables = {
    'p'   : ('unknown field', 'accoustic_pressure', 0),
    'q'   : ('test field',    'accoustic_pressure', 'p'),
}

ebcs = {
}

integrals = {
    'ivol' : ('v', 'gauss_o2_d2'),
    'isurf' : ('s', 'gauss_o2_d1'),
}

equations = {
    'Acoustic pressure' :
    """%s * dw_laplace.ivol.Omega( one.one, q, p )
    - %s * dw_mass_scalar.ivol.Omega( q, p )
    - %s * dw_surface_mass_scalar.isurf.Gamma_out( q, p )
    = %s * dw_surface_integrate.isurf.Gamma_in( q )"""
    % (c*c, w*w, 1j*w*c, 1j*w*c*c*rho*v_n)
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
