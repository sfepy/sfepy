r"""
The transient advection-diffusion equation with a given divergence-free
advection velocity.

Find :math:`u` such that:

.. math::
    \int_{\Omega} s \pdiff{u}{t}
    + \int_{\Omega} s \nabla \cdot \left(\ul{v} u \right)
    + \int_{\Omega} D \nabla s \cdot \nabla u
    = 0
    \;, \quad \forall s \;.

View the results using::

  sfepy-view square_tri2.*.vtk -f u:wu 1:vw
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/square_tri2.mesh'

regions = {
    'Omega' : 'all', # or 'cells of group 6'
    'Gamma_Left' : ('vertices in (x < -0.99999)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.99999)', 'facet'),
}

fields = {
    'concentration' : ('real', 1, 'Omega', 1),
}

variables = {
    'u' : ('unknown field', 'concentration', 0, 1),
    's' : ('test field',    'concentration', 'u'),
}

ebcs = {
    'u1' : ('Gamma_Left', {'u.0' : 2.0}),
    'u2' : ('Gamma_Right', {'u.0' : 0.0}),
}

# Units: D: 0.0001 m^2 / day, v: [0.1, 0] m / day -> time in days.
materials = {
    'm' : ({'D' : 0.0001, 'v' : [[0.1], [0.0]]},),
}

integrals = {
    'i' : 2,
}

equations = {
    'advection-diffusion' :
     """
       dw_dot.i.Omega(s, du/dt)
     + dw_advect_div_free.i.Omega(m.v, s, u)
     + dw_laplace.i.Omega(m.D, s, u)
     = 0
     """
}

solvers = {
    'ts' : ('ts.simple', {
        't0' : 0.0,
        't1' : 10.0,
        'dt' : None,
        'n_step' : 11, # Has precedence over dt.
        'verbose' : 1,
    }),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
    'ls' : ('ls.scipy_direct', {}),
}

options = {
    'ts' : 'ts',
    'nls' : 'newton',
    'ls' : 'ls',
    'save_times' : 'all',
}
