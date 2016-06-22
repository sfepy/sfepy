r"""
Poisson equation with Neumann boundary conditions on a part of the boundary.

Find :math:`t` such that:

.. math::
    \int_{\Omega} c \nabla s \cdot \nabla t
    = \int_{\Gamma_N} s g
    \;, \quad \forall s \;,

where :math:`g` is the given flux, :math:`g = \nabla T \cdot \ul{n}`. See the
tutorial section :ref:`poisson-weak-form-tutorial` for a detailed explanation.
"""
from __future__ import absolute_import
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

materials = {
    'flux' : ({'val' : -50.0},),
    'coef' : ({'val' : 2.0},),
}

regions = {
    'Omega' : 'all', # or 'cells of group 6'
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.099999)', 'facet'),
    'Gamma_N' : ('vertices of surface -f (r.Gamma_Left +v r.Gamma_Right)',
                 'facet'),
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

equations = {
    'Temperature' : """
           dw_laplace.2.Omega(coef.val, s, t)
         = dw_surface_integrate.2.Gamma_N(flux.val, s)
    """
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-10,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}
