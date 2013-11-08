r"""
Transient Laplace equation.

The same example as time_poisson.py, but using the short syntax of keywords,
and explicit time-stepping.

Find :math:`T(t)` for :math:`t \in [0, t_{\rm final}]` such that:

.. math::
    \int_{\Omega} s \pdiff{T}{t}
    + \int_{\Omega} c \nabla s \cdot \nabla T
    = 0
    \;, \quad \forall s \;.
"""
from sfepy import data_dir

from examples.diffusion.time_poisson import get_ic

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

materials = {
    'coef' : ({'val' : 0.01},),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.099999)', 'facet'),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    'T' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 'T'),
}

ebcs = {
    't1' : ('Gamma_Left', {'T.0' : 2.0}),
    't2' : ('Gamma_Right', {'T.0' : -2.0}),
}

ics = {
    'ic' : ('Omega', {'T.0' : 'get_ic'}),
}

functions = {
    'get_ic' : (get_ic,),
}

integrals = {
    'i' : 1,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega( coef.val, s, T ) = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'ts' : ('ts.explicit', {
        't0' : 0.0,
        't1' : 1.0,
        'dt' : 0.00002,
        'n_step' : None,
        'mass' : 'dw_volume_dot.i.Omega( s, T )',
        'lumped' : False, # If True, lump mass matrix so that it is diagonal.
    }),
}

options = {
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : 100,
    'output_format' : 'h5',
}
