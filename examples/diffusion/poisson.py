r"""
Laplace equation using the long syntax of keywords.

See the tutorial section :ref:`poisson-example-tutorial` for a detailed
explanation. See :ref:`diffusion-poisson_short_syntax` for the short syntax
version.

Find :math:`t` such that:

.. math::
    \int_{\Omega} c \nabla s \cdot \nabla t
    = 0
    \;, \quad \forall s \;.
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

material_2 = {
    'name' : 'coef',
    'values' : {'val' : 1.0},
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'cells of group 6',
}

region_03 = {
    'name' : 'Gamma_Left',
    'select' : 'vertices in (x < 0.00001)',
    'kind' : 'facet',
}

region_4 = {
    'name' : 'Gamma_Right',
    'select' : 'vertices in (x > 0.099999)',
    'kind' : 'facet',
}

field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0, # order in the global vector of unknowns
}

variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 't',
}

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Left',
    'dofs' : {'t.0' : 2.0},
}

ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Right',
    'dofs' : {'t.0' : -2.0},
}

integral_1 = {
    'name' : 'i',
    'order' : 2,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega( coef.val, s, t ) = 0"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
    'method' : 'auto',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}
