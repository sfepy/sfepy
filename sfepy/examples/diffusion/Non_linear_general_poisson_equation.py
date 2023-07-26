r"""
Example usage of the non linear diffusion and non linear volume terms.

The example is an adaptation of:
https://sfepy.org/doc-devel/examples/diffusion-poisson_field_dependent_material.html

Find :math:`T(t)` for :math:`t \in [0, t_{\rm final}]` such that:

.. math::
   \int_{\Omega} c(T) \nabla s \cdot \nabla T + \int_{\Omega} g(T) \cdot s
    = 0
    \;, \quad \forall s \;.

where :math:`c(T) and g(T)` are the :math:`T` dependent coefficients.
"""
from __future__ import absolute_import
from sfepy import data_dir
from sfepy.base.base import output

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

def get_conductivity(T):
    """
    Calculates the conductivity as  2+10*(T+2) and returns it.
    """
    val = 2 + 10 * (T + 2)
    return val
def d_get_conductivity(T):
    """
    Calculates the derivative of the conductivity and returns it.
    """
    val = 10
    return val

def nl_src(T):
    """
    Calculates a non linear source term as 10*(T**2 +2) and returns it.
    """
    val = 10 * (T**2 + 2)
    return val

def d_nl_src(T):
    """
    Calculates the derivative of the nl_src term and returns it.
    """
    val = 10*(2*T)
    return val

materials = {
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    'T' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 'T'),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.099999)', 'facet'),
}

ebcs = {
    'T1' : ('Gamma_Left', {'T.0' : 2.0}),
    'T2' : ('Gamma_Right', {'T.0' : -2.0}),
}

functions = {
    'get_conductivity' : (get_conductivity,),
    'd_get_conductivity' : (d_get_conductivity,),
    'nl_src' : (nl_src,),
    'd_nl_src' : (d_nl_src,),
}

ics = {
    'ic' : ('Omega', {'T.0' : 0.0}),
}

integrals = {
    'i' : 1,
}

equations = {
    'Temperature' : """dw_nl_diffusion.i.Omega(get_conductivity, d_get_conductivity, s, T ) + dw_volume_nvf.i.Omega(nl_src, d_nl_src,s,T) = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 10,
        'eps_a' : 1e-10,
        'eps_r' : 1.0,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'save_times' : 'all',
}
