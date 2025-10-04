r"""
Poisson equation solved in a single patch NURBS domain using the isogeometric
analysis (IGA) approach.

Find :math:`t` such that:

.. math::
    \int_{\Omega} c \nabla s \cdot \nabla t
    =  \int_{\Omega_0} f s
    \;, \quad \forall s \;.

Try setting the Dirichlet boundary condition (ebcs) on various sides of the
domain (``'Gamma1'``, ..., ``'Gamma4'``).

View the results using::

  sfepy-view patch2d.vtk -f t:wt:f0.4 1:vw
"""
from sfepy import data_dir

filename_domain = data_dir + '/meshes/iga/patch2d.iga'

materials = {
    'm' : ({'c' : 1.0, 'f' : -10.0},),
}

regions = {
    'Omega' : 'all',
    'Omega_0' : 'vertices in (x > 1.5)',
    'Gamma1' : ('vertices of set xi00', 'facet'),
    'Gamma2' : ('vertices of set xi01', 'facet'),
    'Gamma3' : ('vertices of set xi10', 'facet'),
    'Gamma4' : ('vertices of set xi11', 'facet'),
}

fields = {
    'temperature' : ('real', 1, 'Omega', None, 'H1', 'iga'),
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

ebcs = {
    't1' : ('Gamma3', {'t.0' : 2.0}),
    't2' : ('Gamma4', {'t.0' : -2.0}),
}

integrals = {
    'i' : 3,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega(m.c, s, t)
                       = dw_volume_lvf.i.Omega_0(m.f, s)"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}
