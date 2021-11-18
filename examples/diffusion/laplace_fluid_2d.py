r"""
A Laplace equation that mimics the flow of "dry water" around an obstacle shaped like a Citroen CX.


Description
-----------

As discussed e.g. in the Feynman lectures in Section 12-5 of Volume 2 (https://www.feynmanlectures.caltech.edu/II_12.html#Ch12-S5),
the flow of an irrotational and incompressible fluid can be modelled with a
potential :math:`\ul{v} = \ul{grad}(\phi)` that obeys

.. math::
    div(\ul{grad}\,\phi) = \Delta\,\phi = 0

The weak formulation for this problem is to find :math:`\phi` such that:

.. math::
    \int_{\Omega} \nabla \psi \cdot \nabla \phi
    = \int{\Gamma_{left} \ul{v}_0 \cdot n \psi + \int{\Gamma_{right} \ul{v}_0 \cdot n \psi
    \;, \quad \forall \psi \;,

Usage examples
--------------
Run with simple.py script::

    python3 simple.py examples/diffusion/laplace_fluid_2d.py
    python3 resview.py citroen.vtk

"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/citroen.h5'

materials = {
    'm' : ({'sigma' : 1.0},),
    'f': ({'left': -1.0,
           'right': 1.0},),
}


regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.1)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 1919.9)', 'facet'),
}

fields = {
    'u' : ('real', 1, 'Omega', 1),
}

variables = {
    'phi' : ('unknown field', 'u', 0),
    'psi' : ('test field',    'u', 'phi'),
}

ebcs = {
}

integrals = {
    'i' : 2,
}

equations = {
    'Laplace equation' :
    """dw_laplace.i.Omega( m.sigma, psi, phi )
     = dw_integrate.i.Gamma_Left(f.left, psi)
     + dw_integrate.i.Gamma_Right(f.right, psi)"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}
