"""
A script demonstrating the solution of the scalar Helmholtz equation for a situation inspired by the physical
problem of WiFi propagation in an apartment.

The boundary conditions are defined as perfect radiation conditions.

The PDE for this physical process implies to find :math:`E(x, t)` for :math:`x
\in \Omega` such that:

.. math::
    \left\lbrace
    \begin{aligned}
      \Delta E+ k²n(x)²E = -f_s && \forall x \in \Omega \\
      \partial_n E(x)-ikn(x)E(x)=0 && \forall x \text{ on } \partial \Omega
    \end{aligned}
    \right.
"""
import numpy as nm
from sfepy import data_dir
from sfepy.discrete.fem.utils import refine_mesh

f = 2.4 * 1e9  # change this to 2.4 or 5 if you want to simulate frequencies typically used by your Wifi
c0 = 3e8  # m/s
k = 2 * nm.pi * f / c0

refinement_level = 3
filename_mesh = data_dir + "/meshes/2d/helmholtz_apartment.vtk"
filename_mesh = refine_mesh(filename_mesh, refinement_level)

regions = {
    'Omega': 'all',
    'Walls': ('cells of group 1'),
    'Air': ('cells of group 2'),
    'Source': ('cells of group 3'),
    'Gamma': ('vertices of surface', 'facet'),
}

# air and walls have different material parameters, hence the 1. and 2.4 factors
materials = {
    'air': ({'kn_square': (k * 1.) ** 2},),
    'walls': ({'kn_square': (k * 2.4) ** 2,
               'kn': 1j * k * 2.4},),
    'source_func': 'source_func',
}


def source_func(ts, coors, mode=None, **kwargs):
    """The source function for the antenna."""
    epsilon = 0.1
    x_center = nm.array([-3.23, -1.58])
    if mode == 'qp':
        dists = nm.abs(nm.sum(nm.square(coors - x_center), axis=1))
        val = 3 / nm.pi / epsilon ** 2 * (1 - dists / epsilon)
        return {'val': val[:, nm.newaxis, nm.newaxis]}


functions = {
    'source_func': (source_func,),
}

fields = {
    'electric_field': ('complex', 1, 'Omega', 1),
}

variables = {
    'E': ('unknown field', 'electric_field', 1),
    'q': ('test field', 'electric_field', 'E'),
}

ebcs = {
}

integrals = {
    'i': 2,
}

equations = {
    'Helmholtz equation':
        """- dw_laplace.i.Walls(q, E)
        - dw_laplace.i.Air(q, E)
        - dw_laplace.i.Source(q, E)
        + dw_dot.i.Walls(walls.kn_square, q, E)
        + dw_dot.i.Air(air.kn_square, q, E)
        + dw_dot.i.Source(air.kn_square, q, E)
        + dw_dot.i.Gamma(walls.kn, q, E)
        = - dw_integrate.i.Source(source_func.val, q)
        """
}

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
    }),
}
