r"""
Linear elasticity with comments.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    = 0
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.
"""
#!
#! Linear Elasticity
#! =================
#$ \centerline{Example input file, \today}

#! This file models a cylinder that is fixed at one end while the
#! second end has a specified displacement of 0.01 in the x direction
#! (this boundary condition is named Displaced). There is also a specified
#! displacement of 0.005 in the z direction for points in
#! the region labeled SomewhereTop. This boundary condition is named
#! PerturbedSurface.  The region SomewhereTop is specified as those nodes for
#! which
#!              (z > 0.017) & (x > 0.03) & (x <  0.07).
#! The output is the displacement for each node, saved by default to
#! simple_out.vtk. The material is linear elastic and its properties are
#! specified as Lame parameters (see
#! http://en.wikipedia.org/wiki/Lam%C3%A9_parameters)
#!

#! Mesh
#! ----
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'
#! Regions
#! -------
#! Whole domain 'Omega', left and right ends.
regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < 0.001)', 'facet'),
    'Right' : ('vertices in (x > 0.099)', 'facet'),
    'SomewhereTop' : ('vertices in (z > 0.017) & (x > 0.03) & (x < 0.07)',
                      'vertex'),
}
#! Materials
#! ---------
#! The linear elastic material model is used. Properties are
#! specified as Lame parameters.
materials = {
    'solid' : ({'lam' : 1e1, 'mu' : 1e0},),
}
#! Fields
#! ------
#! A field is used to define the approximation on a (sub)domain
#! A displacement field (three DOFs/node) will be computed on a region
#! called 'Omega' using P1 (four-node tetrahedral) finite elements.
fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}
#! Integrals
#! ---------
#! Define the integral type Volume/Surface and quadrature rule
#! (here: dim=3, order=1).
integrals = {
    'i1' : ('v', 1),
}
#! Variables
#! ---------
#! One field is used for unknown variable (generate discrete degrees
#! of freedom) and the seccond field for the corresponding test variable of
#! the weak formulation.
variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}
#! Boundary Conditions
#! -------------------
#! The left end of the cilinder is fixed (all DOFs are zero) and
#! the 'right' end has non-zero displacements only in the x direction.
ebcs = {
    'Fixed' : ('Left', {'u.all' : 0.0}),
    'Displaced' : ('Right', {'u.0' : 0.01, 'u.[1,2]' : 0.0}),
    'PerturbedSurface' : ('SomewhereTop', {'u.2' : 0.005}),
}
#! Equations
#! ---------
#! The weak formulation of the linear elastic problem.
equations = {
    'balance_of_forces' :
    """dw_lin_elastic_iso.i1.Omega( solid.lam, solid.mu, v, u ) = 0""",
}
#! Solvers
#! -------
#! Define linear and nonlinear solver.
#! Even linear problems are solved by a nonlinear solver (KISS rule) - only one
#! iteration is needed and the final rezidual is obtained for free.
solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                { 'i_max'      : 1,
                  'eps_a'      : 1e-10,
                  'eps_r'      : 1.0,
                  'macheps'   : 1e-16,
                  # Linear system error < (eps_a * lin_red).
                  'lin_red'    : 1e-2,                
                  'ls_red'     : 0.1,
                  'ls_red_warp' : 0.001,
                  'ls_on'      : 1.1,
                  'ls_min'     : 1e-5,
                  'check'     : 0,
                  'delta'     : 1e-6,
                  'is_plot'    : False,
                  # 'nonlinear' or 'linear' (ignore i_max)
                  'problem'   : 'nonlinear'}),
}
