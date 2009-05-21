# 30.04.2009
#!
#! Linear Elasticity
#! =================
#$ \centerline{Example input file, \today}

#! This file models a cylinder that is fixed at one end while the
#! second end has a specified displacement of 0.02 in the x direction
#! (this boundary condition is named PerturbedSurface).
#! The output is the displacement for each node, saved by default to
#! simple_out.vtk. The material is linear elastic and its properties are
#! specified as Lame parameters (see
#! http://en.wikipedia.org/wiki/Lam%C3%A9_parameters)
#!

#! Mesh
#! ----
filename_mesh = '../database/simple.vtk'
#! Regions
#! -------
#! Whole domain 'Omega', left and right ends.
regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001)', {}),
    'Right' : ('nodes in (x > 0.099)', {}),
}
#! Materials
#! ---------
#! The linear elastic material model is used. Properties are
#! specified as Lame parameters.
materials = {
    'solid' : ('here', 'Omega', {'lame' : {'lambda' : 1e1, 'mu' : 1e0}}),
}
#! Fields
#! ------
#! A field is used to define the approximation on a (sub)domain
#! A displacement field (three DOFs/node) will be computed on a region
#! called 'Omega' using P1 (four-node tetrahedral) finite elements.
fields = {
    '3_displacement': ((3,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
}
#! Integrals
#! ---------
#! Define the integral type Volume/Surface and quadrature rule
#! (here: dim=3, order=1).
integrals = {
    'i1' : ('v', 'gauss_o1_d3'),
}
#! Variables
#! ---------
#! One field is used for unknown variable (generate discrete degrees
#! of freedom) and the seccond field for the corresponding test variable of
#! the weak formulation.
variables = {
    'u' : ('unknown field', '3_displacement', 0),
    'v' : ('test field', '3_displacement', 'u'),
}
#! Boundary Conditions
#! -------------------
#! The left end of the cilinder is fixed (all DOFs are zero) and
#! the 'right' end has non-zero displacements only in the x direction.
ebcs = {
    'Fixed' : ('Left', {'u.all' : 0.0}),
    'PerturbedSurface' : ('Right', {'u.0' : 0.02, 'u.1' : 0.0, 'u.2' : 0.0}),
}
#! Equations
#! ---------
#! The weak formulation of the linear elastic problem.
equations = {
    'balance_of_forces' : """dw_lin_elastic_iso.i1.Omega( solid.lame, v, u ) = 0""",
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
#! FE assembling parameters
#! ------------------------
#! 'chunk_size' determines maximum number of elements to assemble in one C
#! function call. Higher values mean faster assembling, but also more memory
#! usage.
fe = {
    'chunk_size' : 1000
}
