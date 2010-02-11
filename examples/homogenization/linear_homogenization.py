# 04.08.2009
#!
#! Homogenization: Linear Elasticity
#! =================================
#$ \centerline{Example input file, \today}

#! Homogenization of heterogeneous linear elastic material

from sfepy.fem.periodic import *
from sfepy.mechanics.matcoefs import stiffness_tensor_youngpoisson
from sfepy.homogenization.utils import define_box_regions
import sfepy.homogenization.coefs_elastic as ce
from sfepy import top_dir
#! Mesh
#! ----
filename_mesh = top_dir + '/meshes/3d/matrix_fiber.mesh'
dim = 3
region_lbn = (0, 0, 0)
region_rtf = (1, 1, 1)
#! Regions
#! -------
#! Regions, edges, ... 
regions = {
    'Y' : ('all', {}),
    'Ym' : ('elements of group 1', {}),
    'Yc' : ('elements of group 2', {}),
}
regions.update( define_box_regions( dim, region_lbn, region_rtf ) )
#! Materials
#! ---------
materials = {
    'matrix' : ('Ym', {'D' : stiffness_tensor_youngpoisson( dim, 0.7e9, 0.4 ) }),
    'reinf' : ('Yc', {'D' : stiffness_tensor_youngpoisson( dim, 70.0e9, 0.2 ) }),
}
#! Fields
#! ------
#! Scalar field for corrector basis functions.
geom = {2 : '2_4', 3 : '3_8'}[dim]
fields = {
    'corrector' : ((dim,1), 'real', 'Y', {'Y' : '%s_Q1' % geom}),
}
#! Variables
#! ---------
#! Unknown and corresponding test variables. Parameter fields
#! used for evaluation of homogenized coefficients.
variables = {
    'u'     : ('unknown field', 'corrector', 0),
    'v'     : ('test field', 'corrector', 'u'),
    'Pi'    : ('parameter field', 'corrector', 'u'),
    'Pi1' : ('parameter field', 'corrector', 'u'),
    'Pi2' : ('parameter field', 'corrector', 'u'),
}
#! Functions
functions = {
    'match_x_plane' : (match_x_plane,),
    'match_y_plane' : (match_y_plane,),
    'match_z_plane' : (match_z_plane,),
}
#! Boundary Conditions
#! -------------------
#! Fixed nodes.
ebcs = {
    'fixed_u' : ('Corners', {'u.all' : 0.0}),
}
if dim == 3:
    epbcs = {
        'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'}, 'match_x_plane'),
        'periodic_y' : (['Near', 'Far'], {'u.all' : 'u.all'}, 'match_y_plane'),
        'periodic_z' : (['Top', 'Bottom'], {'u.all' : 'u.all'}, 'match_z_plane'),
}
else:
    epbcs = {
        'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'}, 'match_x_plane'),
        'periodic_y' : (['Bottom', 'Top'], {'u.all' : 'u.all'}, 'match_y_plane'),
    }
all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:dim] ]
#! Integrals
#! ---------
#! Define the integral type Volume/Surface and quadrature rule.
integrals = {
    'i1' : ('v', 'gauss_o2_d%d' % dim),
    'i2' : ('s', 'gauss_o2_d%d' % (dim-1)),
}
#! Options
#! -------
#! Various problem-specific options.
options = {
    'coef_info' : 'coefs', # homogenized coefficients to compute
    'coefs' : 'coefs',
    'requirements' : 'requirements',
    'ls' : 'ls', # linear solver to use
    'volume' : { 'variables' : ['u'],
                 'expression' : 'd_volume.i1.Y( u )' },
    'output_dir' : './output',
}
#! Equations
#! ---------
#! Equations for corrector functions.
equation_corrs = {
    'balance_of_forces' :
    """dw_lin_elastic.i1.Ym( matrix.D, v, u )
    + dw_lin_elastic.i1.Yc( reinf.D, v, u ) =
    - dw_lin_elastic.i1.Ym( matrix.D, v, Pi )
    - dw_lin_elastic.i1.Yc( reinf.D, v, Pi )"""
}
#! Expressions for homogenized linear elastic coefficients.
expr_coefs = """dw_lin_elastic.i1.Ym( matrix.D, Pi1, Pi2 )
+ dw_lin_elastic.i1.Yc( reinf.D, Pi1, Pi2 )"""
#! Coefficients
#! ------------
#! Definition of homogenized acoustic coefficients.
coefs = {
    'D' : {'requires' : ['pis', 'corrs_rs'],
           'variables' : ['Pi1', 'Pi2'],
           'expression' : expr_coefs,
           'class' : ce.ElasticCoef,}
}
# requirements for elastic homog. coefficients
requirements = {
    'pis' : { 'variables' : ['u'],
              'class' : ce.ShapeDimDim,
              },
    'corrs_rs' : { 'requires' : ['pis'],
                   'variables' : ['u', 'v', 'Pi'],
                   'ebcs' : ['fixed_u'],
                   'epbcs' : all_periodic,
                   'equations' : equation_corrs,
                   'class' : ce.CorrectorsElasticRS,
                   'save_name' : 'corrs_le',
                   'dump_variables' : ['u'],
                   },
}
#! Solvers
#! -------
#! Define linear and nonlinear solver.
solvers = {
    'ls' : ('ls.umfpack', {}),
    'newton' : ('nls.newton', {'i_max' : 1,
                               'eps_a' : 1e-4,
                               'problem' : 'nonlinear',
                               })
}
#! FE assembling parameters
#! ------------------------
#! 'chunk_size' determines maximum number of elements to assemble in one C
#! function call. Higher values mean faster assembling, but also more memory
#! usage.
fe = {
    'chunk_size' : 1000
}
