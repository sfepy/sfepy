# mixed formulation

# 07.08.2009
#!
#! Homogenization: Linear Elasticity
#! =================================
#$ \centerline{Example input file, \today}

#! Homogenization of heterogeneous linear elastic material - mixed formulation

from sfepy.fem.periodic import *
from sfepy.mechanics.matcoefs import stiffness_tensor_youngpoisson, stiffness_tensor_youngpoisson_mixed, bulk_modulus_youngpoisson
from sfepy.homogenization.utils import define_box_regions, get_box_volume
import sfepy.homogenization.coefs_elastic as ce
from sfepy import data_dir
from sfepy.base.base import Struct, debug
from sfepy.homogenization.recovery import compute_micro_u, compute_stress_strain_u, compute_mac_stress_part, add_stress_p

def recovery_le( pb, corrs, macro ):
    
    out = {}
    dim = corrs['corrs_le']['u_00'].shape[1]
    mic_u = - compute_micro_u( corrs['corrs_le'], macro['strain'], 'u', dim )
    mic_p = - compute_micro_u( corrs['corrs_le'], macro['strain'], 'p', dim )
    
    out['u_mic'] = Struct( name = 'output_data',
                           mode = 'vertex', data = mic_u,
                           var_name = 'u', dofs = None )
    out['p_mic'] = Struct( name = 'output_data',
                           mode = 'cell', data = mic_p[:,nm.newaxis,
                                                       :,nm.newaxis],
                           var_name = 'p', dofs = None )

    stress_Ym, strain_Ym = compute_stress_strain_u( pb, 'i1', 'Ym', 'matrix.D', 'u', mic_u )
    stress_Ym += compute_mac_stress_part( pb, 'i1', 'Ym', 'matrix.D', 'u', macro['strain'] )
    add_stress_p( stress_Ym, pb, 'i1', 'Ym', 'p', mic_p )    
    stress_Yc, strain_Yc = compute_stress_strain_u( pb, 'i1', 'Yc', 'reinf.D', 'u', mic_u )
    stress_Yc += compute_mac_stress_part( pb, 'i1', 'Yc', 'reinf.D', 'u', macro['strain'] )
    add_stress_p( stress_Yc, pb, 'i1', 'Yc', 'p', mic_p )

    strain = macro['strain'] + strain_Ym + strain_Yc
    stress = stress_Ym + stress_Yc

    out['cauchy_strain'] = Struct( name = 'output_data',
                                   mode = 'cell', data = strain,
                                   dofs = None )
    out['cauchy_stress'] = Struct( name = 'output_data',
                                   mode = 'cell', data = stress,
                                   dofs = None )
    return out

#! Mesh
#! ----
dim = 3
filename_mesh = data_dir + '/meshes/3d/matrix_fiber.mesh'
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
    'matrix' : ('Ym', {'D' : stiffness_tensor_youngpoisson_mixed( dim, 0.7e9, 0.4 ),
                       'gamma' : bulk_modulus_youngpoisson( 0.7e9, 0.4 )}),
    'reinf' : ('Yc', {'D' : stiffness_tensor_youngpoisson_mixed( dim, 70.0e9, 0.2 ),
                      'gamma' : bulk_modulus_youngpoisson( 70.0e9, 0.2 )}),
}
gamma_m = bulk_modulus_youngpoisson( 0.7e9, 0.4 )
gamma_c = bulk_modulus_youngpoisson( 70.0e9, 0.2 )
#! Fields
#! ------
#! Scalar field for corrector basis functions.
geom = {2 : '2_4', 3 : '3_8'}[dim]
fields = {
    'corrector_u' : ((dim,1), 'real', 'Y', {'Y' : '%s_Q1' % geom}),
    'corrector_p' : ((1,1), 'real', 'Y', {'Y' : '%s_Q0' % geom}),
}
#! Variables
#! ---------
#! Unknown and corresponding test variables. Parameter fields
#! used for evaluation of homogenized coefficients.
variables = {
    'u'     : ('unknown field', 'corrector_u', 0),
    'v'     : ('test field', 'corrector_u', 'u'),
    'p'     : ('unknown field', 'corrector_p', 1),
    'q'     : ('test field', 'corrector_p', 'p'),
    'Pi'    : ('parameter field', 'corrector_u', 'u'),
    'Pi1u' : ('parameter field', 'corrector_u', 'u'),
    'Pi2u' : ('parameter field', 'corrector_u', 'u'),
    'Pi1p' : ('parameter field', 'corrector_p', 'p'),
    'Pi2p' : ('parameter field', 'corrector_p', 'p'),
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
    'coefs' : 'coefs',
    'requirements' : 'requirements',
    'ls' : 'ls', # linear solver to use
    'volume' : { #'variables' : ['u'],
                 #'expression' : 'd_volume.i1.Y( u )',
                 'value' : get_box_volume( dim, region_lbn, region_rtf ),
                 },    
    'output_dir' : 'output',
    'coefs_filename' : 'output/coefs_le_up.h5',
    'recovery_hook' : 'recovery_le',
}
#! Equations
#! ---------
#! Equations for corrector functions.
equation_corrs = {
    'balance_of_forces' :
    """  dw_lin_elastic.i1.Ym( matrix.D, v, u )
       + dw_lin_elastic.i1.Yc( reinf.D, v, u )
       - dw_stokes.i1.Ym( v, p )
       - dw_stokes.i1.Yc( v, p ) =
       - dw_lin_elastic.i1.Ym( matrix.D, v, Pi )
       - dw_lin_elastic.i1.Yc( reinf.D, v, Pi )""",
    'pressure constraint' :
    """- dw_stokes.i1.Ym( u, q )
       - dw_stokes.i1.Yc( u, q )
       - %e * dw_mass_scalar.i1.Ym( q, p )
       - %e * dw_mass_scalar.i1.Yc( q, p ) = 
         dw_stokes.i1.Ym( Pi, q )
       + dw_stokes.i1.Yc( Pi, q )""" % (1.0/gamma_m, 1.0/gamma_c),
}
#! Expressions for homogenized linear elastic coefficients.
expr_coefs = {
    'Q1' : 'dw_lin_elastic.i1.Ym( matrix.D, Pi1u, Pi2u ) +\
      dw_lin_elastic.i1.Yc( reinf.D, Pi1u, Pi2u )',
    'Q2' : '%e * dw_mass_scalar.i1.Ym( Pi1p, Pi2p ) +\
      %e * dw_mass_scalar.i1.Yc( Pi1p, Pi2p )' % (1.0/gamma_m, 1.0/gamma_c),
}
#! Coefficients
#! ------------
#! Definition of homogenized acoustic coefficients.
coefs = {
    'elastic_u' : {'requires' : ['pis', 'corrs_rs'],
                   'variables' : ['Pi1u', 'Pi2u'],
                   'expression' : expr_coefs['Q1'],
                   'class' : ce.ElasticCoef,
                   },
    'elastic_p' : {'requires' : ['corrs_rs'],
                   'variables' : ['Pi1p', 'Pi2p'],
                   'expression' : expr_coefs['Q2'],
                   'class' : ce.ElasticPCoef,
                   },
    'D' : {'requires' : ['elastic_u', 'elastic_p'],
           'class' : 'sum',
           },
    'filenames' : {},
}
# requirements for elastic homog. coefficients
requirements = {
    'pis' : { 'variables' : ['u'],
              'class' : ce.ShapeDimDim,
              },
    'corrs_rs' : { 'requires' : ['pis'],
                   'variables' : ['u', 'v', 'p', 'q', 'Pi'],
                   'ebcs' : ['fixed_u'],
                   'epbcs' : all_periodic,
                   'equations' : equation_corrs,
                   'class' : ce.CorrectorsElasticRS,
                   'save_name' : 'corrs_le',
                   'dump_variables' : ['u', 'p'],
                   },
}
#! Solvers
#! -------
#! Define linear and nonlinear solver.
solvers = {
    'ls' : ('ls.umfpack', {}),
    'newton' : ('nls.newton', {'i_max' : 1,
                               'eps_a' : 1e-4,
                               'problem' : 'nonlinear', })
}
#! FE assembling parameters
#! ------------------------
#! 'chunk_size' determines maximum number of elements to assemble in one C
#! function call. Higher values mean faster assembling, but also more memory
#! usage.
fe = {
    'chunk_size' : 1000
}
