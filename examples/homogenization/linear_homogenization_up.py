# mixed formulation

# 07.08.2009
#!
#! Homogenization: Linear Elasticity
#! =================================
#$ \centerline{Example input file, \today}

#! Homogenization of heterogeneous linear elastic material - mixed formulation
import numpy as nm

import sfepy.fem.periodic as per
from sfepy.mechanics.matcoefs import stiffness_tensor_youngpoisson_mixed, bulk_modulus_youngpoisson
from sfepy.homogenization.utils import define_box_regions, get_box_volume
import sfepy.homogenization.coefs_base as cb

from sfepy import data_dir
from sfepy.base.base import Struct
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
    'matrix' : ({'D' : stiffness_tensor_youngpoisson_mixed( dim, 0.7e9, 0.4 ),
                 'gamma' : 1.0/bulk_modulus_youngpoisson( 0.7e9, 0.4 )},),
    'reinf' : ({'D' : stiffness_tensor_youngpoisson_mixed( dim, 70.0e9, 0.2 ),
                'gamma' : 1.0/bulk_modulus_youngpoisson( 70.0e9, 0.2 )},),
}
#! Fields
#! ------
#! Scalar field for corrector basis functions.
fields = {
    'corrector_u' : ('real', dim, 'Y', 1),
    'corrector_p' : ('real', 1, 'Y', 0),
}
#! Variables
#! ---------
#! Unknown and corresponding test variables. Parameter fields
#! used for evaluation of homogenized coefficients.
variables = {
    'u'     : ('unknown field', 'corrector_u'),
    'v'     : ('test field', 'corrector_u', 'u'),
    'p'     : ('unknown field', 'corrector_p'),
    'q'     : ('test field', 'corrector_p', 'p'),
    'Pi'    : ('parameter field', 'corrector_u', 'u'),
    'Pi1u' : ('parameter field', 'corrector_u', '(set-to-None)'),
    'Pi2u' : ('parameter field', 'corrector_u', '(set-to-None)'),
    'Pi1p' : ('parameter field', 'corrector_p', '(set-to-None)'),
    'Pi2p' : ('parameter field', 'corrector_p', '(set-to-None)'),
}
#! Functions
functions = {
    'match_x_plane' : (per.match_x_plane,),
    'match_y_plane' : (per.match_y_plane,),
    'match_z_plane' : (per.match_z_plane,),
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
    'coefs_filename' : 'coefs_le_up',
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
       - dw_mass_scalar_w.i1.Ym( matrix.gamma, q, p )
       - dw_mass_scalar_w.i1.Yc( reinf.gamma, q, p ) = 
         dw_stokes.i1.Ym( Pi, q )
       + dw_stokes.i1.Yc( Pi, q )""",
}
#! Expressions for homogenized linear elastic coefficients.
expr_coefs = {
    'Q1' : """  dw_lin_elastic.i1.Ym( matrix.D, Pi1u, Pi2u )
              + dw_lin_elastic.i1.Yc( reinf.D, Pi1u, Pi2u )""",
    'Q2' : """  dw_mass_scalar_w.i1.Ym( matrix.gamma, Pi1p, Pi2p )
              + dw_mass_scalar_w.i1.Yc( reinf.gamma, Pi1p, Pi2p )""",
}
#! Coefficients
#! ------------
#! Definition of homogenized acoustic coefficients.
def set_elastic_u(variables, ir, ic, mode, pis, corrs_rs):
    mode2var = {'row' : 'Pi1u', 'col' : 'Pi2u'}

    val = pis.states[ir, ic]['u'] + corrs_rs.states[ir, ic]['u']

    variables[mode2var[mode]].data_from_any(val)

coefs = {
    'elastic_u' : {
        'requires' : ['pis', 'corrs_rs'],
        'expression' : expr_coefs['Q1'],
        'set_variables' : set_elastic_u,
        'class' : cb.CoefSymSym,
    },
    'elastic_p' : {
        'requires' : ['corrs_rs'],
        'expression' : expr_coefs['Q2'],
        'set_variables' : [('Pi1p', 'corrs_rs', 'p'), ('Pi2p', 'corrs_rs', 'p')],
        'class' : cb.CoefSymSym,
    },
    'D' : {
        'requires' : ['c.elastic_u', 'c.elastic_p'],
        'class' : cb.CoefSum,
    },
    'filenames' : {},
}

requirements = {
    'pis' : {
        'variables' : ['u'],
        'class' : cb.ShapeDimDim,
    },
    'corrs_rs' : {
        'requires' : ['pis'],
        'ebcs' : ['fixed_u'],
        'epbcs' : all_periodic,
        'equations' : equation_corrs,
        'set_variables' : [('Pi', 'pis', 'u')],
        'class' : cb.CorrDimDim,
        'save_name' : 'corrs_le',
        'dump_variables' : ['u', 'p'],
        'is_linear' : True,
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
