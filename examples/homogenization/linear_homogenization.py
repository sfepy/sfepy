# 04.08.2009
#!
#! Homogenization: Linear Elasticity
#! =================================
#$ \centerline{Example input file, \today}

#! Homogenization of heterogeneous linear elastic material

import sfepy.fem.periodic as per
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.homogenization.utils import define_box_regions
import sfepy.homogenization.coefs_base as cb
from sfepy import data_dir
from sfepy.base.base import Struct
from sfepy.homogenization.recovery import compute_micro_u, compute_stress_strain_u, compute_mac_stress_part

def recovery_le( pb, corrs, macro ):

    out = {}

    dim = corrs['corrs_le']['u_00'].shape[1]
    mic_u = - compute_micro_u( corrs['corrs_le'], macro['strain'], 'u', dim )

    out['u_mic'] = Struct( name = 'output_data',
                           mode = 'vertex', data = mic_u,
                           var_name = 'u', dofs = None )

    stress_Y, strain_Y = compute_stress_strain_u( pb, 'i1', 'Y', 'mat.D', 'u', mic_u )
    stress_Y += compute_mac_stress_part( pb, 'i1', 'Y', 'mat.D', 'u', macro['strain'] )

    strain = macro['strain'] + strain_Y

    out['cauchy_strain'] = Struct( name = 'output_data',
                                   mode = 'cell', data = strain,
                                   dofs = None )
    out['cauchy_stress'] = Struct( name = 'output_data',
                                   mode = 'cell', data = stress_Y,
                                   dofs = None )
    return out

#! Mesh
#! ----
filename_mesh = data_dir + '/meshes/3d/matrix_fiber.mesh'
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
    'mat' : ({'D' : {'Ym': stiffness_from_youngpoisson(dim, 7.0e9, 0.4),
                     'Yc': stiffness_from_youngpoisson(dim, 70.0e9, 0.2)}},),
}
#! Fields
#! ------
#! Scalar field for corrector basis functions.
fields = {
    'corrector' : ('real', dim, 'Y', 1),
}
#! Variables
#! ---------
#! Unknown and corresponding test variables. Parameter fields
#! used for evaluation of homogenized coefficients.
variables = {
    'u'     : ('unknown field', 'corrector', 0),
    'v'     : ('test field', 'corrector', 'u'),
    'Pi'    : ('parameter field', 'corrector', 'u'),
    'Pi1' : ('parameter field', 'corrector', '(set-to-None)'),
    'Pi2' : ('parameter field', 'corrector', '(set-to-None)'),
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
    'i1' : ('v', 2),
    'i2' : ('s', 2),
}
#! Options
#! -------
#! Various problem-specific options.
options = {
    'coefs' : 'coefs',
    'requirements' : 'requirements',
    'ls' : 'ls', # linear solver to use
    'volume' : { 'variables' : ['u'],
                 'expression' : 'd_volume.i1.Y( u )' },
    'output_dir' : 'output',
    'coefs_filename' : 'coefs_le',
    'recovery_hook' : 'recovery_le',
}
#! Equations
#! ---------
#! Equations for corrector functions.
equation_corrs = {
    'balance_of_forces' :
    """dw_lin_elastic.i1.Y(mat.D, v, u ) =
     - dw_lin_elastic.i1.Y(mat.D, v, Pi )"""
}
#! Expressions for homogenized linear elastic coefficients.
expr_coefs = """dw_lin_elastic.i1.Y(mat.D, Pi1, Pi2 )"""
#! Coefficients
#! ------------
#! Definition of homogenized acoustic coefficients.
def set_elastic(variables, ir, ic, mode, pis, corrs_rs):
    mode2var = {'row' : 'Pi1', 'col' : 'Pi2'}

    val = pis.states[ir, ic]['u'] + corrs_rs.states[ir, ic]['u']

    variables[mode2var[mode]].set_data(val)

coefs = {
    'D' : {
        'requires' : ['pis', 'corrs_rs'],
        'expression' : expr_coefs,
        'set_variables' : set_elastic,
        'class' : cb.CoefSymSym,
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
