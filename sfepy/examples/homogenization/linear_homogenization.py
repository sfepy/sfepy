r"""
Compute homogenized elastic coefficients for a given heterogeneous linear
elastic microstructure.

See [1] for details or [2] and [3] for a quick explanation.

[1] D. Cioranescu,  J.S.J. Paulin: Homogenization in open sets with holes.
Journal of Mathematical Analysis and Applications 71(2), 1979, pages 590-607.
https://doi.org/10.1016/0022-247X(79)90211-7

[2] J. Pinho-da-Cruz, J.A. Oliveira, F. Teixeira-Dias:
Asymptotic homogenisation in linear elasticity.
Part I: Mathematical formulation and finite element modelling.
Computational Materials Science 45(4), 2009, pages 1073-1080.
http://dx.doi.org/10.1016/j.commatsci.2009.02.025

[3] J. Pinho-da-Cruz, J.A. Oliveira, F. Teixeira-Dias:
Asymptotic homogenisation in linear elasticity.
Part II: Finite element procedures and multiscale applications.
Computational Materials Science 45(4), 2009, pages 1081-1096.
http://dx.doi.org/10.1016/j.commatsci.2009.01.027
"""

import sfepy.discrete.fem.periodic as per
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.homogenization.utils import define_box_regions
import sfepy.homogenization.coefs_base as cb
from sfepy import data_dir
from sfepy.base.base import Struct
from sfepy.homogenization.recovery import compute_micro_u,\
    compute_stress_strain_u, compute_mac_stress_part


def recovery_le(pb, corrs, macro):

    out = {}

    dim = corrs['corrs_le']['u_00'].shape[1]
    mic_u = - compute_micro_u(corrs['corrs_le'], macro['strain'], 'u', dim)

    out['u_mic'] = Struct(name='output_data',
                          mode='vertex', data=mic_u)

    stress_Y, strain_Y = \
        compute_stress_strain_u(pb, 'i', 'Y', 'mat.D', 'u', mic_u)
    stress_Y += \
        compute_mac_stress_part(pb, 'i', 'Y', 'mat.D', 'u', macro['strain'])

    strain = macro['strain'] + strain_Y

    out['cauchy_strain'] = Struct(name='output_data',
                                  mode='cell', data=strain)
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress_Y)

    return out


filename_mesh = data_dir + '/meshes/3d/matrix_fiber.mesh'
dim = 3
region_lbn = (0, 0, 0)
region_rtf = (1, 1, 1)

regions = {
    'Y': 'all',
    'Ym': 'cells of group 1',
    'Yc': 'cells of group 2',
}
regions.update(define_box_regions(dim, region_lbn, region_rtf))

materials = {
    'mat': ({'D': {'Ym': stiffness_from_youngpoisson(dim, 7.0e9, 0.4),
                   'Yc': stiffness_from_youngpoisson(dim, 70.0e9, 0.2)}},),
}

fields = {
    'corrector': ('real', dim, 'Y', 1),
}

variables = {
    'u': ('unknown field', 'corrector', 0),
    'v': ('test field', 'corrector', 'u'),
    'Pi': ('parameter field', 'corrector', 'u'),
    'Pi1': ('parameter field', 'corrector', '(set-to-None)'),
    'Pi2': ('parameter field', 'corrector', '(set-to-None)'),
}

functions = {
    'match_x_plane': (per.match_x_plane,),
    'match_y_plane': (per.match_y_plane,),
    'match_z_plane': (per.match_z_plane,),
}

ebcs = {
    'fixed_u': ('Corners', {'u.all': 0.0}),
}

if dim == 3:
    epbcs = {
        'periodic_x': (['Left', 'Right'], {'u.all': 'u.all'},
                       'match_x_plane'),
        'periodic_y': (['Near', 'Far'], {'u.all': 'u.all'},
                       'match_y_plane'),
        'periodic_z': (['Top', 'Bottom'], {'u.all': 'u.all'},
                       'match_z_plane'),
    }
else:
    epbcs = {
        'periodic_x': (['Left', 'Right'], {'u.all': 'u.all'},
                       'match_x_plane'),
        'periodic_y': (['Bottom', 'Top'], {'u.all': 'u.all'},
                       'match_y_plane'),
    }

all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:dim]]

integrals = {
    'i': 2,
}

options = {
    'coefs': 'coefs',
    'requirements': 'requirements',
    'ls': 'ls',  # linear solver to use
    'volume': {'expression': 'ev_volume.i.Y(u)'},
    'output_dir': 'output',
    'coefs_filename': 'coefs_le',
    'recovery_hook': 'recovery_le',
}

equation_corrs = {
    'balance_of_forces':
    """dw_lin_elastic.i.Y(mat.D, v, u) =
     - dw_lin_elastic.i.Y(mat.D, v, Pi)"""
}

expr_coefs = """dw_lin_elastic.i.Y(mat.D, Pi1, Pi2)"""

coefs = {
    'D': {
        'requires': ['pis', 'corrs_rs'],
        'expression': expr_coefs,
        'set_variables': [('Pi1', ('pis', 'corrs_rs'), 'u'),
                          ('Pi2', ('pis', 'corrs_rs'), 'u')],
        'class': cb.CoefSymSym,
    },
    'filenames': {},
}

requirements = {
    'pis': {
        'variables': ['u'],
        'class': cb.ShapeDimDim,
        'save_name': 'corrs_pis',
    },
    'corrs_rs': {
        'requires': ['pis'],
        'ebcs': ['fixed_u'],
        'epbcs': all_periodic,
        'equations': equation_corrs,
        'set_variables': [('Pi', 'pis', 'u')],
        'class': cb.CorrDimDim,
        'save_name': 'corrs_le',
        'is_linear': True,
    },
}

solvers = {
    'ls': ('ls.auto_direct', {'use_presolve' : True}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-4,
    })
}
