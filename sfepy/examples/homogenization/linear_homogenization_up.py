r"""
Compute homogenized elastic coefficients for a given  heterogeneous linear
elastic microstructure, see [1] for details or [2] and [3] for a quick
explanation. The mixed formulation, where displacements and pressures are
as unknowns, is used in this example.

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

import numpy as nm

import sfepy.discrete.fem.periodic as per
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson_mixed,\
    bulk_from_youngpoisson
from sfepy.homogenization.utils import define_box_regions, get_box_volume
import sfepy.homogenization.coefs_base as cb

from sfepy import data_dir
from sfepy.base.base import Struct
from sfepy.homogenization.recovery import compute_micro_u,\
    compute_stress_strain_u, compute_mac_stress_part, add_stress_p


def recovery_le(pb, corrs, macro):
    out = {}
    dim = corrs['corrs_le']['u_00'].shape[1]
    mic_u = - compute_micro_u(corrs['corrs_le'], macro['strain'], 'u', dim)
    mic_p = - compute_micro_u(corrs['corrs_le'], macro['strain'], 'p', dim)

    out['u_mic'] = Struct(name='output_data',
                          mode='vertex', data=mic_u)
    out['p_mic'] = Struct(name='output_data', mode='cell',
                          data=mic_p[:, nm.newaxis, :, nm.newaxis])

    stress_Y, strain_Y = \
        compute_stress_strain_u(pb, 'i', 'Y', 'mat.D', 'u', mic_u)
    stress_Y += \
        compute_mac_stress_part(pb, 'i', 'Y', 'mat.D', 'u', macro['strain'])
    add_stress_p(stress_Y, pb, 'i', 'Y', 'p', mic_p)

    strain = macro['strain'] + strain_Y

    out['cauchy_strain'] = Struct(name='output_data',
                                  mode='cell', data=strain)
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress_Y)

    return out


dim = 3
filename_mesh = data_dir + '/meshes/3d/matrix_fiber.mesh'
region_lbn = (0, 0, 0)
region_rtf = (1, 1, 1)


regions = {
    'Y': 'all',
    'Ym': 'cells of group 1',
    'Yc': 'cells of group 2',
}
regions.update(define_box_regions(dim, region_lbn, region_rtf))

materials = {
    'mat': ({'D': {'Ym': stiffness_from_youngpoisson_mixed(dim, 7.0e9, 0.4),
                   'Yc': stiffness_from_youngpoisson_mixed(dim, 70.0e9, 0.2)},
             'gamma': {'Ym': 1.0/bulk_from_youngpoisson(7.0e9, 0.4),
                       'Yc': 1.0/bulk_from_youngpoisson(70.0e9, 0.2)}},),
}

fields = {
    'corrector_u': ('real', dim, 'Y', 1),
    'corrector_p': ('real', 1, 'Y', 0),
}

variables = {
    'u': ('unknown field', 'corrector_u'),
    'v': ('test field', 'corrector_u', 'u'),
    'p': ('unknown field', 'corrector_p'),
    'q': ('test field', 'corrector_p', 'p'),
    'Pi': ('parameter field', 'corrector_u', 'u'),
    'Pi1u': ('parameter field', 'corrector_u', '(set-to-None)'),
    'Pi2u': ('parameter field', 'corrector_u', '(set-to-None)'),
    'Pi1p': ('parameter field', 'corrector_p', '(set-to-None)'),
    'Pi2p': ('parameter field', 'corrector_p', '(set-to-None)'),
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
    'volume': {'value': get_box_volume(dim, region_lbn, region_rtf), },
    'output_dir': 'output',
    'coefs_filename': 'coefs_le_up',
    'recovery_hook': 'recovery_le',
    'multiprocessing': False,
}

equation_corrs = {
    'balance_of_forces':
    """  dw_lin_elastic.i.Y(mat.D, v, u)
       - dw_stokes.i.Y(v, p) =
       - dw_lin_elastic.i.Y(mat.D, v, Pi)""",
    'pressure constraint':
    """- dw_stokes.i.Y(u, q)
       - dw_dot.i.Y(mat.gamma, q, p) =
       + dw_stokes.i.Y(Pi, q)""",
}

coefs = {
    'elastic_u': {
        'requires': ['pis', 'corrs_rs'],
        'expression': 'dw_lin_elastic.i.Y(mat.D, Pi1u, Pi2u)',
        'set_variables': [('Pi1u', ('pis', 'corrs_rs'), 'u'),
                          ('Pi2u', ('pis', 'corrs_rs'), 'u')],
        'class': cb.CoefSymSym,
    },
    'elastic_p': {
        'requires': ['corrs_rs'],
        'expression': 'dw_dot.i.Y(mat.gamma, Pi1p, Pi2p)',
        'set_variables': [('Pi1p', 'corrs_rs', 'p'),
                          ('Pi2p', 'corrs_rs', 'p')],
        'class': cb.CoefSymSym,
    },
    'D': {
        'requires': ['c.elastic_u', 'c.elastic_p'],
        'class': cb.CoefSum,
    },
    'filenames': {},
}

requirements = {
    'pis': {
        'variables': ['u'],
        'class': cb.ShapeDimDim,
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
    # 'ls': ('ls.auto_iterative', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e2,
    })
}
