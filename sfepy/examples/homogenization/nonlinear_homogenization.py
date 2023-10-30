# -*- coding: utf-8 -*-
import numpy as nm
from sfepy.homogenization.utils import define_box_regions
import sfepy.homogenization.coefs_base as cb
import sfepy.discrete.fem.periodic as per
from sfepy.base.base import Struct
from sfepy.terms.terms_hyperelastic_ul import\
    HyperElasticULFamilyData, NeoHookeanULTerm, BulkPenaltyULTerm
from sfepy.terms.extmods.terms import sym2nonsym
from sfepy.discrete.functions import ConstantFunctionByRegion
from sfepy import data_dir
import sfepy.linalg as la


def recovery_hook(pb, ncoors, region, ts,
                  naming_scheme='step_iel', recovery_file_tag=''):
    from sfepy.base.ioutils import get_print_info
    from sfepy.homogenization.recovery import get_output_suffix
    import os.path as op

    for ii, icell in enumerate(region.cells):
        out = {}
        pb.set_mesh_coors(ncoors[ii], update_fields=True,
                          clear_all=False, actual=True)
        stress = pb.evaluate('ev_integrate_mat.3.Y(mat_he.S, u)',
                             mode='el_avg')

        out['cauchy_stress'] = Struct(name='output_data',
                                      mode='cell',
                                      data=stress,
                                      dofs=None)

        strain = pb.evaluate('ev_integrate_mat.3.Y(mat_he.E, u)',
                             mode='el_avg')

        out['green_strain'] = Struct(name='output_data',
                                     mode='cell',
                                     data=strain,
                                     dofs=None)

        out['displacement'] = Struct(name='output_data',
                                     mode='vertex',
                                     data=ncoors[ii] - pb.get_mesh_coors(),
                                     dofs=None)

        output_dir = pb.conf.options.get('output_dir', '.')
        format = get_print_info(pb.domain.mesh.n_el, fill='0')[1]
        suffix = get_output_suffix(icell, ts, naming_scheme, format,
                                   pb.output_format)

        micro_name = pb.get_output_name(extra='recovered_'
                                        + recovery_file_tag + suffix)
        filename = op.join(output_dir, op.basename(micro_name))
        fpv = pb.conf.options.get('split_results_by', None)
        pb.save_state(filename, out=out, split_results_by=fpv)


def def_mat(ts, mode, coors, term, pb):
    if not (mode == 'qp'):
        return

    if not hasattr(pb, 'family_data'):
        pb.family_data = HyperElasticULFamilyData()

    update_var = pb.conf.options.mesh_update_variable
    if pb.equations is None:
        state_u = pb.create_variables([update_var])[update_var]
    else:
        state_u = pb.get_variables()[update_var]

    if state_u.data[0] is None:
        state_u.init_data()

    state_u.set_data(
        pb.domain.get_mesh_coors(actual=True) - pb.domain.get_mesh_coors())
    state_u.field.clear_mappings()
    family_data = pb.family_data(state_u, term.region, term.integral,
                                 list(term.geometry_types.values())[0])

    if len(state_u.field.mappings0) == 0:
        state_u.field.save_mappings()

    n_el, n_qp, dim, n_en, n_c = state_u.get_data_shape(term.integral,
                                                        term.act_integration,
                                                        term.region.name)

    conf_mat = pb.conf.materials
    solid_key = [key for key in conf_mat.keys() if 'solid' in key][0]
    solid_mat = conf_mat[solid_key].values
    mat = {}
    for mat_key in ['mu', 'K']:
        mat_fun = ConstantFunctionByRegion({mat_key: solid_mat[mat_key]})
        mat[mat_key] = mat_fun.function(ts=ts, coors=coors, mode='qp',
            term=term, problem=pb)[mat_key].reshape((n_el, n_qp, 1, 1))

    shape = family_data.green_strain.shape[:2]
    sym = family_data.green_strain.shape[-2]
    dim2 = dim**2

    fargs = [family_data.get(name)
             for name in NeoHookeanULTerm.family_data_names]
    stress = nm.empty(shape + (sym, 1), dtype=nm.float64)
    tanmod = nm.empty(shape + (sym, sym), dtype=nm.float64)
    NeoHookeanULTerm.stress_function(stress, mat['mu'], *fargs)
    NeoHookeanULTerm.tan_mod_function(tanmod, mat['mu'], *fargs)

    fargs = [family_data.get(name)
             for name in BulkPenaltyULTerm.family_data_names]
    stress_p = nm.empty(shape + (sym, 1), dtype=nm.float64)
    tanmod_p = nm.empty(shape + (sym, sym), dtype=nm.float64)
    BulkPenaltyULTerm.stress_function(stress_p, mat['K'], *fargs)
    BulkPenaltyULTerm.tan_mod_function(tanmod_p, mat['K'], *fargs)

    stress_ns = nm.zeros(shape + (dim2, dim2), dtype=nm.float64)
    tanmod_ns = nm.zeros(shape + (dim2, dim2), dtype=nm.float64)
    sym2nonsym(stress_ns, stress + stress_p)
    sym2nonsym(tanmod_ns, tanmod + tanmod_p)

    npts = nm.prod(shape)
    J = family_data.det_f
    mtx_f = family_data.mtx_f.reshape((npts, dim, dim))

    out = {
        'E': 0.5 * (la.dot_sequences(mtx_f, mtx_f, 'ATB') - nm.eye(dim)),
        'A': ((tanmod_ns + stress_ns) / J).reshape((npts, dim2, dim2)),
        'S': ((stress + stress_p) / J).reshape((npts, sym, 1)),
    }

    return out


filename_mesh = data_dir + '/meshes/2d/special/circle_in_square_small.mesh'
dim = 2

options = {
    'coefs': 'coefs',
    'requirements': 'requirements',
    'volume': {'expression': 'ev_volume.5.Y(u)'},
    'output_dir': './output',
    'coefs_filename': 'coefs_hyper_homog',
    'multiprocessing': True,
    'chunks_per_worker': 2,
    'micro_update': {'coors': [('corrs_rs', 'u', 'mtx_e')]},
    'mesh_update_variable': 'u',
    'recovery_hook': 'recovery_hook',
    'store_micro_idxs': [49, 81],
}

fields = {
    'displacement': ('real', 'vector', 'Y', 1),
}

functions = {
    'match_x_plane': (per.match_x_plane,),
    'match_y_plane': (per.match_y_plane,),
    'mat_fce': (lambda ts, coors, mode=None, term=None, problem=None, **kwargs:
                def_mat(ts, mode, coors, term, problem),),
}

materials = {
    'mat_he': 'mat_fce',
    'solid': ({'K': {'Ym': 1000, 'Yc': 1000},
               'mu': {'Ym': 100, 'Yc': 10},
               },),
}

variables = {
    'u': ('unknown field', 'displacement'),
    'v': ('test field', 'displacement', 'u'),
    'Pi': ('parameter field', 'displacement', 'u'),
    'Pi1u': ('parameter field', 'displacement', '(set-to-None)'),
    'Pi2u': ('parameter field', 'displacement', '(set-to-None)'),
}

regions = {
    'Y': 'all',
    'Ym': 'cells of group 1',
    'Yc': 'cells of group 2',
}

regions.update(define_box_regions(dim, (0., 0.), (1., 1.)))

ebcs = {
    'fixed_u': ('Corners', {'u.all': 0.0}),
}

epbcs = {
    'periodic_ux': (['Left', 'Right'], {'u.all': 'u.all'}, 'match_x_plane'),
    'periodic_uy': (['Bottom', 'Top'], {'u.all': 'u.all'}, 'match_y_plane'),
}

coefs = {
    'A': {
        'requires': ['pis', 'corrs_rs'],
        'expression': 'dw_nonsym_elastic.3.Y(mat_he.A, Pi1u, Pi2u)',
        'set_variables': [('Pi1u', ('pis', 'corrs_rs'), 'u'),
                          ('Pi2u', ('pis', 'corrs_rs'), 'u')],
        'class': cb.CoefNonSymNonSym,
    },
    'S': {
        'expression': 'ev_integrate_mat.3.Y(mat_he.S, u)',
        'class': cb.CoefOne,
    }
}

requirements = {
    'pis': {
        'variables': ['u'],
        'class': cb.ShapeDimDim,
    },
    'corrs_rs': {
        'requires': ['pis'],
        'ebcs': ['fixed_u'],
        'epbcs': ['periodic_ux', 'periodic_uy'],
        'equations': {
            'balance_of_forces':
                """dw_nonsym_elastic.3.Y(mat_he.A, v, u)
               = - dw_nonsym_elastic.3.Y(mat_he.A, v, Pi)"""
        },
        'set_variables': [('Pi', 'pis', 'u')],
        'class': cb.CorrDimDim,
        'save_name': 'corrs_hyper_homog',
    },
}

solvers = {
    'ls': ('ls.auto_direct', {'use_presolve' : True}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-4,
        'problem': 'nonlinear',
    }),
}
