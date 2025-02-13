# -*- coding: utf-8 -*-
r"""
Compressible Neo-Hookean hyperelastic material model. The tangent modulus
and the stress tensor are calculated by a user defined function.
"""
from functools import partial

import numpy as nm
from hyperelastic_ul import (filename_mesh, options, regions, fields,
    variables, ebcs, functions)

mu, K = 20., 1000.

# function which returns tangent modulus and stress tensor
def get_hyperelastic_mat(family_data, mode):
    from sfepy.terms.extmods.terms import sym2nonsym
    from sfepy.terms.terms_hyperelastic_ul import (NeoHookeanULTerm,
        BulkPenaltyULTerm)

    n_el, n_qp, sym, _= family_data.green_strain.shape
    dim = family_data.mtx_f.shape[-1]
    dim2 = dim**2

    ones = nm.ones((n_el, n_qp, 1, 1), dtype=nm.float64)
    mat_mu = ones * mu
    mat_K = ones * K

    fargs = [family_data.get(name)
             for name in NeoHookeanULTerm.family_data_names]
    stress = nm.empty((n_el, n_qp, sym, 1), dtype=nm.float64)
    NeoHookeanULTerm.stress_function(stress, mat_mu, *fargs)

    if mode == 'tan_mod':
        tanmod = nm.empty((n_el, n_qp, sym, sym), dtype=nm.float64)
        NeoHookeanULTerm.tan_mod_function(tanmod, mat_mu, *fargs)


    fargs = [family_data.get(name)
             for name in BulkPenaltyULTerm.family_data_names]
    stress_p = nm.empty((n_el, n_qp, sym, 1), dtype=nm.float64)
    BulkPenaltyULTerm.stress_function(stress_p, mat_K, *fargs)

    if mode == 'tan_mod':
        tanmod_p = nm.empty((n_el, n_qp, sym, sym), dtype=nm.float64)
        BulkPenaltyULTerm.tan_mod_function(tanmod_p, mat_K, *fargs)

        tanmod_ns = nm.zeros((n_el, n_qp, dim2, dim2), dtype=nm.float64)
        sym2nonsym(tanmod_ns, tanmod + tanmod_p)

        stress_ns = nm.zeros((n_el, n_qp, dim2, dim2), dtype=nm.float64)
        sym2nonsym(stress_ns, stress + stress_p)

        out = tanmod_ns + stress_ns

    elif mode == 'stress':
        out = stress + stress_p

    else:
        raise ValueError()

    return out


def stress_strain(out, problem, state, extend = False):
    from sfepy.base.base import Struct

    ev = partial(problem.evaluate, mode='el_avg',
                 get_hyperelastic_mat=get_hyperelastic_mat)
    stress = ev('dw_ul_he_by_fun.i.Omega(get_hyperelastic_mat, v, u)',
                term_mode='stress')
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress)
    strain = ev('dw_ul_he_by_fun.i.Omega(get_hyperelastic_mat, v, u)',
                term_mode='strain')
    out['green_strain'] = Struct(name='output_data',
                                 mode='cell', data=strain)
    return out


integrals = {
    'i': 2,
}

materials = {
}

equations = {
    'balance': 'dw_ul_he_by_fun.i.Omega(get_hyperelastic_mat, v, u) = 0',
}

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 20,
        'eps_a': 1e-5,
        'eps_r': 1e-3,
    }),
    'ts': ('ts.simple', {
        't0': 0,
        't1': 1,
        'n_step': 11,
    }),
}
