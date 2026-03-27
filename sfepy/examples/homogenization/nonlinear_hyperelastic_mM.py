"""
Homogenized nonlinear hyperelastic material with evolving microstructure
deformation in each macroscopic quadrature point.
"""
import numpy as nm
from functools import partial
from sfepy import data_dir, base_dir
from sfepy.base.base import Struct, output
from sfepy.homogenization.micmac import (get_homog_coefs_nonlinear,
    get_homogen_app_form_cache)
import sfepy.linalg as la
from sfepy.discrete.evaluate import Evaluator

hyperelastic_data = {'coefs': {}}

def post_process(out, pb, state, extend=False):
    ev = partial(pb.evaluate, mode='el_avg',
                 get_homog_mat=pb.conf.get_homog_mat)
    stress = ev('dw_ul_he_by_fun.1.Omega(get_homog_mat, v, u)',
                term_mode='stress')
    strain = ev('dw_ul_he_by_fun.1.Omega(get_homog_mat, v, u)',
                term_mode='strain')

    out['cauchy_stress'] = Struct(name='output_data', mode='cell',
                                  data=stress)
    out['green_strain'] = Struct(name='output_data', mode='cell',
                                 data=strain)

    if pb.conf.options.get('recover_micro', False):
        happ = get_homogen_app_form_cache(pb.conf.options['micro_filename'])
        if pb.ts.step == 0:
            rname = pb.conf.options.recovery_region
            rcells = pb.domain.regions[rname].get_cells()
            happ.app_options.store_micro_idxs = strain.shape[1] * rcells
        else:
            hpb = happ.problem
            recovery_hook = hpb.conf.options.get('recovery_hook', None)
            if recovery_hook is not None:
                recovery_hook = hpb.conf.get_function(recovery_hook)
                rname = pb.conf.options.recovery_region
                rcoors = []
                for ii in happ.app_options.store_micro_idxs:
                    key = happ.get_micro_cache_key('coors', ii, pb.ts.step)
                    if key in happ.micro_state_cache:
                        rcoors.append(happ.micro_state_cache[key])

                recovery_hook(hpb, rcoors, pb.domain.regions[rname], pb.ts)

    return out


def get_homog_mat(family_data, mode):
    pb = hyperelastic_data['problem']
    ts = pb.get_timestepper()
    output(f'macro mat. fun: step={ts.step}, iiter={pb.iiter}')

    ckey = (ts.step, pb.iiter)
    ccache = hyperelastic_data['coefs']
    n_el, n_qp, dim, _ = family_data.mtx_f.shape
    sym = family_data.green_strain.shape[2]
    dim2 = dim**2

    if ckey not in ccache:
        ccache.clear()

        mtx_f = family_data.mtx_f.reshape((n_el * n_qp, dim, dim))
        if hasattr(pb, 'mtx_f_prev'):
            rel_mtx_f = la.dot_sequences(mtx_f, nm.linalg.inv(pb.mtx_f_prev),
                                         'AB')
        else:
            rel_mtx_f = mtx_f

        pb.mtx_f_prev = mtx_f.copy()

        macro_data = {'mtx_e': rel_mtx_f - nm.eye(dim)}
        ccache[ckey] = get_homog_coefs_nonlinear(ts, mtx_f, 'qp', macro_data,
                                                 problem=pb,
                                                 iteration=pb.iiter)

    coefs = ccache[ckey]

    if mode == 'tan_mod':
        out = coefs['A'].reshape((n_el, n_qp, dim2, dim2))
    elif mode == 'stress':
        out = coefs['S'].reshape((n_el, n_qp, sym, 1))
    else:
        raise ValueError()

    return out


def ulf_iteration_hook(pb, nls, vec, it, err, err0):
    Evaluator.new_ulf_iteration(pb, nls, vec, it, err, err0)
    pb.iiter = it


def ulf_init(pb):
    pb.domain.mesh.coors_act = pb.domain.mesh.coors.copy()
    pb.iiter = 0
    hyperelastic_data['problem'] = pb


options = {
    'output_dir': 'output',
    'mesh_update_variables': ['u'],
    'nls_iter_hook': ulf_iteration_hook,
    'pre_process_hook': ulf_init,
    'micro_filename': (base_dir +
                       '/examples/homogenization/nonlinear_homogenization.py'),
    'post_process_hook': post_process,
    'recover_micro': True,
    'recovery_region': 'Recovery',
}

materials = {}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

variables = {
    'u': ('unknown field', 'displacement'),
    'v': ('test field', 'displacement', 'u'),
}

filename_mesh = data_dir + '/meshes/2d/its2D.mesh'

regions = {
    'Omega': 'all',
    'Left': ('vertices in (x < 0.001)', 'facet'),
    'Bottom': ('vertices in (y < 0.001 )', 'facet'),
    'Recovery': ('cell 49, 81', 'cell'),
}

ebcs = {
    'l': ('Left', {'u.all': 0.0}),
    'b': ('Bottom', {'u.all': 'move_bottom'}),
}


centre = nm.array([0, 0], dtype=nm.float64)


def move_bottom(ts, coor, **kwargs):
    from sfepy.linalg import rotation_matrix2d

    vec = coor[:, 0:2] - centre
    angle = 3 * ts.step
    print('angle:', angle)
    mtx = rotation_matrix2d(angle)
    out = nm.dot(vec, mtx) - vec

    return out


functions = {
    'move_bottom': (move_bottom,),
}

equations = {
    'balance_of_forces': 'dw_ul_he_by_fun.1.Omega(get_homog_mat, v, u) = 0'
}

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'eps_a': 1e-3,
        'eps_r': 1e-3,
        'i_max': 20,
    }),
    'ts': ('ts.simple', {
        't0': 0,
        't1': 1,
        'n_step': 3 + 1,
        'verbose': 1,
    })
}
