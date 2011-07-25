# -*- coding: utf-8 -*-
import numpy as nm

from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

vf_matrix = 0.5
vf_fibres1 = 0.2
vf_fibres2 = 0.3

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
    'post_process_hook' : 'stress_strain',
}


fields = {
    'displacement': (nm.float64, 3, 'Omega', 1),
}

materials = {
    'solid' : ({
        'K'  : vf_matrix * 1e3, # bulk modulus
        'mu' : vf_matrix * 20e0, # shear modulus of neoHookean term
    },),
    'f1' : 'get_pars_fibres1',
    'f2' : 'get_pars_fibres2',
}

def get_pars_fibres(ts, coors, mode=None, which=0, vf=1.0, **kwargs):
    """
    Parameters
    ----------
    ts : TimeStepper
        Time stepping info.
    coors : array_like
        The physical domain coordinates where the parameters shound be defined.
    mode : 'qp' or 'special'
        Call mode.
    which : int
        Fibre system id.
    vf : float
        Fibre system volume fraction.
    """
    if mode != 'qp': return

    fmax = 10.0
    eps_opt = 0.01
    s = 1.0

    tt = ts.nt * 2.0 * nm.pi

    if which == 0: # system 1
        fdir = nm.array([1.0, 0.0, 0.0], dtype=nm.float64)
        act = 0.5 * (1.0 + nm.sin(tt - (0.5 * nm.pi)))

    elif which == 1: # system 2
        fdir = nm.array([0.0, 1.0, 0.0], dtype=nm.float64)
        act = 0.5 * (1.0 + nm.sin(tt + (0.5 * nm.pi)))

    else:
        raise ValueError('unknown fibre system! (%d)' % which)

    fdir.shape = (3, 1)
    fdir /= nm.linalg.norm(fdir)

    print act

    shape = (coors.shape[0], 1, 1)
    out = {
        'fmax' : vf * nm.tile(fmax, shape),
        'eps_opt' : nm.tile(eps_opt, shape),
        's' : nm.tile(s, shape),
        'fdir' : nm.tile(fdir, shape),
        'act' : nm.tile(act, shape),
    }

    return out

functions = {
    'get_pars_fibres1' : (lambda ts, coors, mode=None, **kwargs:
                          get_pars_fibres(ts, coors, mode=mode, which=0,
                                          vf=vf_fibres1, **kwargs),),
    'get_pars_fibres2' : (lambda ts, coors, mode=None, **kwargs:
                          get_pars_fibres(ts, coors, mode=mode, which=1,
                                          vf=vf_fibres2, **kwargs),),
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001)', {}),
    'Right' : ('nodes in (x > 0.099)', {}),
}

##
# Dirichlet BC.
ebcs = {
    'l' : ('Left', {'u.all' : 0.0}),
}

##
# Balance of forces.
integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'balance'
        : """dw_tl_he_neohook.i1.Omega( solid.mu, v, u )
           + dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )
           + dw_tl_fib_a.i1.Omega( f1.fmax, f1.eps_opt, f1.s, f1.fdir, f1.act,
                                   v, u )
           + dw_tl_fib_a.i1.Omega( f2.fmax, f2.eps_opt, f2.s, f2.fdir, f2.act,
                                   v, u )
           = 0""",
}

def stress_strain(out, problem, state, extend=False):
    from sfepy.base.base import Struct, debug

    ev = problem.evaluate
    strain = ev('dw_tl_he_neohook.i1.Omega( solid.mu, v, u )',
                mode='el_avg', term_mode='strain')
    out['green_strain'] = Struct(name='output_data',
                                 mode='cell', data=strain, dofs=None)

    stress = ev('dw_tl_he_neohook.i1.Omega( solid.mu, v, u )',
                mode='el_avg', term_mode='stress')
    out['neohook_stress'] = Struct(name='output_data',
                                   mode='cell', data=stress, dofs=None )

    stress = ev('dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )',
                mode='el_avg', term_mode= 'stress')
    out['bulk_stress'] = Struct(name='output_data',
                                mode='cell', data=stress, dofs=None)

    return out

##
# Solvers etc.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 7,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'    : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp': 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'      : 0,
    'delta'      : 1e-6,
    'is_plot'    : False,
    'problem'    : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

solver_2 = {
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'    : 0,
    't1'    : 1,
    'dt'    : None,
    'n_step' : 21, # has precedence over dt!
}
