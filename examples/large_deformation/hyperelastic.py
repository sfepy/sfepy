# -*- coding: utf-8 -*-
import numpy as nm

from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
    'post_process_hook' : 'stress_strain',
}


field_1 = {
    'name' : 'displacement',
    'dim' : (3,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

material_1 = {
    'name' : 'solid',
    'values' : {
        'K'  : 1e3, # bulk modulus
        'mu' : 20e0, # shear modulus of neoHookean term
        'kappa' : 10e0, # shear modulus of Mooney-Rivlin term
    }
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
# Dirichlet BC + related functions.
ebcs = {
    'l' : ('Left', {'u.all' : 0.0}),
    'r' : ('Right', {'u.0' : 0.0, 'u.[1,2]' : 'rotate_yz'}),
}

centre = nm.array( [0, 0], dtype = nm.float64 )

def rotate_yz(ts, coor, bc):
    from sfepy.base.la import rotation_matrix2d
    from sfepy.base.base import debug
    
    vec = coor[:,1:3] - centre

    angle = 10.0 * ts.step
    print 'angle:', angle

    mtx = rotation_matrix2d( angle )
    vec_rotated = nm.dot( vec, mtx )

    displacement = vec_rotated - vec
    
    return displacement.T.flat

functions = {
    'rotate_yz' : (rotate_yz,),
}

def stress_strain( out, problem, state, extend = False ):
    from sfepy.base.base import Struct, debug
    from sfepy.fem import eval_term_op as ev

    strain = ev( state, 'dw_tl_he_neohook.i1.Omega( solid.mu, v, u )',
                 problem, call_mode = 'de_strain' )
    out['green_strain'] = Struct( name = 'output_data',
                                  mode = 'cell', data = strain,
                                  dof_types = None )

    stress = ev( state, 'dw_tl_he_neohook.i1.Omega( solid.mu, v, u )',
                 problem, call_mode = 'de_stress' )
    out['neohook_stress'] = Struct( name = 'output_data',
                                    mode = 'cell', data = stress,
                                    dof_types = None )

    stress = ev( state, 'dw_tl_he_mooney_rivlin.i1.Omega( solid.kappa, v, u )',
                 problem, call_mode = 'de_stress' )
    out['mooney_rivlin_stress'] = Struct( name = 'output_data',
                                          mode = 'cell', data = stress,
                                          dof_types = None )

    stress = ev( state, 'dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )',
                 problem, call_mode = 'de_stress' )
    out['bulk_stress'] = Struct( name = 'output_data',
                                 mode = 'cell', data = stress,
                                 dof_types = None )

    return out

##
# Balance of forces.
integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'balance' : """dw_tl_he_neohook.i1.Omega( solid.mu, v, u )
                 + dw_tl_he_mooney_rivlin.i1.Omega( solid.kappa, v, u )
                 + dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )
                 = 0""",
}

##
# Solvers etc.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 5,
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
    'n_step' : 11, # has precedence over dt!
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 1000,
    'cache_override' : False,
}
