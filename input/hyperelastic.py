# -*- coding: utf-8 -*-
import numpy as nm

filename_mesh = 'database/simple.vtk'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
}


field_1 = {
    'name' : 'displacement',
    'dim' : (3,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

material_1 = {
    'name' : 'solid',
    'mode' : 'here',
    'region' : 'Omega',

    'K'  : 1e3, # bulk modulus
    'mu' : 20e0, # shear modulus
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

# Whole domain.
region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

# EBC regions.
region_1 = {
    'name' : 'Left',
    'select' : 'nodes in (x < 0.001)'
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.099)'
}

##
# Dirichlet BC + related functions.
ebcs = {
    'l' : ('Left', {'u.all' : 0.0}),
    'r' : ('Right', {'u.0' : 0.0, 'u.[1,2]' : 'rotate_yz'}),
}

centre = nm.array( [0, 0], dtype = nm.float64 )

def rotate_yz( bc, ts, coor ):
    from sfepy.base.la import rotation_matrix2d
    from sfepy.base.base import debug
    
    vec = coor[:,1:3] - centre

    angle = 10.0 * ts.step
    print 'angle:', angle

    mtx = rotation_matrix2d( angle )
    vec_rotated = nm.dot( vec, mtx )

    displacement = vec_rotated - vec
    
    return displacement.T.flat

##
# Balance of forces.
integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'balance' : """dw_tl_he_neohook.i1.Omega( solid.mu, v, u )
                 + dw_tl_bulk_penalty.i1.Omega( solid.K, v, u )
                 = 0""",
}

##
# Solvers etc.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
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
