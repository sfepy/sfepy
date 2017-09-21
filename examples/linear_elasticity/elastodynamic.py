"""
"""
from __future__ import absolute_import

import numpy as nm

from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_lame

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < 0.001)', 'facet'),
    'Right' : ('vertices in (x > 0.099)', 'facet'),
}

materials = {
    'solid' : ({'D': stiffness_from_lame(dim=3, lam=1e1, mu=1e0),
                'rho': 10.0},),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i' : 1,
}

variables = {
    'u' : ('unknown field', 'displacement', 0, 1),
    'du' : ('unknown field', 'displacement', 1, 1),
    'ddu' : ('unknown field', 'displacement', 2, 1),
    # 'du' : ('parameter field', 'displacement', {'ic' : 'get_ic_du'}, 1),
    # 'ddu' : ('parameter field', 'displacement', {'ic' : 'get_ic_ddu'}, 1),
    'v' : ('test field', 'displacement', 'u'),
    'dv' : ('test field', 'displacement', 'du'),
    'ddv' : ('test field', 'displacement', 'ddu'),
}

ebcs = {
    'Fixed' : ('Left', {'u.all' : 0.0, 'du.all' : 0.0, 'ddu.all' : 0.0}),
    'Displaced' : ('Right', {'u.0' : 0.01, 'u.[1,2]' : 0.0,
                             'du.all' : 0.0, 'ddu.all' : 0.0}),
}

def get_ic(coor, ic, mode='u'):
    val = nm.zeros_like(coor)
    if mode == 'u':
        val[:, 0] = 0.0

    elif mode == 'du':
        val[:, 1] = 0.1

    return val

functions = {
    'get_ic_u' : (get_ic,),
    'get_ic_du' : (lambda coor, ic: get_ic(coor, None, mode='du'),),
}

ics = {
    'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
}

equations = {
    'balance_of_forces' :
    """dw_volume_dot.i.Omega(solid.rho, ddv, ddu)
     + dw_zero.i.Omega(dv, du)
     + dw_lin_elastic.i.Omega(solid.D, v, u) = 0""",
    # 'balance_of_forces' :
    # """dw_volume_dot.i.Omega(solid.rho, v, d2u/dt2)
    #  + dw_lin_elastic.i.Omega(solid.D, v, u) = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
    # 'tss' : ('ts.simple', {
    #     't0' : 0.0,
    #     't1' : 1.0,
    #     'dt' : 0.1,
    #     'n_step' : None,
    # }),
    # 'tsb' : ('ts.bathe', {
    #     't0' : 0.0,
    #     't1' : 1.0,
    #     'dt' : 0.1,
    #     'n_step' : None,
    # }),
    'tsn' : ('ts.newmark', {
        't0' : 0.0,
        't1' : 1.0,
        'dt' : 0.1,
        'n_step' : None,

        'u' : 'u',
        'v' : 'du',
        'a' : 'ddu',
    }),
    # 'tsp' : ('ts.petsc', {
    #     'method' : 'rk',
    #     't0' : 0.0,
    #     't1' : 1.0,
    #     'dt' : 0.1,
    #     'n_step' : None,
    # }),
}
