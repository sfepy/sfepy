"""
Contact of two 2D elastic bodies.
"""
from sfepy.mechanics.matcoefs import stiffness_from_lame
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/two_bodies.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

fields = {
    'displacement': ('real', 2, 'Omega', 1),
}

materials = {
    'solid' : ({'D': stiffness_from_lame(2, lam=5.769, mu=3.846)},),
    'load' : ({'val' : 1.0},)
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

regions = {
    'Omega' : 'all',
    'Omega0' : 'cells of group 0',
    'Omega1' : 'cells of group 1',
    'Bottom' : ('vertices in (y < -0.499)', 'facet'),
    'Top' : ('vertices in (y > 0.249)', 'facet'),
    'Contact0' : ('vertices in ((y > -0.0001) & (y < 0.0001))',
                  'facet', 'Omega0'),
    'Contact1' : ('vertices in ((y > -0.0001) & (y < 0.0001))',
                  'facet', 'Omega1'),
}

ebcs = {
    'fixb' : ('Bottom', {'u.all' : 0.0}),
}

equations = {
    'elasticity' :
    """dw_lin_elastic.2.Omega(solid.D, v, u)
     + dw_contact.10.Contact0(v, u)
     - dw_contact.10.Contact1(v, u)
    = - dw_surface_ltr.2.Top(load.val, v)""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
            'i_max' : 5,
            'eps_a' : 1e-10,
            'eps_r' : 1.0,
            'macheps' : 1e-16,
            # Linear system error < (eps_a * lin_red).
            'lin_red' : 1e-2,
            'ls_red' : 0.1,
            'ls_red_warp' : 0.001,
            'ls_on' : 1.1,
            'ls_min' : 1e-5,
            'check' : 0,
            'delta' : 1e-6,
    })
}
