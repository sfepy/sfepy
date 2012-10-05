r"""
Example explaining how to change Dirichlet boundary conditions depending
on time. It is shown on the stationary Laplace equation for temperature,
so there is no dynamics, only the conditions change with time.

Five time steps are solved on a cube domain, with the temperature fixed
to zero on the bottom face, and set to other values on the left, right
and top faces in different time steps.

Find :math:`t` such that:

.. math::
    \int_{\Omega} c \nabla s \cdot \nabla t
    = 0
    \;, \quad \forall s \;.
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cube_medium_tetra.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < -0.499)', {}),
    'Right' : ('nodes in (x > 0.499)', {}),
    'Bottom' : ('nodes in (z < -0.499)', {}),
    'Top' : ('nodes in (z > 0.499)', {}),
}

materials = {
    'one' : ({'val' : 1.0},),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    't'   : ('unknown field', 'temperature', 0),
    's'   : ('test field',    'temperature', 't'),
}

ebcs = {
    'fixed' : ('Bottom', {'t.all' : 0}),
    't_t02' : ('Left', [(-0.5, 0.5), (2.5, 3.5)], {'t.all' : 1.0}),
    't_t1' : ('Right', [(0.5, 1.5)], {'t.all' : 2.0}),
    't_t4' : ('Top', 'is_ebc', {'t.all' : 3.0}),
}

def is_ebc(ts):
    if ts.step in (2, 4):
        return True

    else:
        return False

functions = {
    'is_ebc' : (is_ebc,),
}

integrals = {
    'ivol' : ('v', 2),
    'isurf' : ('s', 2),
}

equations = {
    'eq' : """dw_laplace.2.Omega( one.val, s, t ) = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
        'problem'   : 'nonlinear',
    }),
    'ts' : ('ts.simple', {
        't0'     : 0.0,
        't1'     : 4.0,
        'dt'     : None,
        'n_step' : 5, # has precedence over dt!

        'quasistatic' : True,
    }),
}
