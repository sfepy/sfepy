r"""
Laplace equation using general term along with the short syntax of keywords.

See examples with dw_laplace term using short syntax
:ref:`diffusion-poisson_short_syntax`
or using long syntax :ref:`diffusion-poisson`.

Find :math:`t` such that:

.. math::
    \int_{\Omega} c \nabla s \cdot \nabla t
    = 0
    \;, \quad \forall s \;.
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

materials = {
    'coef' : ({'val' : 1.0},),
}

regions = {
    'Omega' : 'all', # or 'cells of group 6'
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.099999)', 'facet'),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integrals = {
    'i' : 2,
}

equations = {
    'Temperature' : """intFE.i.Omega.einsum(',i,i', coef.val, s.grad,
                                            t.grad) = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

options = {'post_process_hook_final': 'post_process_hook_final'}

def post_process_hook_final(pb, state):
    ev = pb.evaluate
    norms = []
    norms.append(ev('dw_laplace.i.Omega(coef.val, t, t)', mode='eval'))
    norms.append(ev("intFE.i.Omega.einsum(',i,i', coef.val, t.grad, t.grad)",
                    mode='eval'))
    print "Energetic norms for dw_laplace term (%s) and intFE term (%s)." \
        % (str(norms[0]), str(norms[1]))
