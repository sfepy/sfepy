r"""
Nonlinear diffusion with a field-dependent coefficient and parametric sweep.
Solve the problem:
.. math::
   -\nabla \cdot \left( (1 + \alpha u^2)\nabla u \right)
   = \sin(\pi x)\sin(\pi y)
   \quad \text{in } \Omega,
with homogeneous Dirichlet boundary conditions:
.. math::
   u = 0 \quad \text{on } \partial \Omega.
The diffusion coefficient depends on the solution:
.. math::
   c(u) = 1 + \alpha u^2.
The problem is solved for several values of :math:`\alpha`, and the
solutions are compared to the linear case :math:`\alpha = 0`.
"""
from sfepy import data_dir
from sfepy.base.base import output
import numpy as nm

filename_mesh = data_dir + '/meshes/2d/square_unit_tri.mesh'

alphas = [0.0, 0.25, 0.5, 1.0, 2.0]

# shared state for the material function callback
current_alpha = 0.0

def conductivity(u):
    val = 1.0 + current_alpha * u**2
   #output('conductivity: min:', val.min(), 'max:', val.max())
    return val

def d_conductivity(u):
    # derivative of (1 + alpha * u^2) w.r.t. u
    return 2.0 * current_alpha * u

def get_rhs(ts, coors, mode=None, **kwargs):
    if mode == 'qp':
        x = coors[:, 0]
        y = coors[:, 1]
        val = nm.sin(nm.pi * x) * nm.sin(nm.pi * y)
        val = val.reshape((-1, 1, 1))
        return {'val': val}

def vary_alpha(problem):
    global current_alpha
    baseline = None
    # Newton handles the nonlinearity; we just sweep alpha values
    for alpha in alphas:
        current_alpha = alpha
        yield problem, []
        vec = problem.get_variables()['u']().copy()
        if baseline is None:
            baseline = vec.copy()
        diff_norm = nm.linalg.norm(vec - baseline)
        sol_norm = nm.linalg.norm(vec)
        output('---')
        output('alpha:', alpha, '||u||:', sol_norm, '||u - u0||:', diff_norm)
        yield

materials = {
    'rhs': 'get_rhs',
}

fields = {
    'fu': ('real', 1, 'Omega', 1),
}

variables = {
    'u': ('unknown field', 'fu', 0),
    'v': ('test field', 'fu', 'u'),
}

regions = {
    'Omega': 'all',
    'Gamma': ('vertices of surface', 'facet'),
}

ebcs = {
    'u0': ('Gamma', {'u.0': 0.0}),
}

functions = {
    'get_rhs': (get_rhs,),
}

integrals = {
    'i': 1,
}

equations = {
    'Balance': """
        dw_nl_diffusion.i.Omega( conductivity, d_conductivity, v, u )
        = dw_volume_lvf.i.Omega( rhs.val, v )
    """
}

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 20,
        'eps_a': 1e-10,
        'eps_r': 1.0,
    }),
}

options = {
    'nls': 'newton',
    'ls': 'ls',
    'parametric_hook': 'vary_alpha',
}