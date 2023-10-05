r"""
An example demonstrating the usage of the truss structural elements in 3D.
The mixed (solid and structural elements) bridge structure is fixed
on the left and supported on the right.

Running the simulation::

  sfepy-run sfepy/examples/linear_elasticity/truss_bridge3d.py

Viewing the results::

  sfepy-view bridge3d_S*.vtk -f u_solid:s0:wu_solid:f1e3:p0 u_struct:s1:wu_struct:f1e3:p0
"""
import numpy as nm
from sfepy.base.base import Struct
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy import data_dir


def post_process(out, pb, state, extend=False):
    valP = pb.evaluate('ev_lin_truss_force.0.Struct(truss.EA, u)',
                       mode='el_avg')
    out['P'] = Struct(name='output_data', mode='cell',
                      region_name='Struct', data=valP)

    valS = pb.evaluate('ev_cauchy_stress.i.Solid(solid.D, u)', mode='el_avg')
    out['S'] = Struct(name='output_data', mode='cell',
                      region_name='Solid', data=valS)

    u = out.pop('u').data
    u_solid = u[pb.domain.regions['Solid'].vertices, :]
    out['u_solid'] = Struct(name='output_data', mode='vertex',
                            region_name='Solid', data=u_solid)
    u_truss = u[pb.domain.regions['Struct'].vertices, :]
    out['u_struct'] = Struct(name='output_data', mode='vertex',
                             region_name='Struct', data=u_truss)

    return out


filename_mesh = data_dir + "/meshes/3d/bridge3d.vtk"

options = {
    'post_process_hook': 'post_process',
    'split_results_by': 'region',
}

regions = {
    'Solid': 'cells of group 1',
    'Struct': 'cells of group 2',
    'Omega': ('r.Solid +v r.Struct', 'cell', None, {'finalize': False}),
    'Left': ('vertices in (x < 0.01)', 'vertex'),
    'Right': ('vertices in (x > 11.99)', 'vertex'),
    'Bottom': ('vertices in (z < -0.199)', 'vertex'),
    'Top': ('vertices in (z > -0.001)', 'facet'),
    'LeftBottom': ('r.Left *v r.Bottom', 'vertex'),
    'RightBottom': ('r.Right *v r.Bottom', 'vertex'),
}

materials = {
    'force': ({'val': nm.array([[0, 0, -1000.]]).T},),
    'truss': ({'EA': 210e9 * (0.02**2 * nm.pi)},),
    'solid': ({'D': stiffness_from_youngpoisson(3, 10e9, 0.3)},),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

variables = {
    'u': ('unknown field', 'displacement'),
    'v': ('test field', 'displacement', 'u'),
}

integrals = {
    'i': 2,
}

ebcs = {
    'Fixed_Left': ('LeftBottom', {'u.all': 0.0}),
    'Support_Right': ('RightBottom', {'u.2': 0.0}),
}

equations = {
    'balance_of_forces':
    """dw_lin_elastic.i.Solid(solid.D, v, u)
     + dw_lin_truss.0.Struct(truss.EA, v, u)
     = dw_surface_ltr.i.Top(force.val, v)""",
}

solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {'eps_a': 1e-6, 'eps_r': 1e-9}),
}
