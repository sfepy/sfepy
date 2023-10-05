r"""
An example demonstrating the usage of the truss elements in 2D.
The bridge structure is fixed on the left and supported on the right.

Running the simulation::

  sfepy-run sfepy/examples/linear_elasticity/truss_bridge.py

Viewing the results::

  sfepy-view bridge.vtk -f u:wu:p0 1:vw:p0 S:e:p1 --2d-view
"""
import numpy as nm
from sfepy.base.base import Struct
from sfepy import data_dir


def post_process(out, pb, state, extend=False):
    valS = pb.evaluate('ev_lin_truss_force.0.Truss(truss.EA, u)', mode='el_avg')
    out['S'] = Struct(name='output_data', mode='cell',
                      region_name='Truss', data=valS)

    print('### nodal displacements: ')
    print(out['u'].data)

    return out


def get_pars(ts, coors, mode=None, **kwargs):
    if mode == 'qp':
        EA = nm.array([2, 2, 2, 2, 2, 2,
                       10, 10, 10, 10, 10, 10,
                       3, 3, 3, 3, 3,
                       1, 1, 1, 1]).reshape(-1, 1, 1) * 1000

        return {'EA': EA}


filename_mesh = data_dir + "/meshes/2d/bridge.vtk"

options = {
    'post_process_hook': 'post_process',
}

regions = {
    'Truss': 'cells of group 1',
    'Left': ('vertices in (x < 0.01)', 'vertex'),
    'Right': ('vertices in (x > 59.99)', 'vertex'),
    'Bottom': ('vertex 2, 4, 6, 8, 10', 'vertex'),
}

materials = {
    'force': ({'.val': nm.array([[0, -10],
                                 [0, -10],
                                 [0, -16],
                                 [0, -10],
                                 [0, -10]])},),
    'truss': 'get_pars',
}

functions = {
    'get_pars': (get_pars,),
}

fields = {
    'displacement': ('real', 'vector', 'Truss', 1),
}

variables = {
    'u': ('unknown field', 'displacement'),
    'v': ('test field', 'displacement', 'u'),
}

ebcs = {
    'Fixed_Left': ('Left', {'u.all': 0.0}),
    'Support_Right': ('Right', {'u.1': 0.0}),
}

equations = {
    'balance_of_forces':
    """dw_lin_truss.0.Truss(truss.EA, v, u)
     = dw_point_load.0.Bottom(force.val, v)""",
}

solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {}),
}
