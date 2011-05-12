"""
Example comparing new and old-style terms with scalar variables.
"""
import os

import numpy as nm

from sfepy import data_dir
from sfepy.fem import MeshIO

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'
#filename_mesh = data_dir + '/meshes/3d/cube_big_tetra.mesh'

conf_dir = os.path.dirname(__file__)
io = MeshIO.any_from_filename(filename_mesh, prefix_dir=conf_dir)
bbox = io.read_bounding_box()

dd = bbox[1] - bbox[0]

xmin, ymin, zmin = bbox[0, :] + 1e-4 * dd
xmax, ymax, zmax = bbox[1, :] - 1e-4 * dd

def post_process(out, pb, state, extend=False):

    for vn in ['p', 'r']:
        try:
            dd = pb.evaluate('dw_new_diffusion.2.Omega(m.c, %s, %s)'
                             % (vn, vn))
            print dd

            mass = pb.evaluate('dw_new_mass_scalar.2.Omega(%s, %s)'
                               % (vn, vn))
            print mass

        except:
            pass

    return out

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'post_process_hook' : 'post_process',
}

materials = {
    'm' : ({'c' : 0.0001 * nm.eye(3)},),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_Left' : ('nodes in (x < %f)' % xmin, {}),
    'Gamma_Right' : ('nodes in (x > %f)' % xmax, {}),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 2),
}

variables = {
    'p' : ('unknown field', 'temperature', 0),
    'q' : ('test field',    'temperature', 'p'),
    'r' : ('unknown field', 'temperature', 1),
    's' : ('test field',    'temperature', 'r'),
}

ebcs = {
    'p1' : ('Gamma_Left', {'p.0' : 2.0}),
    'p2' : ('Gamma_Right', {'p.0' : -2.0}),
    'r1' : ('Gamma_Left', {'r.0' : 2.0}),
    'r2' : ('Gamma_Right', {'r.0' : -2.0}),
}

equations = {
    'new equation' :
    """dw_new_diffusion.2.Omega(m.c, q, p)
     + dw_new_mass_scalar.2.Omega(q, p)
     = 0""",
    'equation' :
    """dw_diffusion.2.Omega(m.c, s, r)
     + dw_mass_scalar.2.Omega(s, r)
     = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
    }),
}
