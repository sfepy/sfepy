"""
Example usage::

  sfepy-run sfepy/examples/linear_elasticity/multi_point_constraints.py
  sfepy-run sfepy/examples/linear_elasticity/multi_point_constraints.py -d "dofs=(0,1)"
  sfepy-run sfepy/examples/linear_elasticity/multi_point_constraints.py -d "dofs=2"
  sfepy-view annulus-c.vtk -2
  sfepy-view annulus-c.vtk -2e -f u:wu:f1:p0 1:vw:p0 u:gu:p0
"""
import numpy as nm

from sfepy.mesh.mesh_generators import gen_cylinder_mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.discrete.fem.mesh import Mesh
from sfepy.linalg import get_coors_in_ball
from sfepy.mechanics.matcoefs import stiffness_from_lame

def define(dims=(1, 1, 2, 2, 0), shape=(3, 12, 0), order=1, dofs=(0, 1, 2),
           output_dir='.'):

    if isinstance(dofs, tuple):
        dofs = ','.join([f'{ii}' for ii in dofs])

    def mesh_hook(mesh, mode):
        if mode == 'read':
            _mesh = gen_cylinder_mesh(dims, shape, (0, 0, 0),
                                     make_2d=True)
            # Add constraint vertices (and cells) to the base mesh.
            coors, vgroups, conns, mat_ids, descs = _mesh._get_io_data()
            nc = coors.shape[0]
#            coors = nm.r_[coors, [[0, 0], [0.0, 0.001]]]
            coors = nm.r_[coors, [[0, 0], [0.001, 0.0]]]
            vgroups = nm.r_[vgroups, [1, 2]]
            conns.append([[nc, nc + 1]])
            mat_ids.append([1])
            descs.append('1_2')
            mesh = Mesh.from_data('annulus-c',
                                  coors, vgroups, conns, mat_ids, descs)
            return mesh

        elif mode == 'write':
            pass

    filename_mesh = UserMeshIO(mesh_hook)

    regions = {
        'Omega' : 'all',
        'Gamma1' : ('vertices by get_gamma1', 'facet'),
        'Gamma2' : ('vertices by get_gamma2', 'facet'),
        'OmegaC' : ('vertices by get_omegac', 'cell', None, {'cell_tdim': 1}),
        'Gamma1C' : ('vertices of group 1', 'vertex'),
        'Gamma2C' : ('vertices of group 2', 'vertex'),
    }

    centre = [0, 0]
    functions = {
        'get_omegac' : (lambda coors, domain:
                        get_coors_in_ball(coors, centre, 0.1),),
        'get_gamma1' : (lambda coors, domain:
                        get_coors_in_ball(coors, centre, 1.1),),
        'get_gamma2' : (lambda coors, domain:
                        get_coors_in_ball(coors, centre, 1.9, inside=False),),
    }

    fields = {
        'fu' : ('real', 2, 'Omega', order),
        'fc' : ('real', 3, 'OmegaC', 1),
    }

    variables = {
        'u' : ('unknown field', 'fu', 0),
        'v' : ('test field',    'fu', 'u'),
        'uc' : ('unknown field', 'fc', 1),
        'vc' : ('test field',    'fc', 'uc'),
    }

    ebcs = {
        'u0' : ('Gamma2', {'u.all' : 0.0}),
        'uc' : ('Gamma2C', {f'uc.[{dofs}]' : 0.1},),
    }

    lcbcs = {
        'rigid' : (('Gamma1', 'Gamma1C'), {'u.all' : 'uc.all'}, None, 'rigid2'),
    }

    materials = {
        'm' : ({
            'D' : stiffness_from_lame(dim=2, lam=1e1, mu=1e0),
            'ks' : [[1e+5], [1e+5], [1e+5]],
        },),
    }

    integrals = {
        'i' : 2 * order,
    }

    equations = {
        'eq1' :
        """dw_lin_elastic.i.Omega(m.D, v, u)
         + dw_lin_dspring_rot.0.OmegaC(m.ks, vc, uc)
         = 0
        """,
    }

    solvers = {
        'ls' : ('ls.auto_direct', {}),
        'newton' : ('nls.newton', {
            'i_max'      : 1,
            'eps_a'      : -1e-10,
            'check' : 0,
        }),
    }

    options = {
        'output_dir' : output_dir,
        'nls' : 'newton',
        'ls' : 'ls',
    }

    return locals()
