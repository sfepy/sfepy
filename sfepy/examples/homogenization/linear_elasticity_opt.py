import numpy as nm

from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy import data_dir

def mesh_hook(mesh, mode):
    if mode == 'read':
        mesh = gen_block_mesh([0.0098, 0.0011, 0.1], [5, 3, 17],
                              [0, 0, 0.05], name='specimen',
                              verbose=False)
        return mesh

    elif mode == 'write':
        pass

def optimization_hook(pb):
    cnf = pb.conf
    out = []
    yield pb, out

    state = out[-1][1].get_state_parts()
    coors = pb.domain.cmesh.coors
    displ = state['u'].reshape((coors.shape[0],3))
    # elongation
    mcoors = coors[cnf.mnodes, 2]
    mdispl = displ[cnf.mnodes, 2]
    dl = (mdispl[1] - mdispl[0]) / (mcoors[1] - mcoors[0])

    if hasattr(cnf, 'opt_data'):
        # compute slope of the force-elongation curve
        cnf.opt_data['k'] = cnf.F / dl

    yield None

def get_mat(coors, mode, pb):
    if mode == 'qp':
        # get material data
        if hasattr(pb.conf, 'opt_data'):
            # from homogenization
            D = pb.conf.opt_data['D_homog']
        else:
            # given values
            D = stiffness_from_youngpoisson(3, 150.0e9, 0.3)

        nqp = coors.shape[0]
        return {'D': nm.tile(D, (nqp, 1, 1))}

def define(is_opt=False):
    filename_mesh = UserMeshIO(mesh_hook)
    mnodes = (107, 113) # nodes for elongation eval.

    regions = {
        'Omega': 'all',
        'Bottom': ('vertices in (z < 0.001)', 'facet'),
        'Top': ('vertices in (z > 0.099)', 'facet'),
    }

    functions = {
        'get_mat': (lambda ts, coors, mode=None, problem=None, **kwargs:
                    get_mat(coors, mode, problem),),
    }

    S = 1.083500e-05    # specimen cross-section
    F = 5.0e3           # force
    materials = {
        'solid': 'get_mat',
        'load': ({'val': F / S},),
    }

    fields = {
        'displacement': ('real', 'vector', 'Omega', 1),
    }

    variables = {
        'u': ('unknown field', 'displacement', 0),
        'v': ('test field', 'displacement', 'u'),
    }

    ebcs = {
        'FixedBottom': ('Bottom', {'u.all': 0.0}),
        'FixedTop': ('Top', {'u.0': 0.0, 'u.1': 0.0}),
    }

    equations = {
        'balance_of_forces' :
        """dw_lin_elastic.5.Omega(solid.D, v, u)
         = dw_surface_ltr.5.Top(load.val, v)""",
    }

    solvers = {
        'ls': ('ls.scipy_direct', {}),
        'newton': ('nls.newton', {'eps_a': 1e-6, 'eps_r': 1.e-6,
                                  'check': 0, 'problem': 'nonlinear'}),
    }

    options = {
        'parametric_hook': 'optimization_hook',
        'output_dir' : 'output',
    }

    return locals()
