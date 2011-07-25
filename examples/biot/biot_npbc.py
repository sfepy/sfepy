"""Biot problem with the no-penetration BC on the Walls boundary region."""
# c: 10.03.2009
import os
import numpy as nm

from sfepy.linalg import get_coors_in_tube
from sfepy.mechanics.matcoefs import stiffness_tensor_lame

def define():
    from sfepy import data_dir

    filename = data_dir + '/meshes/3d/cylinder.mesh'
    output_dir = 'output'
    return define_input(filename, output_dir)

def cinc_simple(coors, mode):
    axis = nm.array([1, 0, 0], nm.float64)
    if mode == 0: # In
        centre = nm.array([-0.00001, 0.0, 0.0], nm.float64)
        radius = 0.019
        length = 0.00002
    elif mode == 1: # Out
        centre = nm.array([0.09999, 0.0, 0.0], nm.float64)
        radius = 0.019
        length = 0.00002
    elif mode == 2: # Rigid
        centre = nm.array([0.05, 0.0, 0.0], nm.float64)
        radius = 0.015
        length = 0.03
    else:
        raise ValueError('unknown mode %s!' % mode)

    return get_coors_in_tube(coors,
                             centre, axis, -1, radius, length)

def define_regions(filename):
    if filename.find('simple.mesh'):
        dim = 3
        regions = {
            'Omega' : ('all', {}),
            'Walls' : ('nodes of surface -n (r.Outlet +n r.Inlet)',
                       {'can_cells' : True}),
            'Inlet' : ('nodes by cinc_simple0', {'can_cells' : False}),
            'Outlet' : ('nodes by cinc_simple1', {'can_cells' : False}),
            'Rigid' : ('nodes by cinc_simple2', {}),
        }

    else:
        raise ValueError('unknown mesh %s!' % filename)

    return regions, dim

def get_pars(ts, coor, mode, output_dir='.', **kwargs):
    if mode == 'qp':
        n_nod, dim = coor.shape
        sym = (dim + 1) * dim / 2

        out = {}
        out['D'] = nm.tile(stiffness_tensor_lame(dim, lam=1.7, mu=0.3),
                           (coor.shape[0], 1, 1))

        aa = nm.zeros((sym, 1), dtype=nm.float64)
        aa[:dim] = 0.132
        aa[dim:sym] = 0.092
        out['alpha'] = nm.tile(aa, (coor.shape[0], 1, 1))

        perm = nm.eye(dim, dtype=nm.float64)
        out['K'] = nm.tile(perm, (coor.shape[0], 1, 1))

        return out

def post_process(out, pb, state, extend=False):
    from sfepy.base.base import Struct

    dvel = pb.evaluate('de_diffusion_velocity.i1.Omega( m.K, p )',
                       mode='el_avg')
    out['dvel'] = Struct(name='output_data',
                         mode='cell', data=dvel, dofs=None)

    stress = pb.evaluate('de_cauchy_stress.i1.Omega( m.D, u )',
                         mode='el_avg')
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress, dofs=None)
    return out

def define_input(filename, output_dir):
    
    filename_mesh = filename
    options = {
        'output_dir' : output_dir,
        'output_format' : 'vtk',
        'post_process_hook' : 'post_process',

        'ls' : 'ls',
        'nls' : 'newton',
    }

    functions = {
        'cinc_simple0' : (lambda coors, domain:
                          cinc_simple(coors, 0),),
        'cinc_simple1' : (lambda coors, domain:
                          cinc_simple(coors, 1),),
        'cinc_simple2' : (lambda coors, domain:
                          cinc_simple(coors, 2),),
        'get_pars' : (lambda ts, coors, mode=None, **kwargs:
                      get_pars(ts, coors, mode,
                               output_dir=output_dir, **kwargs),),
    }
    regions, dim = define_regions(filename_mesh)

    field_1 = {
        'name' : 'displacement',
        'dtype' : nm.float64,
        'shape' : dim,
        'region' : 'Omega',
        'approx_order' : 1,
    }
    field_2 = {
        'name' : 'pressure',
        'dtype' : nm.float64,
        'shape' : 1,
        'region' : 'Omega',
        'approx_order' : 1,
    }

    variables = {
        'u'       : ('unknown field',   'displacement', 0),
        'v'       : ('test field',      'displacement', 'u'),
        'p'       : ('unknown field',   'pressure', 1),
        'q'       : ('test field',      'pressure', 'p'),
    }

    ebcs = {
        'inlet' : ('Inlet', {'p.0' : 1.0, 'u.all' : 0.0}),
        'outlet' : ('Outlet', {'p.0' : -1.0}),
    }

    lcbcs = {
        'rigid' : ('Outlet', {'u.all' : 'rigid'}),
        'no_penetration' : ('Walls', {'u.all' : 'no_penetration'}),
    }

    material_1 = {
        'name' : 'm',
        'function' : 'get_pars',
    }

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o2_d%d' % dim,
    }

    equations = {
        'eq_1' : 
        """dw_lin_elastic.i1.Omega( m.D, v, u )
         - dw_biot.i1.Omega( m.alpha, v, p ) 
         = 0""",
        'eq_2' :
        """dw_biot.i1.Omega( m.alpha, u, q )
         + dw_diffusion.i1.Omega( m.K, q, p )
         = 0""",
    }

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct', # Direct solver.
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',
    }

    return locals()
