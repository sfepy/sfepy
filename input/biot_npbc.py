"""Biot problem with the no-penetration BC on the Walls boundary region."""
# c: 10.03.2009
import os
import numpy as nm
from sfepy.mechanics.matcoefs import stiffness_tensor_lame

def define():
    filename = 'database/simple.mesh'
    output_dir = 'output'
    return define_input(filename, output_dir)

def coors_in_cylinder(x, y, z, centre, axis, radius, length, inside=True):
    """
    Select coordinates in a cylinder given by centre, axis and length.
    """
    vec = nm.vstack((x, y, z)) - centre
    
    drv = nm.cross(axis, vec, axisb=0)
    dr = nm.sqrt(nm.sum(drv * drv, 1))
    dl = nm.dot(axis, vec)

    if inside:
        out = nm.where((dl >= 0.0) & (dl <= length) & (dr <= radius), 1, 0)
    else:
        out = nm.where((dl >= 0.0) & (dl <= length) & (dr >= radius), 1, 0)

    return out

def cinc_simple(x, y, z, mode):
    axis = nm.array([1, 0, 0], nm.float64)
    if mode == 0: # In
        centre = nm.array([-0.00001, 0.0, 0.0], nm.float64).reshape((3,1))
        radius = 0.019
        length = 0.00002
    elif mode == 1: # Out
        centre = nm.array([0.09999, 0.0, 0.0], nm.float64).reshape((3,1))
        radius = 0.019
        length = 0.00002
    elif mode == 2: # Rigid
        centre = nm.array([0.05, 0.0, 0.0], nm.float64).reshape((3,1))
        radius = 0.015
        length = 0.03
    else:
        raise ValueError('unknown mode %s!' % mode)

    return coors_in_cylinder(x, y, z, centre, axis, radius, length)

def define_regions(filename):
    if filename.find('simple.mesh'):
        dim, geom = 3, '3_4_P'
        regions = {
            'Omega' : ('all', {}),
            'Walls' : ('nodes of surface -n (r.Outlet +n r.Inlet)',
                       {'can_cells' : True}),
            'Inlet' : ('nodes by cinc_simple( x, y, z, 0 )',
                       {'can_cells' : False}),
            'Outlet' : ('nodes by cinc_simple( x, y, z, 1 )',
                        {'can_cells' : False}),
            'Rigid' : ('nodes by cinc_simple( x, y, z, 2 )', {}),
        }

    else:
        raise ValueError('unknown mesh %s!' % filename)

    return regions, dim, geom

def get_pars(ts, coor, region, ig, output_dir='.'):
    n_nod, dim = coor.shape
    sym = (dim + 1) * dim / 2

    out = {}
    out['D'] = stiffness_tensor_lame(dim, lam=1.7, mu=0.3)

    aa = nm.zeros((sym,), dtype=nm.float64)
    aa[:dim] = 0.132
    aa[dim:sym] = 0.092
    out['alpha'] = aa
    
    perm = nm.eye(dim, dtype=nm.float64)
    out['K'] = perm
    
    return out

def post_process(out, pb, state, extend=False):
    from sfepy.base.base import Struct
    from sfepy.fem.evaluate import eval_term_op

    dvel = eval_term_op(state,
                        'de_diffusion_velocity.i1.Omega( m.K, p )', pb)
    out['dvel'] = Struct(name='output_data',
                         mode='cell', data=dvel,
                         dof_types=None)

    stress = eval_term_op(state, 'de_cauchy_stress.i1.Omega( m.D, u )', pb)
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress,
                                  dof_types=None)
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

    regions, dim, geom = define_regions(filename_mesh)

    field_1 = {
        'name' : 'displacement',
        'dim' : (dim,1),
        'domain' : 'Omega',
        'bases' : {'Omega' : '%s1' % geom}
    }
    field_2 = {
        'name' : 'pressure',
        'dim' : (1,1),
        'domain' : 'Omega',
        'bases' : {'Omega' : '%s1' % geom}
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
        'mode' : 'function',
        'region' : 'Omega',
        'function' : 'get_pars',
        'extra_args' : {'output_dir' : output_dir},
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

    fe = {
        'chunk_size' : 100000,
        'cache_override' : True,
    }

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.umfpack', # Direct solver.
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',
    }

    return locals()
