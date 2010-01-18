# -*- coding: utf-8 -*-
import numpy as nm

filename_mesh = '../database/simple.vtk'

# Time-stepping parameters.
t0 = 0.0
t1 = 1.0
n_step = 21

from sfepy.solvers.ts import TimeStepper
ts = TimeStepper(t0, t1, None, n_step)

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
    'post_process_hook' : 'post_process',
}


fields = {
    'displacement': ((3,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
    'pressure'    : ((1,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
}

materials = {
    # Perfused solid.
    'ps' : ('Omega', {
        'mu' : 20e0, # shear modulus of neoHookean term
        'k'  : ts.dt * nm.eye(3, dtype=nm.float64), # reference permeability
        'N_f' : 1.0, # reference porosity
    }),
    # Surface pressure traction.
    'traction' : ('Right', None, 'get_traction'),
}

variables = {
    'u' : ('unknown field', 'displacement', 0, 'previous'),
    'v' : ('test field', 'displacement', 'u'),
    'p' : ('unknown field', 'pressure', 1),
    'q' : ('test field', 'pressure', 'p'),
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001)', {}),
    'Right' : ('nodes in (x > 0.099)', {}),
}

##
# Dirichlet BC.
ebcs = {
    'l' : ('Left', {'u.all' : 0.0, 'p.0' : 'get_pressure'}),
}

##
# Balance of forces.
integrals = {
    'i1' : ('v', 'gauss_o1_d3'),
    'i2' : ('s3', 'gauss_o2_d2'),
}

equations = {
    'force_balance'
        : """dw_tl_he_neohook.i1.Omega( ps.mu, v, u )
           + dw_tl_bulk_pressure.i1.Omega( v, u, p )
           + dw_tl_surface_traction.i2.Right( traction.pressure, v, u )
           = 0""",
    'mass_balance'
        : """dw_tl_volume.i1.Omega( q, u )
           + dw_tl_diffusion.i1.Omega( ps.k, ps.N_f, q, p, u[-1])
           = dw_tl_volume.i1.Omega( q, u[-1] )"""
}

def post_process(out, problem, state, extend=False):
    from sfepy.base.base import Struct, debug

    val = problem.evaluate('dw_tl_he_neohook.i1.Omega( ps.mu, v, u )',
                           state, call_mode='de_strain')
    out['green_strain'] = Struct(name = 'output_data',
                                 mode = 'cell', data = val,
                                 dof_types = None)

    val = problem.evaluate('dw_tl_he_neohook.i1.Omega( ps.mu, v, u )',
                           state, call_mode='de_stress')
    out['neohook_stress'] = Struct(name = 'output_data',
                                   mode = 'cell', data = val,
                                   dof_types = None)

    val = problem.evaluate('dw_tl_bulk_pressure.i1.Omega( v, u, p )',
                           state, call_mode='de_stress')
    out['bulk_pressure'] = Struct(name = 'output_data',
                                  mode = 'cell', data = val,
                                  dof_types = None)

    val = problem.evaluate('dw_tl_diffusion.i1.Omega( ps.k, ps.N_f, q, p, u[-1] )',
                           state, call_mode='de_diffusion_velocity')
    out['diffusion_velocity'] = Struct(name = 'output_data',
                                       mode = 'cell', data = val,
                                       dof_types = None)

    return out

##
# Solvers etc.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 7,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'    : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp': 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'      : 0,
    'delta'      : 1e-6,
    'is_plot'    : False,
    'problem'    : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

solver_2 = {
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'    : t0,
    't1'    : t1,
    'dt'    : None,
    'n_step' : n_step, # has precedence over dt!
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 100000,
    'cache_override' : False,
}

##
# Functions.
def get_traction(ts, coors, mode=None):
    """
    Pressure traction.
    
    Parameters
    ----------
    ts : TimeStepper
        Time stepping info.
    coors : array_like
        The physical domain coordinates where the parameters shound be defined.
    mode : 'qp' or 'special'
        Call mode.
    """
    if mode != 'qp': return

    tt = ts.nt * 2.0 * nm.pi

    dim = coors.shape[1]
    val = 1e-1 * nm.sin(tt) * nm.eye(dim, dtype=nm.float64)

    shape = (coors.shape[0], 1, 1)
    out = {
        'pressure' : nm.tile(val, shape),
    }

    return out

def get_pressure(ts, coor, bc):
    """Internal pressure Dirichlet boundary condition."""
    tt = ts.nt * 2.0 * nm.pi

    val = nm.zeros((coor.shape[0],), dtype=nm.float64)

    val[:] = 1e-2 * nm.sin(tt)

    return val

functions = {
    'get_traction' : (lambda ts, coors, mode=None, region=None, ig=None:
                      get_traction(ts, coors, mode=mode),),
    'get_pressure' : (get_pressure,),
}

