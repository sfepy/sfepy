"""Stabilized Navier-Stokes problem with grad-div, SUPG and PSPG stabilization
solved by a custom Oseen solver, see [1].

[1] G. Matthies and G. Lube. On streamline-diffusion methods of inf-sup stable
discretisations of the generalised Oseen problem. Number 2007-02 in Preprint
Series of Institut fuer Numerische und Angewandte Mathematik,
Georg-August-Universitaet Goettingen, 2007.
"""
from sfepy import top_dir

filename_mesh = top_dir + '/meshes/3d/elbow2.mesh'

options = {
    'solution' : 'steady',
    'nls' : 'oseen',
    'ls' : 'ls',
}

regions = {
    'Omega' : ('all', {}),
    'Walls' : ('nodes of surface -n (r.Outlet +n r.Inlet)',
               {'can_cells' : False}),
    'Inlet' : ('nodes by cinc0', {'can_cells' : False}),
    'Outlet' : ('nodes by cinc1', {'can_cells' : False}),
}

fields = {
    'velocity' : ((3,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
    'pressure' : ((1,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
}

variables = {
    'u'   : ('unknown field',   'velocity', 0),
    'v'   : ('test field',      'velocity', 'u'),
    'b'   : ('parameter field', 'velocity', 'u'),
    'p'   : ('unknown field',   'pressure', 1),
    'q'   : ('test field',      'pressure', 'p'),
}

ebcs = {
    'Walls_velocity' : ('Walls', {'u.all' : 0.0}),
    'Inlet_velocity' : ('Inlet', {'u.1' : 1.0, 'u.[0,2]' : 0.0}),
}

materials = {
    'fluid' : ('Omega',
               {'viscosity' : 1.25e-5,
                'density' : 1e0}),
    'stabil' : ('Omega',
                {'.gamma' : None,
                 '.delta' : None,
                 '.tau'   : None,
                 '.tau_red' : 1.0e-0, # <= 1.0; if tau is None: tau = tau_red *
                                     # delta
                 '.tau_mul'   : 1.0,
                 '.delta_mul' : 1.0e-0,
                 '.gamma_mul' : 1.0e0,
                 # 'edge': longest edge, 'volume': volume-based, 'max': max. of
                 # previous
                 '.diameter_mode' : 'max'}), # 'edge', 'volume', 'max'
}

integrals = {
    'i1' : ('v', 'gauss_o2_d3'),
    'i2' : ('v', 'gauss_o3_d3'),
}

##
# Stationary Navier-Stokes equations with grad-div, SUPG and PSPG stabilization.
equations = {
    'balance' :
    """  dw_div_grad.i2.Omega( fluid.viscosity, v, u )
       + dw_lin_convect.i2.Omega( v, b, u )
       - dw_stokes.i1.Omega( v, p )
       + dw_st_grad_div.i2.Omega( stabil.gamma, v, u )
       + dw_st_supg_c.i1.Omega( stabil.delta, v, b, u )
       + dw_st_supg_p.i1.Omega( stabil.delta, v, b, p )
       = 0""",
    'incompressibility' :
    """  dw_stokes.i1.Omega( u, q )
       + dw_st_pspg_c.i1.Omega( stabil.tau, q, b, u )
       + dw_st_pspg_p.i1.Omega( stabil.tau, q, p )
       = 0""",
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 10000
}

solver_1 = {
    'name' : 'oseen',
    'kind' : 'nls.oseen',

    'needs_problem_instance' : True,

    'adimensionalize' : False,
    'check_navier_stokes_rezidual' : False,

    'fluid_mat_name' : 'fluid',
    'stabil_mat_name' : 'stabil',
    'lin_convect_eq_name' : 'balance',
    'div_eq_name' : 'incompressibility',

    'i_max'      : 10,
    'eps_a'      : 1e-8,
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'is_plot'    : False,

    # Uncomment the following to get a convergence log.
    ## 'log'        : {'text' : 'oseen_log.txt',
    ##                 'plot' : 'oseen_log.png'},
}

solver_2 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

##
# Functions.
import os.path as op
import utils

cinc_name = 'cinc_' + op.splitext(op.basename(filename_mesh))[0]
cinc = getattr(utils, cinc_name)

functions = {
    'cinc0' : (lambda coors, domain=None: cinc(coors, 0),),
    'cinc1' : (lambda coors, domain=None: cinc(coors, 1),),
}
