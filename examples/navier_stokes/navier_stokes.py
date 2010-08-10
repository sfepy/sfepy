# 26.02.2007, c
# last revision: 25.02.2008
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/elbow2.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'post_process_hook' : 'verify_incompressibility',
}

field_1 = {
    'name' : '3_velocity',
    'dtype' : 'real',
    'shape' : (3,),
    'region' : 'Omega',
    'approx_order' : '1B',
}

field_2 = {
    'name' : 'pressure',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

# Can use logical operations '&' (and), '|' (or).
region_1000 = {
    'name' : 'Omega',
    'select' : 'elements of group 6',
}

region_0 = {
    'name' : 'Walls',
    'select' : 'nodes of surface -n (r.Outlet +n r.Inlet)',
    'can_cells' : False,
}
region_1 = {
    'name' : 'Inlet',
    'select' : 'nodes by cinc0', # In
    'can_cells' : False,
}
region_2 = {
    'name' : 'Outlet',
    'select' : 'nodes by cinc1', # Out
    'can_cells' : False,
}

ebc_1 = {
    'name' : 'Walls',
    'region' : 'Walls',
    'dofs' : {'u.all' : 0.0},
}
ebc_2 = {
    'name' : 'Inlet',
    'region' : 'Inlet',
    'dofs' : {'u.1' : 1.0, 'u.[0,2]' : 0.0},
}

material_1 = {
    'name' : 'fluid',
    'values' : {
        'viscosity' : 1.25e-3,
        'density' : 1e0,
    },
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_velocity',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_velocity',
    'dual' : 'u',
}
variable_3 = {
    'name' : 'p',
    'kind' : 'unknown field',
    'field' : 'pressure',
    'order' : 1,
}
variable_4 = {
    'name' : 'q',
    'kind' : 'test field',
    'field' : 'pressure',
    'dual' : 'p',
}
variable_5 = {
    'name' : 'pp',
    'kind' : 'parameter field',
    'field' : 'pressure',
    'like' : 'p',
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}
integral_2 = {
    'name' : 'i2',
    'kind' : 'v',
    'quadrature' : 'gauss_o3_d3',
}

##
# Stationary Navier-Stokes equations.
equations = {
    'balance' :
    """+ dw_div_grad.i2.Omega( fluid.viscosity, v, u )
       + dw_convect.i2.Omega( v, u )
       - dw_stokes.i1.Omega( v, p ) = 0""",
    'incompressibility' :
    """dw_stokes.i1.Omega( u, q ) = 0""",
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 1000
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 5,
    'eps_a'      : 1e-8,
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 0.99999,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

def verify_incompressibility( out, problem, state, extend = False ):
    """This hook is normally used for post-processing (additional results can
    be inserted into `out` dictionary), but here we just verify the weak
    incompressibility condition."""
    from sfepy.base.base import Struct, debug, nm

    vv = problem.get_variables()
    one = nm.ones( (vv['p'].field.n_nod,), dtype = nm.float64 )
    vv['p'].data_from_any( one )
    zero = problem.evaluate('dw_stokes.i1.Omega( u, p )', p=one, u=vv['u'](),
                            call_mode='d_eval')
    print 'div( u ) = %.3e' % zero

    return out

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
