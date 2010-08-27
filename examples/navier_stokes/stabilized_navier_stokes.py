"""Stabilized Navier-Stokes problem with grad-div, SUPG and PSPG stabilization
solved by a custom Oseen solver, see [1].

[1] G. Matthies and G. Lube. On streamline-diffusion methods of inf-sup stable
discretisations of the generalised Oseen problem. Number 2007-02 in Preprint
Series of Institut fuer Numerische und Angewandte Mathematik,
Georg-August-Universitaet Goettingen, 2007.
"""
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/elbow2.mesh'

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
    'velocity' : ('real', 3, 'Omega', 1),
    'pressure' : ('real', 1, 'Omega', 1),
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
    'fluid' : ({'viscosity' : 1.25e-5,
                'density' : 1e0},),
    'stabil' : ({'.gamma' : None,
                 '.delta' : None,
                 '.tau'   : None,
                 '.tau_red' : 1.0e-0, # <= 1.0; if tau is None: tau = tau_red *
                                     # delta
                 '.tau_mul'   : 1.0,
                 '.delta_mul' : 1.0e-0,
                 '.gamma_mul' : 1.0e0,
                 # 'edge': longest edge, 'volume': volume-based, 'max': max. of
                 # previous
                 '.diameter_mode' : 'max'},), # 'edge', 'volume', 'max'
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
    'stabilization_hook' : 'create_stabil_mat',

    'adimensionalize' : False,
    'check_navier_stokes_rezidual' : False,

    'i_max'      : 10,
    'eps_a'      : 1e-8,
    'eps_r'      : 1.0,
    'macheps'    : 1e-16,
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
import numpy as nm

import utils

cinc_name = 'cinc_' + op.splitext(op.basename(filename_mesh))[0]
cinc = getattr(utils, cinc_name)

def create_stabil_mat(problem):
    """Using the stabilization material stub make it the true material."""
    from sfepy.base.base import dict_to_struct, debug
    from sfepy.fem.functions import Function

    # Identity map...
    ns = {'p' : 'p', 'q' : 'q',
          'u' : 'u', 'b' : 'b', 'v' : 'v',
          'fluid' : 'fluid', 'omega' : 'omega', 'i1' : 'i1', 'i2' : 'i2'}

    variables = problem.get_variables()

    # Indices to the state vector.
    ii = {}
    ii['u'] = variables.get_indx(ns['u'])
    ii['us'] = variables.get_indx(ns['u'], stripped=True)
    ii['ps'] = variables.get_indx(ns['p'], stripped=True)

    stabil_mat = problem.materials['stabil']
    stabil = dict_to_struct(stabil_mat.datas['special'], flag=(1,))

    # The viscosity.
    fluid_mat = problem.materials[ns['fluid']]
    viscosity = fluid_mat.function()['viscosity']

    # The Friedrich's constant.
    c_friedrichs = problem.domain.get_diameter()
    sigma = 1e-12 # 1 / dt.

    # Element diameter modes.
    diameter_modes = {'edge' : 0, 'volume' : 1, 'max' : 2}

    def mat_fun(ts, coor, mode=None, region=None, ig=None,
                b_norm=1.0):
        if mode != 'qp': return

        print '|b|_max (mat_fun):', b_norm
        gamma = viscosity + b_norm * c_friedrichs

        data = {}
        if stabil.gamma is None:
            data['gamma'] = stabil.gamma_mul * gamma
        else:
            data['gamma'] = nm.asarray( stabil.gamma_mul * stabil.gamma,
                                        dtype = nm.float64 )
        data['gamma'] = nm.tile(data['gamma'], (coor.shape[0], 1, 1))

        if stabil.delta is None:
            term = problem.equations['balance'].terms['dw_lin_convect']
            for ig in term.iter_groups():
                # This sets term.ig - for 1 group only!!!
                break
            var = variables[ns['u']]
            ap, vg = term.get_approximation(var)
            delta = 1.0
            mode = diameter_modes[stabil.diameter_mode]
            cells = region.get_cells( ig )
            diameters2 = problem.domain.get_element_diameters( ig, cells, vg,
                                                             mode )
            val1 = min( 1.0, 1.0 / sigma )
            val2 = sigma * c_friedrichs**2
            val3 = (b_norm**2) * min( (c_friedrichs**2) / viscosity, 1.0 / sigma )
#            print val1, gamma, val2, val3
            delta = stabil.delta_mul * val1 * diameters2 / (gamma + val2 + val3)

            n_qp = coor.shape[0] / diameters2.shape[0]
            data['diameters2'] = nm.repeat(diameters2, n_qp)
            data['diameters2'].shape = data['diameters2'].shape + (1, 1)

            data['delta'] = nm.repeat(delta, n_qp)
            data['delta'].shape = data['delta'].shape + (1, 1)
        else:
            val = stabil.delta_mul * stabil.delta
            data['delta'] = nm.tile(data['delta'], (coor.shape[0], 1, 1))
        
        if stabil.tau is None:
            data['tau'] = stabil.tau_red * data['delta']
        else:
            data['tau'] = nm.asarray( stabil.tau_mul * stabil.tau,
                                      dtype = nm.float64 )
            data['tau'] = nm.tile(data['tau'], (coor.shape[0], 1, 1))

        return data

    stabil_mat.set_function(Function('stabil', mat_fun))

    return stabil_mat, ns, ii

functions = {
    'cinc0' : (lambda coors, domain=None: cinc(coors, 0),),
    'cinc1' : (lambda coors, domain=None: cinc(coors, 1),),
    'create_stabil_mat' : (create_stabil_mat,),
}
