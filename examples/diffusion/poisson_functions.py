"""
Example demonstrating use of functions for defining material parameters,
regions, parameter variables or boundary conditions.
"""
import numpy as nm
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

materials = {
    'm' : ({'c' : 1.0},),
    'load' : 'get_pars',
}

regions = {
    'Omega' : ('all', {}),
    'Omega_L' : ('nodes by get_middle_ball', {}),
    'Gamma_Left' : ('nodes in (x < 0.00001)', {}),
    'Gamma_Right' : ('nodes in (x > 0.099999)', {}),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    'u' : ('unknown field', 'temperature', 0),
    'v' : ('test field',    'temperature', 'u'),
    'p' : ('parameter field', 'temperature',
           {'setter' : 'get_load_variable'}),
}

ebcs = {
    'u1' : ('Gamma_Left', {'u.0' : 'get_ebc'}),
    'u2' : ('Gamma_Right', {'u.0' : -2.0}),
}

integrals = {
    'i1' : ('v', 'gauss_o1_d3'),
}

equations = {
    'Laplace equation' :
    """dw_laplace.i1.Omega( m.c, v, u )
     = - dw_mass_scalar.i1.Omega_L( load.val, v, p )"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
                 }),
}

def get_pars(ts, coors, mode=None, **kwargs):
    """
    We can define the coefficient `load.val` as a function of space.

    For scalar parameters, the shape has to be set to `(coors.shape[0], 1, 1)`.
    """
    if mode == 'qp':
        x = coors[:,0]

        val = 55.0 * (x - 0.05)

        val.shape = (coors.shape[0], 1, 1)
        return {'val' : val}

def get_middle_ball(coors, domain=None):
    x, y, z = coors[:,0], coors[:,1], coors[:,2]

    r1 = nm.sqrt((x - 0.025)**2.0 + y**2.0 + z**2)
    r2 = nm.sqrt((x - 0.075)**2.0 + y**2.0 + z**2)
    flag = nm.where((r1 < 2.3e-2) | (r2 < 2.3e-2))[0]

    return flag

def get_load_variable(ts, coors, region=None):
    y = coors[:,1]

    val = 5e5 * y

    return val

def get_ebc(coor, amplitude):
    z = coor[:,2]
    val = amplitude * nm.sin(z * 2.0 * nm.pi)
    return val

functions = {
    'get_pars' : (get_pars,),
    'get_load_variable' : (get_load_variable,),
    'get_middle_ball' : (get_middle_ball,),
    'get_ebc' : (lambda ts, coor, bc, problem, **kwargs: get_ebc(coor, 5.0),),
}
