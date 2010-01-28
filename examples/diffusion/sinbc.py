# c: 20.03.2008, r: 20.03.2008
from sfepy import top_dir

filename_mesh = top_dir + '/meshes/various_formats/comsol_tri.txt'

material_1 = {
    'name' : 'm',
    'region' : 'Omega',
    'values' : {'val' : 1.0},
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}
region_3 = {
    'name' : 'Gamma_Bottom',
    'select' : 'nodes in (y < 0.00001)',
}
region_4 = {
    'name' : 'Gamma_Top',
    'select' : 'nodes in (y > 0.59999)',
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P2'}
}

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0,
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 't',
}

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Top',
    'dofs' : {'t.0' : 'ebc_sin'},
}
ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Bottom',
    'dofs' : {'t.0' : 'ebc_sin2'},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

equations = {
    'Temperature' : """dw_laplace.i1.Omega( m.val, s, t )
                       = 0"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
}

fe = {
    'chunk_size' : 1000
}

import numpy as nm

amplitude = 2.0
def ebc_sin(ts, coor, bc):
    x0 = 0.5 * (coor[:,0].min() + coor[:,0].max())
    val = amplitude * nm.sin( (coor[:,0] - x0) * 2. * nm.pi )
    return val

def ebc_sin2(ts, coor, bc):
    x0 = 0.5 * (coor[:,0].min() + coor[:,0].max())
    val = amplitude * nm.sin( (coor[:,0] - x0) * 3. * nm.pi )
    return val

functions = {
    'ebc_sin' : (ebc_sin,),
    'ebc_sin2' : (ebc_sin2,),
}
