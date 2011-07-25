# c: 15.02.2008, r: 28.03.2008
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/various_formats/octahedron.node'

material_2 = {
    'name' : 'coef',
    'values' : {'K' : [[1.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 1.0]],
                'val' : 1.0},
}

field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
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

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

def get_line(x, y, mode):
    import numpy as nm

    if mode == 0:
        val = nm.where((x + y) >= 3.5)[0]
    elif mode == 1:
        val = nm.where((x + y) <= -3.5)[0]
    print mode, val
    return val

functions = {
    'get_line0' : (lambda coors, domain=None:
                   get_line(coors[:,0], coors[:,1], 0),),
    'get_line1' : (lambda coors, domain=None:
                   get_line(coors[:,0], coors[:,1], 1),),
}

region_03 = {
    'name' : 'Gamma_Left',
    'select' : 'nodes by get_line0',
}
region_4 = {
    'name' : 'Gamma_Right',
    'select' : 'nodes by get_line1',
}

ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Left',
    'dofs' : {'t.0' : 2.0},
}

ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Right',
    'dofs' : {'t.0' : -2.0},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'Temperature' : """dw_diffusion.i1.Omega( coef.K, s, t ) = 0"""
#    'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
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
