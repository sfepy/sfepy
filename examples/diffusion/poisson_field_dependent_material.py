r"""
Stationary Laplace equation. 
Find :math:`T(t)` for :math:`t \in [0, t_{\rm final}]` such that:

.. math::
   \int_{\Omega} c(T) \nabla s \cdot \nabla T
    = 0
    \;, \quad \forall s \;.
	
where :math:`c(T)` is the :math:`T` dependent diffusion coefficient.
Each iterations calculates :math:`T` and adjusts :math:`c(T)`.
"""
from sfepy import data_dir
import numpy as np
filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'
print data_dir
t0 = 0.0
t1 = 0.1
n_step = 11

def get_conductivity(ts, coors, problem, equations=None,mode=None, **kwargs):
	"""
	Calculates the conductivity as 2+10*T and returns it. 
	This relation results in larger T gradients where T is small.
	"""
	if mode == 'qp':
		try:
			T_values=equations.variables._objs.find('T').evaluate_at(coors) # T-field values interpolated at the cell coordinates
			val=(2+10*(T_values+2))
		except:
			val=(22*np.ones(len(coors[:,0]))) # in case no T values can be found use T=0 everywhere to calculate conductivity.
		val.shape=(coors.shape[0], 1, 1)
		return {'val':val}

material_2 = {
    'name' : 'coef',
    'function':'get_conductivity',
    'kind' : 'time-dependent', # 'stationary' or 'time-dependent'
}

field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

variable_1 = {
    'name' : 'T',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0,
    'history' : 1,
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 'T',
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_Left' : ('nodes in (x < 0.00001)', {}),
    'Gamma_Right' : ('nodes in (x > 0.099999)', {}),
}

ebcs = {
    'T1': ('Gamma_Left', {'T.0' : 2.0}),
    'T2': ('Gamma_Right', {'T.0' : -2.0}),
}


functions = {
    'get_conductivity' : (get_conductivity,),
}
    
ics = {
    'ic' : ('Omega', {'T.0' : 0}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'order' : 1,
}

equations = {
    'Temperature' :
    """ dw_laplace.i1.Omega( coef.val, s, T ) = 0"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',

    'presolve' : True,
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
    'problem'   : 'linear', # 'nonlinear' or 'linear' (ignore i_max)
}

solver_2 = {
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'    : t0,
    't1'    : t1,
    'dt'    : None,
    'n_step' : n_step, # has precedence over dt!
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
}
