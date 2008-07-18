##
# c: 22.02.2008, r: 26.03.2008

#fileName_mesh = 'database/simple.mesh'
#fileName_mesh = 'database/phono/cube_sphere.mesh'
#fileName_mesh = 'database/t.1.node'
file_name_mesh = 'tmp/t.1.node'

options = {
    'save_eig_vectors' : None,
    'squared' : False,
    'n_eigs' : 10,
    'eigen_solver' : 'eigen1',
}

if file_name_mesh.find( 'cube_' ) >= 0:
    # Domain $Y_1$.
    region_1 = {
        'name' : 'Y1',
        'select' : 'elements of group 1',
    }

    # Domain $Y_2$.
    region_2 = {
        'name' : 'Omega',
        'select' : 'elements of group 2',
    }

    # Surface of $Y_2$.
    region_100 = {
        'name' : 'Surface',
        'select' : 'r.Y1 *n r.Omega',
        'can_cells' : False,
    }
else:
    # Whole domain $Y$.
    region_1000 = {
        'name' : 'Omega',
        'select' : 'all',
    }

    # Domain $Y_2$.
    region_2 = {
        'name' : 'Surface',
        'select' : 'nodes of surface',
    }

def get_sphere( x, y, z, mode ):
    import numpy as nm
    r = x**2 + y**2 + z**2
    val = nm.where( (1 <= r) & (r <= 2), 1, 0 )
    return val

region_03 = {
    'name' : 'sphere',
    'select' : 'nodes by get_sphere( x, y, z, 0 )',
}
    
material_1 = {
    'name' : 'm',
    'mode' : 'here',
    'region' : 'Omega',

    'val' : 0.5,
    'one' : 1.0,
}

material_2 = {
    'name' : 'mat_v',
    'mode' : 'function',
    'region' : 'Omega',

    'function' : 'fun_v',
    'extra_args' : {'mode' : 'r^2'},
}

field_0 = {
    'name' : 'field_Psi',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}

variable_1 = {
    'name' : 'Psi',
    'kind' : 'unknown field',
    'field' : 'field_Psi',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : 'field_Psi',
    'dual' : 'Psi',
}
variable_3 = {
    'name' : 'V',
    'kind' : 'parameter field',
    'field' : 'field_Psi',
    'like' : 'Psi',
}
variable_4 = {
    'name' : 'n',
    'kind' : 'parameter field',
    'field' : 'field_Psi',
    'like' : 'Psi',
}

ebc_1 = {
    'name' : 'ZeroSurface',
    'region' : 'Surface',
    'dofs' : {'Psi.0' : 0.0},
}

equations = {
    'lhs' : """  dw_laplace.i1.Omega( m.val, v, Psi )
               + dw_mass_scalar_variable.i1.Omega( mat_v.V, v, Psi )""",
    'rhs' : """dw_mass_scalar.i1.Omega( v, Psi )""",
    'sphere' : """ d_volume_integrate.i1.sphere( n )""",
}

equations_vh = {
    'poisson' : """dw_laplace.i1.Omega( m.one, v, Psi )
                   = dw_mass_scalar.i1.Omega( v, n )"""
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_iterative',

    'method' : 'cg',
    'i_max'      : 1000,
    'eps_a'      : 1e-12,
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
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

solver_2 = {
    'name' : 'eigen1',
    'kind' : 'eig.pysparse',

    'tau' : -10.0,
    'eps_a' : 1e-5,
    'i_max' : 150,
    'method' : 'qmrs',
    'verbosity' : 0,
    'strategy' : 1,
}

fe = {
    'chunk_size' : 100000
}

##
# c: 01.02.2008, r: 12.06.2008
def fun_v( ts, coor, region, ig, mode = None, vhxc = None ):
    import numpy as nm

    if vhxc is None:
        vhxc = 0.0

    out = {}
    C = 0.5
    r = nm.sqrt( coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2 )
    vc = - C * 5.0 / r
    V = vhxc + vc
        
    out['v'] = V
    return out
