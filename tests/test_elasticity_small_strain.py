# 10.07.2007, c
# last revision: 25.03.2008
from sfepy import data_dir

filename_meshes = ['/meshes/3d/cube_medium_tetra.mesh',
                   '/meshes/3d/cube_medium_tetra.mesh',
                   '/meshes/3d/cube_medium_hexa.mesh']
filename_meshes = [data_dir + name for name in filename_meshes]

all_your_bases = [1, 2, 1]

filename_mesh = None

field_1 = {
    'name' : '3_displacement',
    'dtype' : 'real',
    'shape' : (3,),
    'region' : 'Omega',
    'approx_order' : None,
}

def get_pars( dim, full = False ):
    import numpy as nm
    sym = (dim + 1) * dim / 2
    lam = 1e1
    mu = 1e0
    o = nm.array( [1.] * dim + [0.] * (sym - dim), dtype = nm.float64 )
    oot = nm.outer( o, o )
    if full:
        return lam * oot + mu * nm.diag( o + 1.0 )
    else:
        return lam, mu

material_1 = {
    'name' : 'solid',
    'values' : {
        'lam' : get_pars( 3 )[0],
        'mu' : get_pars( 3 )[1],
        'Dijkl' : get_pars( 3, True ),
    }
}

material_2 = {
    'name' : 'spring',
    'values' : {
        '.pars' : {'stiffness' : 1e0, 'projection' : None},
    }
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_displacement',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_displacement',
    'dual' : 'u',
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

region_1 = {
    'name' : 'Bottom',
    'select' : 'nodes in (z < -0.499)',
}
region_2 = {
    'name' : 'Top',
    'select' : 'nodes in (z > 0.499)',
}

ebc_1 = {
    'name' : 'Load',
    'region' : 'Top',
    'dofs' : {'u.2' : 0.1},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'order' : 2,
}

equations_iso = {
    'balance_of_forces' :
    """dw_lin_elastic_iso.i1.Omega( solid.lam, solid.mu, v, u )
     = dw_point_lspring.i1.Bottom( spring.pars, v, u )""",
}
equations_general = {
    'balance_of_forces' :
    """dw_lin_elastic.i1.Omega( solid.Dijkl, v, u )
     = dw_point_lspring.i1.Bottom( spring.pars, v, u )""",
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

from sfepy.base.testing import TestCommon

##
# 10.07.2007, c
class Test( TestCommon ):
    tests = ['test_get_solution', 'test_linear_terms']

    ##
    # 10.07.2007, c
    def from_conf( conf, options ):
        return Test( conf = conf, options = options )
    from_conf = staticmethod( from_conf )

    ##
    # c: 25.03.2008, r: 25.03.2008
    def test_linear_terms( self ):
        ok = True
        for sols in self.solutions:
            ok = ok and self.compare_vectors( sols[0], sols[1],
                                             label1 = 'isotropic',
                                             label2 = 'general' )
        return ok
        
    ##
    # c: 10.07.2007, r: 25.03.2008
    def test_get_solution( self ):
        from sfepy.applications import solve_pde
        from sfepy.base.base import IndexedStruct
        import os.path as op

        ok = True
        self.solutions = []
        for ii, approx_order in enumerate(all_your_bases):
            fname = filename_meshes[ii]

            self.conf.filename_mesh = fname
            fields = {'field_1' : {
                          'name' : '3_displacement',
                          'dtype' : 'real',
                          'shape' : (3,),
                          'region' : 'Omega',
                          'approx_order' : approx_order,
                    }
            }
            self.conf.edit('fields', fields)
            self.report( 'mesh: %s, base: %s' % (fname, approx_order) )
            status = IndexedStruct()

            self.report( 'isotropic' )
            self.conf.equations = self.conf.equations_iso
            problem, state1 = solve_pde(self.conf, nls_status=status,
                                        save_results=False)
            converged = status.condition == 0
            ok = ok and converged
            self.report( 'converged: %s' % converged )

            self.report( 'general' )
            self.conf.equations = self.conf.equations_general
            problem, state2 = solve_pde(self.conf, nls_status=status,
                                        save_results=False)
            converged = status.condition == 0
            ok = ok and converged
            self.report( 'converged: %s' % converged )

            self.solutions.append((state1(), state2()))

            name = op.join(self.options.out_dir,
                           '_'.join(('test_elasticity_small_strain',
                                     op.splitext(op.basename(fname))[0],
                                     '%d' % approx_order))
                           + '.vtk')
            problem.save_state(name, state1)

##             trunk = op.join( self.options.out_dir,
##                              op.splitext( op.basename( fname ) )[0] )
##             problem.save_field_meshes( trunk )
##             problem.save_regions( trunk )
            
        return ok
