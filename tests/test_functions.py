# c: 14.04.2008, r: 14.04.2008
import numpy as nm
filename_mesh = '../database/tests/plane.mesh'

def get_pars(ts, coor, region=None, ig=None, extra_arg=None):
    print extra_arg
    if extra_arg == 'hello!':
        ic = 0
    else:
        ic = 1
    return {('x_%s' % ic) : coor[:,ic]}

def get_p_edge(ts, coor, bc=None):
    if bc.name == 'p_left':
        return nm.sin(nm.pi * coor[:,1])
    else:
        return nm.cos(nm.pi * coor[:,1])

functions = {
    'get_pars1' : (lambda ts, coor, region=None, ig=None:
                   get_pars(ts, coor, region, ig, extra_arg='hello!'),),
    'get_p_edge' : (get_p_edge,),
}

function_1 = {
    'name' : 'get_pars2',
    'function' : lambda ts, coor, region=None, ig=None:
        get_pars(ts, coor, region, ig, extra_arg='hi!'),
}

materials = {
    'mf1' : ('Omega', None, 'get_pars1'),
    'mf2' : ('Omega', None, 'get_pars2'),
    'mf3' : ('Omega', {'a' : 1.0, 'b' : 2.0}),
}

fields = {
    'pressure' : ((1,1), 'real', 'Omega', {'Omega' : '2_3_P2'}),
}

variables = {
    'p'   : ('unknown field', 'pressure', 0),
    'q'   : ('test field',    'pressure', 'p'),
}

wx = 0.499
regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < -%.3f)' % wx, {}),
    'Right' : ('nodes in (x > %.3f)' % wx, {}),
}

integrals = {
    'i1' : ('v', 'gauss_o2_d2'),
}

ebcs = {
    'p_left' : ('Left', {'p.all' : 'get_p_edge'}),
    'p_right' : ('Right', {'p.all' : 'get_p_edge'}),
}

equations = {
    'e1' : """dw_laplace.i1.Omega( mf3.a, q, p ) = 0""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',
}

fe = {
    'chunk_size' : 1000
}

from sfepy.base.testing import TestCommon, assert_
from sfepy.base.base import pause, debug

class Test( TestCommon ):

    def from_conf( conf, options ):
        from sfepy.fem import ProblemDefinition

        problem = ProblemDefinition.from_conf(conf, init_variables=False)
        test = Test(problem = problem, conf = conf, options = options)
        return test
    from_conf = staticmethod( from_conf )


    def test_material_functions(self):
        problem = self.problem
        ts = problem.get_default_ts(step=0)
        problem.materials.time_update(ts, problem)

        coor = problem.domain.get_mesh_coors()
        
        mat1 = problem.materials['mf1']
        assert_(nm.all(coor[:,0] == mat1.datas[0]['x_0']))

        mat2 = problem.materials['mf2']
        assert_(nm.all(coor[:,1] == mat2.datas[0]['x_1']))

        mat3 = problem.materials['mf3']
        assert_(mat3.datas[0]['a'] == 1.0)
        assert_(mat3.datas[0]['b'] == 2.0)

        return True
#        mat.time_update(ts, problem)

    def test_ebc_functions(self):
        import os.path as op
        problem = self.problem

        problem.set_variables(self.conf.variables) 
        problem.set_equations(self.conf.equations) 

        problem.time_update()
        vec = problem.solve()
        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0] + '_ebc.vtk')
        problem.save_state(name, vec)

        ok = True
        domain = problem.domain

        iv = domain.regions['Left'].get_vertices(0)
        coor = domain.get_mesh_coors()[iv]
        ok = ok and self.compare_vectors(vec[iv], nm.sin(nm.pi * coor[:,1]),
                                         label1='state_left', label2='bc_left')

        iv = domain.regions['Right'].get_vertices(0)
        coor = domain.get_mesh_coors()[iv]
        ok = ok and self.compare_vectors(vec[iv], nm.cos(nm.pi * coor[:,1]),
                                         label1='state_right', label2='bc_right')

        return ok
