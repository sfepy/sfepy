# c: 14.04.2008, r: 14.04.2008
import numpy as nm

from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/square_unit_tri.mesh'

def get_pars(ts, coors, mode=None, region=None, ig=None, extra_arg=None):
    if mode == 'special':
        if extra_arg == 'hello!':
            ic = 0
        else:
            ic = 1
        return {('x_%s' % ic) : coors[:,ic]}

def get_p_edge(ts, coors, bc=None):
    if bc.name == 'p_left':
        return nm.sin(nm.pi * coors[:,1])
    else:
        return nm.cos(nm.pi * coors[:,1])

def get_circle(coors, domain=None):
    r = nm.sqrt(coors[:,0]**2.0 + coors[:,1]**2.0)
    return nm.where(r < 0.2)[0]

functions = {
    'get_pars1' : (lambda ts, coors, mode=None, region=None, ig=None:
                   get_pars(ts, coors, mode, region, ig, extra_arg='hello!'),),
    'get_p_edge' : (get_p_edge,),
    'get_circle' : (get_circle,),
}

# Just another way of adding a function, besides 'functions' keyword.
function_1 = {
    'name' : 'get_pars2',
    'function' : lambda ts, coors,mode=None,  region=None, ig=None:
        get_pars(ts, coors, mode, region, ig, extra_arg='hi!'),
}

materials = {
    'mf1' : (None, 'get_pars1'),
    'mf2' : 'get_pars2',
    # Dot denotes a special value, that is not propagated to all QP.
    'mf3' : ({'a' : 10.0, 'b' : 2.0, '.c' : 'ahoj'},),
}

fields = {
    'pressure' : (nm.float64, 1, 'Omega', 2),
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
    'Circle' : ('nodes by get_circle', {}),
}

ebcs = {
    'p_left' : ('Left', {'p.all' : 'get_p_edge'}),
    'p_right' : ('Right', {'p.all' : 'get_p_edge'}),
}

equations = {
    'e1' : """dw_laplace.2.Omega( mf3.a, q, p ) = 0""",
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

        problem = ProblemDefinition.from_conf(conf)
        test = Test(problem = problem, conf = conf, options = options)
        return test
    from_conf = staticmethod( from_conf )


    def test_material_functions(self):
        problem = self.problem
        ts = problem.get_default_ts(step=0)

        problem.materials.time_update(ts,
                                      problem.domain,
                                      problem.equations)

        coors = problem.domain.get_mesh_coors()
        mat1 = problem.materials['mf1']
        assert_(nm.all(coors[:,0] == mat1.get_data(None, None, 'x_0')))

        mat2 = problem.materials['mf2']
        assert_(nm.all(coors[:,1] == mat2.get_data(None, None, 'x_1')))

        mat3 = problem.materials['mf3']
        key = mat3.get_keys(region_name='Omega')[0]
        assert_(nm.all(mat3.get_data(key, 0, 'a') == 10.0))
        assert_(nm.all(mat3.get_data(key, 0, 'b') == 2.0))
        assert_(mat3.get_data(None, None, 'c') == 'ahoj')

        return True
#        mat.time_update(ts, problem)

    def test_ebc_functions(self):
        import os.path as op
        problem = self.problem

        problem.set_equations(self.conf.equations) 

        problem.time_update()
        vec = problem.solve()
        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0] + '_ebc.vtk')
        problem.save_state(name, vec)

        ok = True
        domain = problem.domain

        iv = domain.regions['Left'].get_vertices(0)
        coors = domain.get_mesh_coors()[iv]
        ok = ok and self.compare_vectors(vec[iv], nm.sin(nm.pi * coors[:,1]),
                                         label1='state_left', label2='bc_left')

        iv = domain.regions['Right'].get_vertices(0)
        coors = domain.get_mesh_coors()[iv]
        ok = ok and self.compare_vectors(vec[iv], nm.cos(nm.pi * coors[:,1]),
                                         label1='state_right', label2='bc_right')

        return ok

    def test_region_functions(self):
        import os.path as op
        problem = self.problem

        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0])
        problem.save_regions(name, ['Circle'])

        return True
