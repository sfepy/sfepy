# c: 14.04.2008, r: 14.04.2008
from __future__ import absolute_import
import numpy as nm

from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/square_unit_tri.mesh'

def get_pars(ts, coors, mode=None, extra_arg=None,
             equations=None, term=None, problem=None, **kwargs):
    if mode == 'special':
        if extra_arg == 'hello!':
            ic = 0
        else:
            ic = 1
        coors = problem.get_mesh_coors()
        return {('x_%s' % ic) : coors[:,ic]}

def get_p_edge(ts, coors, bc=None, **kwargs):
    if bc.name == 'p_left':
        return nm.sin(nm.pi * coors[:,1])
    else:
        return nm.cos(nm.pi * coors[:,1])

def get_u_edge(ts, coors, bc=None, **kwargs):
    out = nm.zeros_like(coors)
    out[:, 1] = nm.arange(out.shape[0]) + 1.0
    return out

def get_circle(coors, domain=None):
    r = nm.sqrt(coors[:,0]**2.0 + coors[:,1]**2.0)
    return nm.where(r < 0.2)[0]

functions = {
    'get_pars1' : (lambda ts, coors, mode=None, **kwargs:
                   get_pars(ts, coors, mode, extra_arg='hello!', **kwargs),),
    'get_p_edge' : (get_p_edge,),
    'get_u_edge' : (get_u_edge,),
    'get_circle' : (get_circle,),
}

# Just another way of adding a function, besides 'functions' keyword.
function_1 = {
    'name' : 'get_pars2',
    'function' : lambda ts, coors, mode=None, **kwargs:
        get_pars(ts, coors, mode, extra_arg='hi!', **kwargs),
}

materials = {
    'mf1' : (None, 'get_pars1'),
    'mf2' : 'get_pars2',
    # Dot denotes a special value, that is not propagated to all QP.
    'mf3' : ({'a' : 10.0, 'b' : 2.0, '.c' : 'ahoj'},),
}

fields = {
    'pressure' : (nm.float64, 1, 'Omega', 2),
    'displacement' : (nm.float64, 2, 'Omega', 2),
}

variables = {
    'p'   : ('unknown field', 'pressure', 0),
    'q'   : ('test field',    'pressure', 'p'),
    'u'   : ('unknown field', 'displacement', 1),
    'v'   : ('test field',    'displacement', 'u'),
}

wx = 0.499
regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < -%.3f)' % wx, 'facet'),
    'Right' : ('vertices in (x > %.3f)' % wx, 'facet'),
    'Circle' : 'vertices by get_circle',
}

ebcs = {
    'p_left' : ('Left', {'p.all' : 'get_p_edge'}),
    'p_right' : ('Right', {'p.all' : 'get_p_edge'}),
    'u_right' : ('Right', {'u.all' : 'get_u_edge'}),
}

equations = {
    'e1' : """dw_laplace.2.Omega( mf3.a, q, p ) = 0""",
    'e2' : """dw_div_grad.2.Omega( mf3.b, v, u ) = 0""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',
}

from sfepy.base.base import assert_
from sfepy.base.testing import TestCommon

class Test( TestCommon ):

    def from_conf( conf, options ):
        from sfepy.discrete import Problem

        problem = Problem.from_conf(conf)
        test = Test(problem = problem, conf = conf, options = options)
        return test
    from_conf = staticmethod( from_conf )


    def test_material_functions(self):
        from sfepy.discrete import Material

        problem = self.problem
        conf = problem.conf

        ts = problem.get_default_ts(step=0)

        conf_mat1 = conf.get_item_by_name('materials', 'mf1')
        mat1 = Material.from_conf(conf_mat1, problem.functions)
        mat1.time_update(ts, None, mode='normal', problem=problem)

        coors = problem.domain.get_mesh_coors()
        assert_(nm.all(coors[:,0] == mat1.get_data(None, 'x_0')))

        conf_mat2 = conf.get_item_by_name('materials', 'mf2')
        mat2 = Material.from_conf(conf_mat2, problem.functions)
        mat2.time_update(ts, None, mode='normal', problem=problem)

        assert_(nm.all(coors[:,1] == mat2.get_data(None, 'x_1')))

        materials = problem.get_materials()
        materials.time_update(ts, problem.equations, mode='normal',
                              problem=problem)
        mat3 = materials['mf3']
        key = mat3.get_keys(region_name='Omega')[0]

        assert_(nm.all(mat3.get_data(key, 'a') == 10.0))
        assert_(nm.all(mat3.get_data(key, 'b') == 2.0))
        assert_(mat3.get_data(None, 'c') == 'ahoj')

        return True

    def test_ebc_functions(self):
        import os.path as op
        problem = self.problem

        problem.set_equations(self.conf.equations)

        problem.time_update()
        state = problem.solve()
        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0] + '_ebc.vtk')
        problem.save_state(name, state)

        ok = True
        domain = problem.domain

        vecs = state.get_parts()
        vec = vecs['p']

        iv = domain.regions['Left'].vertices
        coors = domain.get_mesh_coors()[iv]
        _ok = self.compare_vectors(vec[iv], nm.sin(nm.pi * coors[:,1]),
                                   label1='p_state_left',
                                   label2='p_bc_left')
        ok = _ok and ok

        iv = domain.regions['Right'].vertices
        coors = domain.get_mesh_coors()[iv]
        _ok = self.compare_vectors(vec[iv], nm.cos(nm.pi * coors[:,1]),
                                   label1='p_state_right',
                                   label2='p_bc_right')
        ok = _ok and ok

        vec = vecs['u']
        vec.shape = (-1, 2)
        ok = self.compare_vectors(vec[iv, 0], nm.zeros(len(iv)),
                                  label1='u_0_state_right',
                                  label2='u_0_bc_right')
        ok = _ok and ok

        ok = self.compare_vectors(vec[iv, 1], nm.arange(len(iv)) + 1.0,
                                  label1='u_1_state_right',
                                  label2='u_1_bc_right')
        ok = _ok and ok

        return ok

    def test_region_functions(self):
        import os.path as op
        problem = self.problem

        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0])
        problem.save_regions(name, ['Circle'])

        return True
