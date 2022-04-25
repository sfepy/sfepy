import numpy as nm
import pytest

from sfepy import data_dir
from sfepy.base.base import assert_
import sfepy.base.testing as tst

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

    elif mode == 'qp':
        return {'a' : nm.tile(-2 + 1j, (coors.shape[0], 1, 1))}

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
    'get_pars' : (get_pars,),
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
    # Complex values.
    'mf4' : 'get_pars',
    'mf5' : ({'a' : -2 - 1j},),
    'mf6' : ({'a' : {'Circle' : 1 + 1j, 'Rest' : 3j}},),
}

fields = {
    'pressure' : (nm.float64, 1, 'Omega', 2),
    'displacement' : (nm.float64, 2, 'Omega', 2),
    'complex' : (nm.complex128, 1, 'Omega', 2),
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
    'Rest' : 'r.Omega -c r.Circle',
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

variables2 = {
    'r'   : ('unknown field', 'complex', 2),
    's'   : ('test field',    'complex', 'r'),
}

equations2 = {
    'e3' : """dw_laplace.2.Omega(mf4.a, s, r) = 0""",
    'e4' : """dw_laplace.2.Omega(mf5.a, s, r) = 0""",
    'e5' : """dw_laplace.2.Circle(mf6.a, s, r)
            + dw_laplace.2.Rest(mf6.a, s, r) = 0""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',
}

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf)
    return problem

def test_material_functions(problem):
    from sfepy.discrete import Material
    from sfepy.base.conf import transform_variables

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

    pb = problem.copy()
    pb.set_variables(transform_variables(conf.variables2))
    pb.set_equations(conf.equations2)
    materials = pb.get_materials()
    materials.time_update(ts, pb.equations, mode='normal', problem=pb)

    mat4 = materials['mf4']
    key = mat4.get_keys(region_name='Omega')[0]
    assert_(nm.all(mat4.get_data(key, 'a') == -2 + 1j))

    mat5 = materials['mf5']
    key = mat5.get_keys(region_name='Omega')[0]
    assert_(nm.all(mat5.get_data(key, 'a') == -2 - 1j))

    mat6 = materials['mf6']
    key = mat6.get_keys(region_name='Circle')[0]
    assert_(nm.all(mat6.get_data(key, 'a') == 1 + 1j))
    key = mat6.get_keys(region_name='Rest')[0]
    assert_(nm.all(mat6.get_data(key, 'a') == 3j))

def test_ebc_functions(problem, output_dir):
    import os.path as op

    problem.set_equations(problem.conf.equations)

    problem.time_update()
    state = problem.solve(save_results=False)
    name = op.join(output_dir, 'test_ebc_functions.vtk')
    problem.save_state(name, state)

    ok = True
    domain = problem.domain

    vecs = state.get_state_parts()
    vec = vecs['p']

    iv = domain.regions['Left'].vertices
    coors = domain.get_mesh_coors()[iv]
    _ok = tst.compare_vectors(vec[iv], nm.sin(nm.pi * coors[:,1]),
                              label1='p_state_left',
                              label2='p_bc_left')
    ok = _ok and ok

    iv = domain.regions['Right'].vertices
    coors = domain.get_mesh_coors()[iv]
    _ok = tst.compare_vectors(vec[iv], nm.cos(nm.pi * coors[:,1]),
                              label1='p_state_right',
                              label2='p_bc_right')
    ok = _ok and ok

    vec = vecs['u']
    vec.shape = (-1, 2)
    ok = tst.compare_vectors(vec[iv, 0], nm.zeros(len(iv)),
                             label1='u_0_state_right',
                             label2='u_0_bc_right')
    ok = _ok and ok

    ok = tst.compare_vectors(vec[iv, 1], nm.arange(len(iv)) + 1.0,
                             label1='u_1_state_right',
                             label2='u_1_bc_right')
    ok = _ok and ok

    assert ok

def test_region_functions(problem, output_dir):
    import os.path as op

    name = op.join(output_dir, 'test_region_functions.vtk')
    problem.save_regions_as_groups(name, ['Circle'])
