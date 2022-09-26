import os.path as op
import numpy as nm
import pytest

import sfepy.base.testing as tst

def fix_u_fun(ts, coors, bc=None, problem=None, extra_arg=None):
    return nm.zeros_like(coors)

@pytest.fixture(scope='module')
def data():
    import sfepy
    from sfepy.base.base import Struct
    from sfepy.discrete.fem import Mesh, FEDomain, Field
    mesh = Mesh.from_file('meshes/2d/rectangle_tri.mesh',
                          prefix_dir=sfepy.data_dir)
    domain = FEDomain('domain', mesh)
    dim = domain.shape.dim

    min_x, max_x = domain.get_mesh_bounding_box()[:,0]
    eps = 1e-8 * (max_x - min_x)

    omega = domain.create_region('Omega', 'all')
    gamma1 = domain.create_region('Gamma1',
                                  'vertices in x < %.10f' % (min_x + eps),
                                  'facet')
    gamma2 = domain.create_region('Gamma2',
                                  'vertices in x > %.10f' % (max_x - eps),
                                  'facet')

    field = Field.from_args('fu', nm.float64, 'vector', omega,
                            approx_order=2)

    return Struct(dim=dim, omega=omega, gamma1=gamma1, gamma2=gamma2,
                  field=field)

def test_term_evaluation(data):
    from sfepy.discrete import Integral, FieldVariable
    from sfepy.terms.terms import Term

    integral = Integral('i', order=3)

    u = FieldVariable('u', 'parameter', data.field,
                      primary_var_name='(set-to-None)')

    term = Term.new('ev_volume(u)', integral, data.omega, u=u)
    term *= 10.0

    term.setup()

    vol = term.evaluate()

    tst.report('volume: %.8f == 2000.0' % vol)
    assert nm.allclose(vol, 2000.0, rtol=1e-15, atol=0)

def test_term_arithmetics(data):
    from sfepy.discrete import FieldVariable, Integral
    from sfepy.terms.terms import Term

    integral = Integral('i', order=3)

    u = FieldVariable('u', 'parameter', data.field,
                      primary_var_name='(set-to-None)')

    t1 = Term.new('ev_volume(u)', integral, data.omega, u=u)
    t2 = Term.new('ev_volume(u)', integral, data.gamma1, u=u)

    expr = 2.2j * (t1 * 5.5 - 3j * t2) * 0.25

    ok = len(expr) == 2
    if not ok:
        tst.report('wrong expression length!')

    _ok = nm.allclose(expr[0].sign, 3.025j, rtol=1e-15, atol=0)
    if not _ok:
        tst.report('wrong sign of the first term!')
    ok = ok and _ok

    _ok = nm.allclose(expr[1].sign, 1.65, rtol=1e-15, atol=0)
    if not _ok:
        tst.report('wrong sign of the second term!')
    ok = ok and _ok

    assert ok

def test_variables(data):
    from sfepy.discrete import FieldVariable, Integral

    ok = True

    u = FieldVariable('u', 'parameter', data.field,
                      primary_var_name='(set-to-None)')
    u.set_constant(1.0)
    vec = u() # Nodal values.

    _ok = nm.allclose(vec, 1.0)
    tst.report('set constant:', _ok)
    ok = _ok and ok

    def fun(coors):
        val = nm.empty_like(coors)
        val[:, 0] = 2 * coors[:, 1] - coors[:, 0]
        val[:, 1] = coors[:, 0] + 3 * coors[:, 1]
        return val
    u.set_from_function(fun)

    coors = u.field.get_coor()
    eu = u.evaluate_at(coors)[..., 0]

    _ok = nm.allclose(eu, fun(coors), rtol=0.0, atol=1e-13)
    tst.report('set from function:', _ok)
    ok = _ok and ok

    integral = Integral('i', order=2)
    gu_qp = u.evaluate(mode='grad', integral=integral)

    # du_i/dx_j, i = column, j = row.
    gu = nm.array([[-1.,  1.],
                   [ 2.,  3.]])
    _ok = nm.allclose(gu_qp, gu[None, None, ...], rtol=0.0, atol=1e-13)
    tst.report('set from function - gradient:', _ok)
    ok = _ok and ok

    u_qp = gu_qp[..., :, :1]
    u.set_from_qp(u_qp, integral)
    vu = u()

    _ok = (nm.allclose(vu[::2], -1, rtol=0.0, atol=1e-13) and
           nm.allclose(vu[1::2], 2, rtol=0.0, atol=1e-13))
    tst.report('set from qp:', _ok)
    ok = _ok and ok

    assert ok

def test_solving(data, output_dir):
    from sfepy.base.base import IndexedStruct
    from sfepy.discrete import (FieldVariable, Material, Problem, Function,
                                Equation, Equations, Integral)
    from sfepy.discrete.conditions import Conditions, EssentialBC
    from sfepy.terms import Term
    from sfepy.solvers.ls import ScipyDirect
    from sfepy.solvers.nls import Newton
    from sfepy.mechanics.matcoefs import stiffness_from_lame

    u = FieldVariable('u', 'unknown', data.field)
    v = FieldVariable('v', 'test', data.field, primary_var_name='u')

    m = Material('m', D=stiffness_from_lame(data.dim, 1.0, 1.0))
    f = Material('f', val=[[0.02], [0.01]])

    bc_fun = Function('fix_u_fun', fix_u_fun,
                      extra_args={'extra_arg' : 'hello'})

    fix_u = EssentialBC('fix_u', data.gamma1, {'u.all' : bc_fun})
    shift_u = EssentialBC('shift_u', data.gamma2, {'u.0' : 0.1})

    integral = Integral('i', order=3)

    t1 = Term.new('dw_lin_elastic(m.D, v, u)',
                  integral, data.omega, m=m, v=v, u=u)

    t2 = Term.new('dw_volume_lvf(f.val, v)', integral, data.omega, f=f, v=v)

    eq = Equation('balance', t1 + t2)
    eqs = Equations([eq])

    ls = ScipyDirect({})

    nls_status = IndexedStruct()
    nls = Newton({}, lin_solver=ls, status=nls_status)

    pb = Problem('elasticity', equations=eqs)
    ## pb.save_regions_as_groups('regions')

    pb.set_bcs(ebcs=Conditions([fix_u, shift_u]))
    pb.set_solver(nls)

    state = pb.solve(save_results=False)

    name = op.join(output_dir, 'test_high_level_solving.vtk')
    pb.save_state(name, state)

    ok = nls_status.condition == 0
    if not ok:
        tst.report('solver did not converge!')

    _ok = state.has_ebc()
    if not _ok:
        tst.report('EBCs violated!')

    ok = ok and _ok

    assert ok
