import os.path as op

import numpy as nm
import pytest

import sfepy.base.testing as tst
from sfepy import data_dir

def init_vec(variables):
    return nm.random.rand(variables.di.n_dof_total)

def check_vec(vec, ii, ok, conds, variables):
    from sfepy.discrete.common.dof_info import expand_nodes_to_equations

    for var_name, var_conds in conds.group_by_variables().items():
        var = variables[var_name]
        for cond in var_conds:
            cond.canonize_dof_names(var.dofs)
            tst.report('%d: %s %s: %s %s'
                       % (ii, var.name,
                          cond.name, cond.region.name, cond.dofs[0]))
            nods = var.field.get_dofs_in_region(cond.region)
            eq = expand_nodes_to_equations(nods, cond.dofs[0], var.dofs)

            off = variables.di.indx[var_name].start
            n_nod = len(nods)
            for cdof, dof_name in enumerate(cond.dofs[0]):
                idof = var.dofs.index(dof_name)
                eqs = eq[n_nod * cdof : n_nod * (cdof + 1)]

                _ok = nm.allclose(vec[off + eqs], idof,
                                  atol=1e-14, rtol=0.0)
                if not _ok:
                    tst.report(' %s: failed! (all of %s == %f)'
                               % (dof_name, vec[off + eqs], idof))
                ok = ok and _ok

    return ok

@pytest.fixture(scope='module')
def data():
    from sfepy.base.base import Struct
    from sfepy.discrete import FieldVariable, Variables, Problem
    from sfepy.discrete.fem import Mesh, FEDomain, Field

    mesh = Mesh.from_file(data_dir + '/meshes/2d/square_unit_tri.mesh')
    domain = FEDomain('domain', mesh)

    omega = domain.create_region('Omega', 'all')
    domain.create_region('Left',
                         'vertices in (x < -0.499)',
                         'facet')
    domain.create_region('LeftStrip',
                         'vertices in (x < -0.499)'
                         ' & (y > -0.199) & (y < 0.199)',
                         'facet')
    domain.create_region('LeftFix',
                         'r.Left -v r.LeftStrip',
                         'facet')
    domain.create_region('Right',
                         'vertices in (x > 0.499)',
                         'facet')
    domain.create_region('RightStrip',
                         'vertices in (x > 0.499)'
                         ' & (y > -0.199) & (y < 0.199)',
                         'facet')
    domain.create_region('RightFix',
                         'r.Right -v r.RightStrip',
                         'facet')

    fu = Field.from_args('fu', nm.float64, 'vector', omega, approx_order=2)
    u = FieldVariable('u', 'unknown', fu)

    fp = Field.from_args('fp', nm.float64, 'scalar', omega, approx_order=2)
    p = FieldVariable('p', 'unknown', fp)

    pb = Problem('test', domain=domain, fields=[fu, fp],
                 auto_conf=False)

    return Struct(problem=pb, variables=Variables([u, p]))

def test_ics(data):
    from sfepy.discrete.conditions import Conditions, InitialCondition

    variables = data.variables
    omega = data.problem.domain.regions['Omega']

    all_ics = []
    all_ics.append(InitialCondition('ic0', omega,
                                    {'p.all' : 0.0}))
    all_ics.append(InitialCondition('ic1', omega,
                                    {'u.1' : 1.0}))
    all_ics.append(InitialCondition('ic2', omega,
                                    {'u.all' : nm.array([0.0, 1.0])}))
    all_ics.append(InitialCondition('ic3', omega,
                                    {'p.0' : 0.0,
                                     'u.0' : 0.0, 'u.1' : 1.0}))

    ok = True
    for ii, ics in enumerate(all_ics):
        if not isinstance(ics, list): ics = [ics]

        ics = Conditions(ics)
        variables.setup_initial_conditions(ics, functions=None)

        vec = init_vec(variables)
        variables.apply_ic(vec)

        ok = check_vec(vec, ii, ok, ics, variables)

    assert ok

def test_ebcs(data):
    from sfepy.discrete.conditions import Conditions, EssentialBC

    variables = data.variables
    regions = data.problem.domain.regions

    all_ebcs = []
    all_ebcs.append(EssentialBC('fix_u1', regions['LeftFix'],
                                {'u.all' : nm.array([0.0, 1.0])}))
    all_ebcs.append(EssentialBC('fix_u2', regions['LeftStrip'],
                                {'u.0' : 0.0, 'u.1' : 1.0}))
    all_ebcs.append(EssentialBC('fix_p1', regions['RightFix'],
                                {'p.all' : 0.0}))
    all_ebcs.append(EssentialBC('fix_p2', regions['RightStrip'],
                                {'p.0' : 0.0}))
    all_ebcs.append([EssentialBC('fix_p3', regions['Right'],
                                 {'p.0' : 0.0}),
                     EssentialBC('fix_u3', regions['Left'],
                                 {'u.0' : 0.0, 'u.1' : 1.0})])

    ok = True
    for ii, bcs in enumerate(all_ebcs):
        if not isinstance(bcs, list): bcs = [bcs]

        ebcs = Conditions(bcs)
        variables.equation_mapping(ebcs=ebcs, epbcs=None,
                                   ts=None, functions=None)
        vec = init_vec(variables)
        variables.apply_ebc(vec)

        ok = check_vec(vec, ii, ok, ebcs, variables)

    assert ok

def test_epbcs(data):
    from sfepy.discrete import Function, Functions
    from sfepy.discrete.conditions import Conditions, PeriodicBC
    from sfepy.discrete.common.dof_info import expand_nodes_to_equations
    from sfepy.discrete.fem.periodic import match_y_line

    variables = data.variables
    regions = data.problem.domain.regions

    match_y_line = Function('match_y_line', match_y_line)
    pbc = PeriodicBC('pbc', [regions['LeftStrip'], regions['RightStrip']],
                     {'u.[1,0]' : 'u.[0,1]'}, match='match_y_line')

    functions = Functions([match_y_line])

    epbcs = Conditions([pbc])
    variables.equation_mapping(ebcs=None, epbcs=epbcs,
                               ts=None, functions=functions)

    vec = init_vec(variables)
    variables.apply_ebc(vec)

    var = variables['u']
    var_bcs = epbcs.group_by_variables()['u']
    bc = var_bcs['pbc']
    bc.canonize_dof_names(var.dofs)

    nods0 = var.field.get_dofs_in_region(bc.regions[0])
    nods1 = var.field.get_dofs_in_region(bc.regions[1])

    coors0 = var.field.get_coor(nods0)
    coors1 = var.field.get_coor(nods1)
    i0, i1 = match_y_line(coors0, coors1)

    eq0 = expand_nodes_to_equations(nods0[i0], bc.dofs[0], var.dofs)
    eq1 = expand_nodes_to_equations(nods1[i1], bc.dofs[1], var.dofs)

    ok = True

    _ok = variables.has_ebc(vec)
    if not _ok:
        tst.report('EPBCs were not applied correctly!')
    ok = ok and _ok

    _ok = len(nm.setdiff1d(eq0, var.eq_map.master)) == 0
    if not _ok:
        tst.report('master equations mismatch! (set(%s) == set(%s))'
                   % (eq0, var.eq_map.master))
    ok = ok and _ok

    _ok = len(nm.setdiff1d(eq1, var.eq_map.slave)) == 0
    if not _ok:
        tst.report('slave equations mismatch! (set(%s) == set(%s))'
                   % (eq1, var.eq_map.slave))
    ok = ok and _ok

    off = variables.di.indx['u'].start
    _ok = nm.allclose(vec[off + eq0], vec[off + eq1], atol=1e-14, rtol=0.0)
    if not _ok:
        tst.report('periodicity test failed! (%s == %s)'
                   % (vec[off + eq0], vec[off + eq1]))
    ok = ok and _ok

    assert ok

def test_save_ebc(data, output_dir):
    from sfepy.discrete import (FieldVariable, Integral,
                                Equation, Equations, Problem)
    from sfepy.discrete.conditions import Conditions, EssentialBC
    from sfepy.terms import Term

    name = op.join(output_dir,
                   op.splitext(op.basename(__file__))[0])

    integral = Integral('i', order=1)

    u = data.variables['u']
    v = FieldVariable('v', 'test', u.field, primary_var_name='u')

    p = data.variables['p']
    q = FieldVariable('q', 'test', p.field, primary_var_name='p')

    regions = data.problem.domain.regions
    omega = regions['Omega']

    # Problem.save_ebc() requires to have equations defined.
    t1 = Term.new('dw_lin_elastic(v, u)',
                  integral, omega, v=v, u=u)
    t2 = Term.new('dw_laplace(q, p)', integral, omega, q=q, p=p)
    eq = Equation('aux', t1 + t2)
    eqs = Equations([eq])

    pb = Problem('test', equations=eqs)

    all_ebcs = []
    all_ebcs.append(EssentialBC('fix_u1', regions['RightFix'],
                                {'u.all' : nm.array([0.0, 1.0])}))
    all_ebcs.append(EssentialBC('fix_u2', regions['LeftStrip'],
                                {'u.0' : 0.0, 'u.1' : 1.0}))
    all_ebcs.append(EssentialBC('fix_p1', regions['LeftFix'],
                                {'p.all' : 0.0}))
    all_ebcs.append(EssentialBC('fix_p2', regions['RightStrip'],
                                {'p.0' : 0.0}))
    ebcs = Conditions(all_ebcs)

    pb.time_update(ebcs=ebcs)

    pb.save_ebc(name + '_ebcs_f.vtk', ebcs=ebcs, force=True)
    pb.save_ebc(name + '_ebcs.vtk', ebcs=ebcs, default=-1, force=False)

    assert True
