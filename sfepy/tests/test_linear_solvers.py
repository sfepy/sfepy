import pytest

from sfepy import data_dir
import sfepy.base.testing as tst

filename_mesh = data_dir + '/meshes/3d/special/cube_cylinder.mesh'

if 0:
    from sfepy.discrete.fem.utils import refine_mesh
    refinement_level = 1
    filename_mesh = refine_mesh(filename_mesh, refinement_level)

material_2 = {
    'name' : 'coef',
    'values' : {'val' : 1.0},
}

field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.0001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.999)', 'facet'),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integral_1 = {
    'name' : 'i',
    'order' : 1,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega(coef.val, s, t) = 0"""
}

class DiagPC:
    """
    Diagonal (Jacobi) preconditioner.

    Equivalent to setting `'precond' : 'jacobi'`.
    """

    def setUp(self, pc):
        A = pc.getOperators()[0]

        self.idiag = 1.0 / A.getDiagonal()

    def apply(self, pc, x, y):
        y.pointwiseMult(x, self.idiag)

def setup_petsc_precond(mtx, problem):
    return DiagPC()

solvers = {
    'd00' : ('ls.scipy_direct',
             {}
    ),
    'd01' : ('ls.scipy_direct',
             {'method' : 'umfpack',
              'warn' : True,}
    ),
    'd02' : ('ls.scipy_direct',
             {'method' : 'superlu',
              'warn' : True,}
    ),
    'd10' : ('ls.mumps', {}),
    'd20' : ('ls.pypardiso', {}),
    'i00' : ('ls.pyamg',
             {'method' : 'ruge_stuben_solver',
              'accel' : 'cg',
              'eps_r'   : 1e-12,
              'method:max_levels' : 5,
              'solve:cycle' : 'V',}
    ),
    'i01' : ('ls.pyamg',
             {'method' : 'smoothed_aggregation_solver',
              'accel' : 'cg',
              'eps_r'   : 1e-12,}
    ),
    'i02' : ('ls.pyamg_krylov',
             {'method' : 'cg',
              'eps_r'   : 1e-12,
              'i_max' : 1000,}
    ),
    'i10' : ('ls.petsc',
             {'method' : 'cg', # ksp_type
              'precond' : 'none', # pc_type
              'eps_a' : 1e-12, # abstol
              'eps_r' : 1e-12, # rtol
              'i_max' : 1000,} # maxits
    ),
    'i11' : ('ls.petsc',
             {'method' : 'cg', # ksp_type
              'precond' : 'python', # just for output (unused)
              'setup_precond' : setup_petsc_precond, # user-defined pc
              'eps_a' : 1e-12, # abstol
              'eps_r' : 1e-12, # rtol
              'i_max' : 1000,} # maxits
    ),
    'i12' : ('ls.petsc',
             {'method' : 'cg', # ksp_type
              'precond' : 'jacobi', # pc_type
              'eps_a' : 1e-12, # abstol
              'eps_r' : 1e-12, # rtol
              'i_max' : 1000,} # maxits
    ),
    'i13' : ('ls.petsc',
             {'method' : 'cg', # ksp_type
              'precond' : 'icc', # pc_type
              'eps_a' : 1e-12, # abstol
              'eps_r' : 1e-12, # rtol
              'i_max' : 1000,} # maxits
    ),
    'i20' : ('ls.scipy_iterative',
             {'method' : 'cg',
              'i_max'   : 1000,
              'eps_a'   : 1e-12,
              'eps_r'   : 1e-12,}
    ),
    'i21' : ('ls.scipy_iterative',
             {'method' : 'bicgstab',
              'i_max'   : 1000,
              'eps_a'   : 1e-12,
              'eps_r'   : 1e-12,}
    ),
    'i22' : ('ls.scipy_iterative',
             {'method' : 'qmr',
              'i_max'   : 1000,
              'eps_a'   : 1e-12,
              'eps_r'   : 1e-12,}
    ),

    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}


options = {
    'nls' : 'newton',
}

can_fail = ['ls.pyamg', 'ls.pyamg_krylov', 'ls.petsc', 'ls.mumps',
            'ls.scipy_direct', 'ls.pypardiso']

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf)
    problem.time_update()

    return problem

def _list_linear_solvers(confs):
    d = []
    for key, val in confs.items():
        if val.kind.find('ls.') == 0:
            d.append(val)
    d.sort(key=lambda a: a.name)

    return d

def test_solvers(problem, output_dir):
    from sfepy.base.base import IndexedStruct
    import os.path as op

    solver_confs = _list_linear_solvers(problem.solver_confs)

    ok = True
    tt = []
    for solver_conf in solver_confs:
        method = solver_conf.get('method', '')
        precond = solver_conf.get('precond', '')
        name = ' '.join((solver_conf.name, solver_conf.kind,
                         method, precond)).rstrip()
        tst.report(name)
        tst.report('matrix size:', problem.mtx_a.shape)
        tst.report('        nnz:', problem.mtx_a.nnz)
        status = IndexedStruct()
        try:
            problem.init_solvers(status=status,
                                 ls_conf=solver_conf,
                                 force=True)
            state = problem.solve(save_results=False)
            failed = status.nls_status.condition != 0
        except Exception as aux:
            failed = True
            status = None
            exc = aux

        ok = ok and ((not failed) or (solver_conf.kind in can_fail))

        if status is not None:
            status = status.nls_status
            for kv in status.time_stats.items():
                tst.report('%10s: %7.2f [s]' % kv)
            tst.report('condition: %d, err0: %.3e, err: %.3e'
                       % (status.condition, status.err0, status.err))
            tt.append([name,
                       status.time_stats['solve'],
                       status.ls_n_iter,
                       status.err])

            aux = name.replace(' ', '_')
            fname = op.join(output_dir, 'test_linear_solvers_%s.vtk') % aux
            problem.save_state(fname, state)
        else:
            tst.report('solver failed:')
            tst.report(exc)
            tt.append([name, -1, 1e10, 1e10])

    tt.sort(key=lambda a: a[1])
    tst.report('solution times / numbers of iterations (residual norms):')
    for row in tt:
        tst.report('%.2f [s] / % 4d' % (row[1], row[2]),
                   '(%.3e)' % row[3], ':', row[0])

    assert ok

def test_ls_reuse(problem):
    import numpy as nm
    from sfepy.solvers import Solver

    problem.init_solvers(ls_conf=problem.solver_confs['d00'])
    nls = problem.get_nls()

    state0 = problem.get_initial_state()
    state0.apply_ebc()
    vec0 = state0.get_state(problem.active_only)

    problem.update_materials()

    rhs = nls.fun(vec0)
    mtx = nls.fun_grad(vec0)

    ok = True
    for name in ['i12', 'i01']:
        solver_conf = problem.solver_confs[name]
        method = solver_conf.get('method', '')
        precond = solver_conf.get('precond', '')
        name = ' '.join((solver_conf.name, solver_conf.kind,
                         method, precond)).rstrip()
        tst.report(name)
        try:
            ls = Solver.any_from_conf(solver_conf)

        except:
            tst.report('skipped!')
            continue

        conf = ls.conf.copy()
        conf.force_reuse = True

        sol00 = ls(rhs, mtx=mtx, conf=conf)
        digest00 = ls.mtx_digest

        sol0 = ls(rhs, mtx=mtx)
        digest0 = ls.mtx_digest

        sol1 = ls(rhs, mtx=2*mtx, conf=conf)
        digest1 = ls.mtx_digest

        sol2 = ls(rhs, mtx=2*mtx)
        digest2 = ls.mtx_digest
        ls(rhs, mtx=2*mtx)
        digest3 = ls.mtx_digest

        _ok = digest00 != digest0
        tst.report(digest00, '!=', digest0, ':', _ok); ok = ok and _ok
        _ok = digest0 == digest1
        tst.report(digest0, '==', digest1, ':', _ok); ok = ok and _ok
        _ok = digest1 != digest2
        tst.report(digest1, '!=', digest2, ':', _ok); ok = ok and _ok
        _ok = digest2[1] == digest3[1]
        tst.report(digest2[1], '==', digest3[1], ':', _ok); ok = ok and _ok
        _ok = nm.allclose(sol00, sol0, atol=1e-12, rtol=0.0)
        tst.report('sol00 == sol0:', _ok); ok = ok and _ok
        _ok = nm.allclose(sol0, sol1, atol=1e-12, rtol=0.0)
        tst.report('sol0 == sol1:', _ok); ok = ok and _ok
        _ok = nm.allclose(sol0, 2 * sol2, atol=1e-12, rtol=0.0)
        tst.report('sol0 == 2 * sol2:', _ok); ok = ok and _ok

    assert ok
