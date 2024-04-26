from functools import partial

import numpy as nm
import pytest

from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh
import sfepy.mechanics.matcoefs as mc
import sfepy.base.testing as tst

def define(t1=15e-6, dt=1e-6, dims=(0.1, 0.02, 0.005), shape=(11, 3, 3),
           young=70e9, poisson=0.3, density=2700, mass_lumping='row_sum',
           mass_beta=0.2):

    def mesh_hook(mesh, mode):
        """
        Generate the block mesh.
        """
        if mode == 'read':
            mesh = gen_block_mesh(dims, shape, 0.5 * nm.array(dims),
                                  name='user_block', verbose=False)
            return mesh

        elif mode == 'write':
            pass

    filename_mesh = UserMeshIO(mesh_hook)
    dim = len(dims)

    def get_sensor(coors, domain=None):
        ii = coors.argmax(axis=0)
        return ii[:1]

    regions = {
        'Omega' : 'all',
        'Left' : ('vertices in (x < 1e-12)', 'facet'),
        'Sensor' : ('vertices by get_sensor', 'vertex'),
    }

    materials = {
        'solid' : ({
            'D': mc.stiffness_from_youngpoisson(
                dim=dim, young=young, poisson=poisson, plane='strain'
            ),
            'rho': density,
            '.lumping' : mass_lumping,
            '.beta' : mass_beta,
         },),
    }

    fields = {
        'displacement': ('real', 'vector', 'Omega', 2),
    }

    integrals = {
        'i' : 4,
    }

    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'du' : ('unknown field', 'displacement', 1),
        'ddu' : ('unknown field', 'displacement', 2),
        'v' : ('test field', 'displacement', 'u'),
        'dv' : ('test field', 'displacement', 'du'),
        'ddv' : ('test field', 'displacement', 'ddu'),
    }
    var_names = {'u' : 'u', 'du' : 'du', 'ddu' : 'ddu'}

    ebcs = {
        'fix' : ('Left', {'u.all' : 0.0, 'du.all' : 0.0, 'ddu.all' : 0.0}),
    }

    def get_ic(coors, ic, mode='u'):
        val = nm.zeros_like(coors)
        if mode == 'u':
            val[:, 0] = 0.0

        elif mode == 'du':
            xmax = coors[:, 0].max()
            val[:, -1] = nm.where((coors[:, 0] > (xmax - 1e-12)), 1.0, 0.0)

        return val

    functions = {
        'get_sensor' : (get_sensor,),
        'get_ic_u' : (get_ic,),
        'get_ic_du' : (lambda coor, ic: get_ic(coor, None, mode='du'),),
    }

    ics = {
        'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
    }

    equations = {
        'balance_of_forces' :
        """dw_dot.i.Omega(solid.rho, ddv, ddu)
         + dw_zero.i.Omega(dv, du)
         + dw_lin_elastic.i.Omega(solid.D, v, u) = 0""",
    }

    solvers = {
        'ls' : ('ls.auto_direct', {
            # Reuse the factorized linear system from the first time step.
            'use_presolve' : True,
            # Speed up the above by omitting the matrix digest check used
            # normally for verification that the current matrix corresponds to
            # the factorized matrix stored in the solver instance. Use with
            # care!
            'use_mtx_digest' : False,
        }),
        'lsrmm' : ('ls.rmm', {
            'rmm_term' : """de_mass.i.Omega(solid.rho, solid.lumping,
                                            solid.beta, ddv, ddu)""",
            'debug' : False,
        }),
        'newton' : ('nls.newton', {
            'i_max'      : 1,
            'eps_a'      : 1e-6,
            'eps_r'      : 1e-6,
        }),
        'tsvv' : ('ts.velocity_verlet', {
            # Explicit method.
            't0' : 0.0,
            't1' : t1,
            'dt' : 0.1 * dt,
            'n_step' : None,

            'is_linear'  : True,

            'var_names' : var_names,
            'verbose' : 1,
        }),
        'tscd' : ('ts.central_difference', {
            # Explicit method.
            't0' : 0.0,
            't1' : t1,
            'dt' : 0.1 * dt,
            'n_step' : None,

            'is_linear'  : True,

            'var_names' : var_names,
            'verbose' : 1,
        }),
        'tsn' : ('ts.newmark', {
            't0' : 0.0,
            't1' : t1,
            'dt' : dt,
            'n_step' : None,

            'is_linear'  : True,

            'beta' : 0.25,
            'gamma' : 0.5,

            'var_names' : var_names,
            'verbose' : 1,
        }),
        'tsga' : ('ts.generalized_alpha', {
            't0' : 0.0,
            't1' : t1,
            'dt' : dt,
            'n_step' : None,

            'is_linear'  : True,

            'rho_inf' : 0.95,
            'alpha_m' : None,
            'alpha_f' : None,
            'beta' : None,
            'gamma' : None,

            'var_names' : var_names,
            'verbose' : 1,
        }),
        'tsb' : ('ts.bathe', {
            't0' : 0.0,
            't1' : t1,
            'dt' : 0.5 * dt,
            'n_step' : None,

            'is_linear'  : True,

            'var_names' : var_names,
            'verbose' : 1,
        }),
        'tscedb' : ('tsc.ed_basic', {
            'eps_r' : (1e-3, 1e-1),
            'eps_a' : (1e-6, 1e-2),
            'fmin' : 0.3,
            'fmax' : 2.5,
            'fsafety' : 0.85,
        }),
        'tscedl' : ('tsc.ed_linear', {
            'eps_r' : (1e-3, 1e-1),
            'eps_a' : (1e-6, 1e-2),
            'fmin' : 0.3,
            'fmax' : 2.5,
            'fsafety' : 0.85,
            'fred' : 0.9,
            'inc_wait' : 10,
            'min_finc' : 1.7,
        }),
        'tscedpid' : ('tsc.ed_pid', {
            'eps_r' : (1e-3, 1e-1),
            'eps_a' : (1e-6, 1e-2),
            'fmin' : 0.3,
            'fmax' : 2.5,
            'fsafety' : 0.85,
            'pcoef' : 0.4,
            'icoef' : 0.3,
            'dcoef' : 0,
        }),
    }

    options = {
        'ts' : 'tsn',
        'nls' : 'newton',
        'ls' : 'ls',

        'save_times' : 31,

        'active_only' : False,

        'output_format' : 'h5',
    }

    return locals()

@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    define_dict = define()
    conf = ProblemConf.from_dict(define_dict, sys.modules[__name__])

    pb = Problem.from_conf(conf)
    pb.update_materials()

    # Get full size matrices for energy calculations.
    pb.init_solvers()
    tss = pb.solver
    ebcs = pb.conf.ebcs
    pb.time_update(ebcs={})
    tss.constant_matrices = None
    pb.Mf, Cf, pb.Kf = tss.get_matrices(tss.nls, pb.set_default_state()())

    # Restore EBCs.
    pb.time_update(ebcs=ebcs)

    return pb

def _list_solvers(confs, kind='ts'):
    d = [val for val in confs.values() if val.kind.startswith(kind+'.')]
    d.sort(key=lambda a: a.name)

    return d

@pytest.mark.slow
def test_ed_solvers(problem, output_dir):
    from scipy.integrate import simpson
    from sfepy.base.base import IndexedStruct

    tss_confs = _list_solvers(problem.solver_confs)
    tsc_confs = _list_solvers(problem.solver_confs, kind='tsc')

    vu = problem.get_variables()['u']
    sensor = problem.domain.regions['Sensor']
    isens = 3 * vu.field.get_dofs_in_region(sensor)[0] + 2

    def store_ths(pb, ts, variables, ths):
        sp = variables.get_state_parts()
        u1, v1, a1 = sp['u'], sp['du'], sp['ddu']

        e_u = 0.5 * u1 @ pb.Kf @ u1
        e_t = 0.5 * v1 @ pb.Mf @ v1
        ths.append((ts.time, u1[isens], v1[isens], a1[isens],
                    e_u, e_t, e_u + e_t))

    all_ths = []
    stats = []
    t1s = []
    for tsc_conf in [None] + tsc_confs:
        problem.tsc_conf = None
        for tss_conf in tss_confs:
            status = IndexedStruct()
            problem.init_solvers(tsc_conf=tsc_conf, ts_conf=tss_conf,
                                 status=status, force=True)
            ths = []
            problem.solve(status=status, save_results=False,
                          step_hook=partial(store_ths, ths=ths))
            all_ths.append(nm.array(ths))
            stats.append((problem.solver.tsc.conf.kind, tss_conf.kind,
                          status.n_step, status.time))
            t1s.append(problem.solver.ts.time)

    kinds = [val[0:2] for val in stats]
    stats.sort(key=lambda x: x[-1])
    tst.report('solution times / numbers of time steps:')
    for row in stats:
        tst.report('%.2f [s] / % 4d' % (row[3], row[2]), ':', row[:2])

    # import matplotlib.pyplot as plt
    # for ii, ths in enumerate(all_ths):
    #     fig, ax = plt.subplots()
    #     ax.plot(ths[:,0], ths[:,4])
    #     ax.plot(ths[:,0], ths[:,5])
    #     ax.plot(ths[:,0], ths[:,6])
    #     ax.set_title(kinds[ii])
    # plt.show()

    all_iths = nm.array(
        [[simpson(val, x=ths[:, 0]) for val in ths.T[1:]] for ths in all_ths]
    )

    tst.report('status, solver: time integral of (u, v, a, e_u, e_t, e_u-e_t)')
    iths_ref = all_iths[0]
    e0 = all_ths[0][0, -1]
    e_rtols = {
        ('tsc.fixed', 'ts.bathe') : 1e-1,
        ('tsc.fixed', 'ts.generalized_alpha') : 1e-2,
        ('tsc.fixed', 'ts.newmark') : 1e-12,
        ('tsc.fixed', 'ts.central_difference') : 1e-2,
        ('tsc.fixed', 'ts.velocity_verlet') : 1e-2,
        ('tsc.ed_basic', 'ts.bathe') : 1e-1,
        ('tsc.ed_basic', 'ts.generalized_alpha') : 1e-3,
        ('tsc.ed_basic', 'ts.newmark') : 1e-12,
        ('tsc.ed_basic', 'ts.central_difference') : 2e-2,
        ('tsc.ed_basic', 'ts.velocity_verlet') : 2e-2,
        ('tsc.ed_linear', 'ts.bathe') : 1e-1,
        ('tsc.ed_linear', 'ts.generalized_alpha') : 1e-3,
        ('tsc.ed_linear', 'ts.newmark') : 1e-12,
        ('tsc.ed_linear', 'ts.central_difference') : 2e-2,
        ('tsc.ed_linear', 'ts.velocity_verlet') : 2e-2,
        ('tsc.ed_pid', 'ts.bathe') : 1e-2,
        ('tsc.ed_pid', 'ts.generalized_alpha') : 1e-4,
        ('tsc.ed_pid', 'ts.newmark') : 1e-12,
        ('tsc.ed_pid', 'ts.central_difference') : 1e-2,
        ('tsc.ed_pid', 'ts.velocity_verlet') : 1e-2,
    }
    ok = True
    for ii, iths in enumerate(all_iths):
        ienergy = e0 * t1s[ii]
        _ok = ((abs(iths[0] - iths_ref[0]) < 2e-9) and
               nm.isclose(iths[-1], ienergy, atol=0, rtol=e_rtols[kinds[ii]]))
        tst.report(('%d % 20s:' + (6 * ' % .2e'))
                   % ((_ok, kinds[ii]) + tuple(iths)))
        ok = _ok and ok

    assert ok
    assert nm.isclose(e0, 1.8e-4, atol=0, rtol=1e-12)

def test_rmm_solver(problem, output_dir):
    from sfepy.base.base import IndexedStruct

    ls_conf = problem.solver_confs['ls']
    lsr_conf = problem.solver_confs['lsrmm']
    tss_conf = problem.solver_confs['tscd']

    vu = problem.get_variables()['u']
    sensor = problem.domain.regions['Sensor']
    isens = 3 * vu.field.get_dofs_in_region(sensor)[0] + 2

    def store_ths(pb, ts, variables, ths):
        sp = variables.get_state_parts()
        u1, v1 = sp['u'], sp['du']

        e_u = 0.5 * u1 @ pb.Kf @ u1
        e_t = 0.5 * v1 @ pb.Mf @ v1
        ths.append((ts.time, u1[isens], v1[isens], e_u, e_t))

    problem.tsc_conf = None

    status = IndexedStruct()
    problem.init_solvers(ls_conf=ls_conf, ts_conf=tss_conf, status=status,
                         force=True)
    ths = []
    problem.solve(status=status, save_results=False,
                  step_hook=partial(store_ths, ths=ths))
    ths = nm.array(ths)

    statusr = IndexedStruct()
    problem.init_solvers(ls_conf=lsr_conf, ts_conf=tss_conf, status=statusr,
                         force=True)
    thsr = []
    problem.solve(status=statusr, save_results=False,
                  step_hook=partial(store_ths, ths=thsr))
    thsr = nm.array(thsr)

    tst.report(f'solution times: CMM: {status.time}, RMM: {statusr.time}')

    ratio = abs(ths[:,2].max() / ths[:,1].max())

    # import matplotlib.pyplot as plt
    # colors = plt.cm.tab10.colors
    # fig, ax = plt.subplots()
    # ax.plot(ths[:,0], ths[:,3], color=colors[0], ls='-')
    # ax.plot(ths[:,0], ths[:,4], color=colors[1], ls='-')
    # ax.plot(thsr[:,0], thsr[:,3], color=colors[0], ls='--')
    # ax.plot(thsr[:,0], thsr[:,4], color=colors[1], ls='--')
    # fig, ax = plt.subplots()
    # ax.plot(ths[:,0], ratio * ths[:,1], color=colors[0], ls='-')
    # ax.plot(ths[:,0], ths[:,2], color=colors[1], ls='-')
    # ax.plot(thsr[:,0], ratio * thsr[:,1], color=colors[0], ls='--')
    # ax.plot(thsr[:,0], thsr[:,2], color=colors[1], ls='--')
    # plt.show()

    dt = ths[1, 0] - ths[0, 0]
    ierrs = nm.linalg.norm(ths[:, 1:] - thsr[:, 1:], axis=0) * dt

    assert ierrs[0] < 1e-13
    assert ierrs[1] < 3 * ratio * 1e-13
    assert ierrs[2] < 1e-11
    assert ierrs[3] < 1e-11

def test_active_only(output_dir):
    """
    Note: with tsc the results would differ, as eval_scaled_norm() depends on
    the vector length.
    """
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    define_dict = define(dims=(0.1, 0.02), shape=(3, 3))
    conf = ProblemConf.from_dict(define_dict, sys.modules[__name__])

    conf.options.active_only = False
    pb = Problem.from_conf(conf)
    pb.tsc_conf = None
    variables_f = pb.solve(save_results=False)

    conf.options.active_only = True
    pb = Problem.from_conf(conf)
    pb.tsc_conf = None
    variables_t = pb.solve(save_results=False)

    assert nm.allclose(variables_f(), variables_t(), atol=1e-7, rtol=0)
