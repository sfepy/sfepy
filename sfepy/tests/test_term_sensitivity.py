import gc
import numpy as nm
import pytest

from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy import data_dir
import sfepy.base.testing as tst

filename_mesh = data_dir + '/meshes/3d/special/cube_sphere.mesh'

fields = {
    'scalar_field' : ('real', 'scalar', 'Omega', 1),
    'vector_field' : ('real', 'vector', 'Omega', 1),
}

variables = {
    'U1': ('parameter field', 'vector_field', '(set-to-None)'),
    'U2': ('parameter field', 'vector_field', '(set-to-None)'),
    'P1': ('parameter field', 'scalar_field', '(set-to-None)'),
    'P2': ('parameter field', 'scalar_field', '(set-to-None)'),
    'V': ('parameter field', 'vector_field', '(set-to-None)'),
}

regions = {
    'Omega': 'all',
    'Omega1': 'cells of group 1',
    'Omega2': 'cells of group 2',
    'Interface': ('r.Omega1 *v r.Omega2', 'facet', 'Omega1'),
}

integrals = {
    'i': 2,
}

piezo_g = nm.array([[0, 0, 0, 0, 10, 0],
                    [0, 0, 0, 0, 0, 10],
                    [-5, -5, 15, 0, 0, 0]])

materials = {
    'mat': ({
        'D': {
            'Omega1': stiffness_from_youngpoisson(3, 100, 0.3),
            'Omega2': stiffness_from_youngpoisson(3, 1, 0.49),
        },
        'K': {'Omega1': nm.eye(3), 'Omega2': nm.eye(3) * 1e-3},
        'c': {'Omega1': 1, 'Omega2': 1000},
        'g': {'Omega1': piezo_g * 1e2, 'Omega2': piezo_g},
        's': {
            'Omega1': nm.array([[1, 2, 3, 0, 0, 0]]).T,
            'Omega2': nm.array([[8, 4, 2, 0, 0, 0]]).T,
        },
    },),
}

equations = {}

test_terms = [
    ('ev_sd_lin_elastic', 'dw_lin_elastic', 'Omega', 'mat.D', 'U1', 'U2'),
    ('de_sd_lin_elastic', 'dw_lin_elastic', 'Omega', 'mat.D', 'U1', 'U2'),
    ('ev_sd_dot', 'dw_dot', 'Omega', None, 'U1', 'U2'),
    ('de_sd_dot', 'de_dot', 'Omega', 'mat.K', 'U1', 'U2'),
    ('de_sd_dot', 'de_dot', 'Omega', 'mat.c', 'P1', 'P2'),
    ('ev_sd_diffusion', 'dw_diffusion', 'Omega', 'mat.K', 'P1', 'P2'),
    ('de_sd_diffusion', 'de_diffusion', 'Omega', 'mat.K', 'P1', 'P2'),
    ('ev_sd_div', 'dw_stokes', 'Omega', None, 'U1', 'P1'),
    ('de_sd_stokes', 'de_stokes', 'Omega', None, 'U1', 'P1'),
    ('ev_sd_div_grad', 'dw_div_grad', 'Omega', 'mat.c', 'U1', 'U2'),
    ('de_sd_div_grad', 'de_div_grad', 'Omega', 'mat.c', 'U1', 'U2'),
    ('ev_sd_piezo_coupling', 'dw_piezo_coupling', 'Omega', 'mat.g', 'U1', 'P1'),
    ('de_sd_piezo_coupling', 'dw_piezo_coupling', 'Omega', 'mat.g', 'U1', 'P1'),
    ('ev_sd_convect', 'de_convect', 'Omega', None, 'U1', 'U2'),
    ('ev_sd_surface_ltr', 'dw_surface_ltr', 'Interface', 'mat.s', 'U1', None),
    ('de_sd_surface_ltr', 'de_surface_ltr', 'Interface', 'mat.s', 'U1', None),
    ('de_sd_surface_ltr', 'de_surface_ltr', 'Interface', 'mat.c', 'U1', None),
    ('de_sd_v_dot_grad_s', 'dw_v_dot_grad_s', 'Omega', 'mat.c', 'U1', 'P1'),
]


def modify_mesh(val, spbox, dv_mode, cp_pos):
    cpoints0 = spbox.get_control_points().copy()
    for pts, dir in dv_mode:
        for pt in pts:
            spbox.move_control_point(cp_pos[pt], nm.array(dir) * val)

    new_coors = spbox.evaluate()
    spbox.set_control_points(cpoints0)

    return new_coors


@pytest.fixture(scope='module')
def problem():
    import sys
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    problem = Problem.from_conf(conf, init_equations=False)
    return problem


@pytest.mark.slow
def test_sensitivity(problem):
    from sfepy.discrete import Variables
    from sfepy.mesh.splinebox import SplineBox

    tolerance = 1e-4
    ok = True

    variables = Variables.from_conf(problem.conf.variables, problem.fields)

    for var_name in variables.names:
        var = variables[var_name]
        n_dof = var.field.n_nod * var.field.shape[0]
        aux = nm.arange(n_dof, dtype=nm.float64)
        var.set_data(aux)

    mesh = problem.domain.mesh
    bbox = nm.array(mesh.get_bounding_box()).T
    spbox = SplineBox(bbox, mesh.coors)

    dvel_modes = [
        # expand inner cylinder, no volume change
        [([20, 21, 22, 23], (-1, -1, 0)),
         ([24, 25, 26, 27], (-1, 1, 0)),
         ([36, 37, 38, 39], (1, -1, 0)),
         ([40, 41, 42, 43], (1, 1, 0))],
        # volume change
        [(range(16, 32), (0.2, 0, 0)),
         (range(32, 48), (0.4, 0, 0)),
         (range(48, 52), (0.6, 0.2, 0.2)),
         (range(52, 56), (0.8, 0.2, 0.3)),
         (range(56, 60), (1.0, 0.2, 0.4)),
         (range(60, 64), (1.2, 0.2, 0.5))],
    ]

    r4 = range(4)
    cp_pos = {i*16 + j*4 + k: (i, j, k)
        for k in r4 for j in r4 for i in r4}

    # compute design velocities
    dvels = []
    for dv_mode in dvel_modes:
        dvel = 0
        for pts, dir in dv_mode:
            for pt in pts:
                dvel += spbox.evaluate_derivative(cp_pos[pt], dir)
        dvels.append(dvel)

    for tname_sa, tname, rname, mat, var1, var2 in test_terms:
        args = [] if mat is None else [mat]
        args += [var1] if var2 is None else [var1, var2]
        term = '%s.i.%s(%s)' % (tname, rname, ', '.join(args))
        term_sa = '%s.i.%s(%s)' % (tname_sa, rname, ', '.join(args + ['V']))

        val = problem.evaluate(term, var_dict=variables.as_dict())
        tst.report('%s: %s' % (tname, val))

        dt = 1e-6
        for ii, dvel in enumerate(dvels):
            val = problem.evaluate(term, var_dict=variables.as_dict())
            variables['V'].set_data(dvel)
            val_sa = problem.evaluate(term_sa, var_dict=variables.as_dict())
            tst.report('%s - mesh_velocity mode %d' % (tname_sa, ii))
            # mesh perturbation +
            new_coors = modify_mesh(dt/2., spbox, dvel_modes[ii], cp_pos)
            problem.set_mesh_coors(new_coors, update_fields=True)
            val1 = problem.evaluate(term, var_dict=variables.as_dict())

            # mesh perturbation -
            new_coors = modify_mesh(-dt/2., spbox, dvel_modes[ii], cp_pos)
            problem.set_mesh_coors(new_coors, update_fields=True)
            val2 = problem.evaluate(term, var_dict=variables.as_dict())

            val_fd = (val1 - val2) / dt
            err = nm.abs(val_sa - val_fd) / nm.linalg.norm(val_sa)
            tst.report('term:               %s' % val)
            tst.report('sensitivity term:   %s' % val_sa)
            tst.report('finite differences: %s' % val_fd)
            tst.report('relative error:     %s' % err)

            _ok = err < tolerance

            ok = ok and _ok

            gc.collect()

    assert ok
