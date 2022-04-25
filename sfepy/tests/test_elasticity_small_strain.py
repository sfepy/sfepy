import pytest

from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_lame
import sfepy.base.testing as tst

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
    sym = (dim + 1) * dim // 2
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
        'Dijkl' : get_pars( 3, True ),
        'D' : stiffness_from_lame(3, get_pars(3)[0], get_pars(3)[1]),
        'lam' : get_pars(3)[0],
        'mu' : get_pars(3)[1],
    }
}

material_2 = {
    'name' : 'spring',
    'values' : {
        '.stiffness' : 1e0,
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
    'select' : 'vertices in (z < -0.499)',
    'kind' : 'facet',
}
region_2 = {
    'name' : 'Top',
    'select' : 'vertices in (z > 0.499)',
    'kind' : 'facet',
}

ebc_1 = {
    'name' : 'Load',
    'region' : 'Top',
    'dofs' : {'u.2' : 0.1},
}

integral_1 = {
    'name' : 'i',
    'order' : 2,
}

equations_getpars = {
    'balance_of_forces' :
    """dw_lin_elastic.i.Omega(solid.Dijkl, v, u)
     = dw_point_lspring.i.Bottom(spring.stiffness, v, u)""",
}

equations_matcoefs = {
    'balance_of_forces' :
    """dw_lin_elastic.i.Omega(solid.D, v, u)
     = dw_point_lspring.i.Bottom(spring.stiffness, v, u)""",
}

equations_iso = {
    'balance_of_forces' :
    """dw_lin_elastic_iso.i.Omega(solid.lam, solid.mu, v, u)
     = dw_point_lspring.i.Bottom(spring.stiffness, v, u)""",
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
}

@pytest.fixture(scope='module')
def solutions(output_dir):
    import sys
    from sfepy.applications import solve_pde
    from sfepy.base.base import IndexedStruct, Struct
    from sfepy.base.conf import ProblemConf
    import os.path as op

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])

    solutions = []
    for ii, approx_order in enumerate(all_your_bases):
        sols = Struct()

        fname = filename_meshes[ii]

        conf.filename_mesh = fname
        fields = {'field_1' : {
                      'name' : '3_displacement',
                      'dtype' : 'real',
                      'shape' : (3,),
                      'region' : 'Omega',
                      'approx_order' : approx_order,
                }
        }
        conf.edit('fields', fields)
        tst.report('mesh: %s, base: %s' % (fname, approx_order))
        status = IndexedStruct()

        tst.report('getpars')
        conf.equations = conf.equations_getpars
        problem, state1 = solve_pde(conf, status=status,
                                    save_results=False)
        converged = status.nls_status.condition == 0
        sols.getpars = (converged, state1().copy())
        tst.report('converged: %s' % converged)

        tst.report('matcoefs')
        conf.equations = conf.equations_matcoefs
        problem, state2 = solve_pde(conf, status=status,
                                    save_results=False)
        converged = status.nls_status.condition == 0
        sols.matcoefs = (converged, state2().copy())
        tst.report('converged: %s' % converged)

        tst.report('iso')
        conf.equations = conf.equations_iso
        problem, state3 = solve_pde(conf, status=status,
                                    save_results=False)
        converged = status.nls_status.condition == 0
        sols.iso = (converged, state3().copy())
        tst.report('converged: %s' % converged)

        solutions.append(sols)

        name = op.join(output_dir,
                       '_'.join(('test_elasticity_small_strain',
                                 op.splitext(op.basename(fname))[0],
                                 '%d' % approx_order))
                       + '.vtk')
        problem.save_state(name, state1)

    return solutions

def test_converged(solutions):
    for sols in solutions:
        assert sols.getpars[0]
        assert sols.matcoefs[0]
        assert sols.iso[0]

def test_linear_terms(solutions):
    for sols in solutions:
        assert tst.compare_vectors(sols.getpars[1], sols.matcoefs[1],
                                   label1='getpars',
                                   label2='matcoefs')
        assert tst.compare_vectors(sols.getpars[1], sols.iso[1],
                                   label1='getpars',
                                   label2='iso')
