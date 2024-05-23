import numpy as nm
import pytest

from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.solvers import Solver
import sfepy.base.testing as tst

# Mesh dimensions.
dims = nm.array([1, 1])

# Mesh resolution.
shape = [21, 21]

def mesh_hook(mesh, mode):
    """
    Generate the block mesh.
    """
    if mode == 'read':
        mesh = gen_block_mesh(dims, shape, [0, 0], name='user_block',
                              verbose=False)
        return mesh

    elif mode == 'write':
        pass

filename_mesh = UserMeshIO(mesh_hook)

regions = {
    'Omega' : 'all',
    'Gamma' : ('vertices of surface', 'facet'),
}

fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

ebcs = {
    't' : ('Gamma', {'t.0' : 0.0}),
}

integrals = {
    'i' : 2,
}

equations = {
    'eq' : 'dw_laplace.i.Omega(s, t) = 0'
}

solvers = {
    'evp0' : ('eig.scipy', {
        'method' : 'eig',
    }),
    'evp0h' : ('eig.scipy', {
        'method' : 'eigh',
    }),
    'evp0s' : ('eig.scipy', {
        'method' : 'eigs',
        'maxiter' : 100,
    }),
    'evp0sh' : ('eig.scipy', {
        'method' : 'eigsh',
        'maxiter' : 100,
    }),
    'evp1' : ('eig.sgscipy', {}),
    'evp2' : ('eig.scipy_lobpcg', {
        'i_max' : 100,
        'largest' : False,
    }),
    'evp4' : ('eig.slepc', {
        'method' : 'arnoldi',
        'problem' : 'hep',
        'i_max' : 100,
        'eps' : 1e-10,
        'which' : 'smallest_real',
    }),
    'evp5' : ('eig.matlab', {
        'maxit' : 100,
        'tol' : 1e-10,
        'which' : 'sr',
    }),
    'evp6' : ('eig.primme', {
        'maxiter' : 200,
        'tol' : 1e-10,
        'which' : 'sa',
    }),
    'evp7' : ('eig.octave', {
        'maxit' : 100,
        'tol' : 1e-10,
        'which' : 'sm',
    }),
}

eigs_expected = [nm.array([0.04904454, 0.12170685, 0.12170685,
                           0.19257998, 0.24082108]),
                 []]

can_fail = ['eig.slepc', 'eig.matlab', 'eig.octave', 'eig.primme']
can_miss = ['evp0'] # Depending on scipy version, evp0 can miss an
                    # eigenvalue.

@pytest.fixture(scope='module')
def data():
    import sys
    from sfepy.base.base import Struct
    from sfepy.discrete import Problem
    from sfepy.base.conf import ProblemConf

    conf = ProblemConf.from_dict(globals(), sys.modules[__name__])
    pb = Problem.from_conf(conf, init_solvers=False)
    pb.time_update()
    mtx = pb.equations.evaluate(mode='weak', dw_mode='matrix',
                                asm_obj=pb.mtx_a)

    return Struct(conf=conf, mtx=mtx)

def _list_eigenvalue_solvers(confs):
    d = []
    for key, val in confs.items():
        if val.kind.find('eig.') == 0:
            d.append(val)
    d.sort(key=lambda a: a.name)

    return d

def test_eigenvalue_solvers(data):
    from sfepy.base.base import IndexedStruct

    eig_confs = _list_eigenvalue_solvers(data.conf.solvers)

    all_n_eigs = [5, 0]

    ok = True
    tt = []
    for ii, n_eigs in enumerate(all_n_eigs):
        for eig_conf in eig_confs:
            tst.report(eig_conf.name)

            try:
                eig_solver = Solver.any_from_conf(eig_conf)

            except (ValueError, ImportError):
                if eig_conf.kind in can_fail:
                    continue

                else:
                    raise

            status = IndexedStruct()
            try:
                eigs, vecs = eig_solver(data.mtx, n_eigs=n_eigs,
                                        eigenvectors=True, status=status)

            except KeyboardInterrupt:
                raise

            except:
                if eig_conf.kind in can_fail:
                    continue

                else:
                    raise
            tt.append([' '.join((eig_conf.name, eig_conf.kind)),
                       status.time, n_eigs])

            tst.report(eigs)

            _ok = nm.allclose(eigs.real, eigs_expected[ii],
                              rtol=0.0, atol=1e-8)
            tt[-1].append(_ok)

            ok = ok and (_ok or (eig_conf.kind in can_fail)
                         or (eig_conf.name in can_miss))

    tt.sort(key=lambda x: x[1])
    tst.report('solution times:')
    for row in tt:
        tst.report('%.2f [s] : %s (%d) (ok: %s)'
                   % (row[1], row[0], row[2], row[3]))

    assert ok
