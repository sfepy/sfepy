from __future__ import absolute_import

import numpy as nm

from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.solvers import Solver
import six

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
}

eigs_expected = [nm.array([0.04904454, 0.12170685, 0.12170685,
                           0.19257998, 0.24082108]),
                 []]

from sfepy.base.testing import TestCommon

class Test(TestCommon):
    can_fail = ['eig.slepc', 'eig.matlab']
    can_miss = ['evp0'] # Depending on scipy version, evp0 can miss an
                        # eigenvalue.

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete import Problem

        pb = Problem.from_conf(conf, init_solvers=False)
        pb.time_update()
        mtx = pb.equations.evaluate(mode='weak', dw_mode='matrix',
                                    asm_obj=pb.mtx_a)

        test = Test(mtx=mtx, conf=conf, options=options)
        return test

    def _list_eigenvalue_solvers(self, confs):
        d = []
        for key, val in six.iteritems(confs):
            if val.kind.find('eig.') == 0:
                d.append(val)
        d.sort(key=lambda a: a.name)

        return d

    def test_eigenvalue_solvers(self):
        from sfepy.base.base import IndexedStruct

        eig_confs = self._list_eigenvalue_solvers(self.conf.solvers)

        all_n_eigs = [5, 0]

        ok = True
        tt = []
        for ii, n_eigs in enumerate(all_n_eigs):
            for eig_conf in eig_confs:
                self.report(eig_conf.name)

                try:
                    eig_solver = Solver.any_from_conf(eig_conf)

                except (ValueError, ImportError):
                    if eig_conf.kind in self.can_fail:
                        continue

                    else:
                        raise

                status = IndexedStruct()
                eigs, vecs = eig_solver(self.mtx, n_eigs=n_eigs,
                                        eigenvectors=True, status=status)
                tt.append([' '.join((eig_conf.name, eig_conf.kind)),
                           status.time, n_eigs])

                self.report(eigs)

                _ok = nm.allclose(eigs.real, eigs_expected[ii],
                                  rtol=0.0, atol=1e-8)
                tt[-1].append(_ok)

                ok = ok and (_ok or (eig_conf.kind in self.can_fail)
                             or (eig_conf.name in self.can_miss))

        tt.sort(key=lambda x: x[1])
        self.report('solution times:')
        for row in tt:
            self.report('%.2f [s] : %s (%d) (ok: %s)'
                        % (row[1], row[0], row[2], row[3]))

        return ok
