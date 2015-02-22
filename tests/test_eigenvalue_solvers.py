import time

import numpy as nm

from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.solvers import Solver

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
        'maxiter' : 100,
    }),
    'evp0h' : ('eig.scipy', {
        'method' : 'eigh',
        'maxiter' : 100,
    }),
    'evp1' : ('eig.sgscipy', {}),
    'evp2' : ('eig.scipy_lobpcg', {
        'i_max' : 100,
        'largest' : False,
    }),
    'evp3' : ('eig.pysparse', {
        'i_max' : 100,
        'eps_a' : 1e-10,
        'strategy' : 0,
    }),
}

eigs_expected = nm.array([0.04904454, 0.12170685, 0.12170685,
                          0.19257998, 0.24082108])

from sfepy.base.testing import TestCommon

class Test(TestCommon):
    can_fail = ['eig.pysparse']

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
        for key, val in confs.iteritems():
            if val.kind.find('eig.') == 0:
                d.append(val)
        d.sort(key=lambda a: a.name)

        return d

    def test_eigenvalue_solvers(self):
        eig_confs = self._list_eigenvalue_solvers(self.conf.solvers)

        n_eigs = 5

        ok = True
        tt = []
        for eig_conf in eig_confs:
            self.report(eig_conf.name)

            eig_solver = Solver.any_from_conf(eig_conf)

            t0 = time.clock()
            eigs, vecs = eig_solver(self.mtx, n_eigs=n_eigs, eigenvectors=True)
            tt.append((' '.join((eig_conf.name, eig_conf.kind)),
                       time.clock() - t0))

            self.report(eigs)

            _ok = nm.allclose(eigs.real, eigs_expected, rtol=0.0, atol=1e-8)
            ok = ok and (_ok or (eig_conf.kind in self.can_fail))

        tt.sort(key=lambda x: x[1])
        self.report('solution times:')
        for row in tt:
            self.report('%.2f [s]' % row[1], ':', row[0])

        return ok
