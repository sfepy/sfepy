import os.path as op

import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy import data_dir
from sfepy.fem import Mesh, Domain

def refine(domain, out_dir, level=3):
    for ii in range(level):
        domain = domain.refine()
        filename = op.join(out_dir, 'refine_' + domain.mesh.name + '.mesh')
        domain.mesh.write(filename, io='auto')

    return domain

expected_coors = {
    '2_3' : nm.array([[0.0, 0.0],
                      [1.0, 0.0],
                      [0.0, 1.0],
                      [0.5, 0.0],
                      [0.0, 0.5],
                      [0.5, 0.5]], dtype=nm.float64),
    '2_4' : nm.array([[0.0, 0.0],
                      [1.0, 0.0],
                      [1.0, 1.0],
                      [0.0, 1.0],
                      [0.5, 0.0],
                      [0.0, 0.5],
                      [1.0, 0.5],
                      [0.5, 1.0],
                      [0.5, 0.5]], dtype=nm.float64),
    '3_4' : nm.array([[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0],
                      [0.5, 0.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.0, 0.0, 0.5],
                      [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5],
                      [0.0, 0.5, 0.5]], dtype=nm.float64),
    '3_8' : nm.array([[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [1.0, 1.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0],
                      [1.0, 0.0, 1.0],
                      [1.0, 1.0, 1.0],
                      [0.0, 1.0, 1.0],
                      [0.5, 0.0, 0.0],
                      [0.0, 0.5, 0.0],
                      [0.0, 0.0, 0.5],
                      [1.0, 0.5, 0.0],
                      [1.0, 0.0, 0.5],
                      [0.5, 1.0, 0.0],
                      [1.0, 1.0, 0.5],
                      [0.0, 1.0, 0.5],
                      [0.5, 0.0, 1.0],
                      [0.0, 0.5, 1.0],
                      [1.0, 0.5, 1.0],
                      [0.5, 1.0, 1.0],
                      [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5],
                      [0.0, 0.5, 0.5],
                      [1.0, 0.5, 0.5],
                      [0.5, 1.0, 0.5],
                      [0.5, 0.5, 1.0],
                      [0.5, 0.5, 0.5]], dtype=nm.float64),
}

expected_conn = {
    '2_3' : nm.array([[0, 3, 4],
                      [3, 5, 4],
                      [1, 5, 3],
                      [2, 4, 5]], dtype=nm.int32),
    '2_4' : nm.array([[0, 4, 8, 5],
                      [1, 6, 8, 4],
                      [2, 7, 8, 6],
                      [3, 5, 8, 7]], dtype=nm.int32),
    '3_4' : nm.array([[0, 4, 5, 6],
                      [4, 1, 7, 8],
                      [5, 7, 2, 9],
                      [6, 8, 9, 3],
                      [4, 5, 6, 8],
                      [4, 5, 8, 7],
                      [5, 6, 8, 9],
                      [5, 7, 9, 8]], dtype=nm.int32),
    '3_8' : nm.array([[0,  8, 20,  9, 10, 21, 26, 22],
                      [1, 11, 20,  8, 12, 23, 26, 21],
                      [2, 13, 20, 11, 14, 24, 26, 23],
                      [3,  9, 20, 13, 15, 22, 26, 24],
                      [4, 17, 25, 16, 10, 22, 26, 21],
                      [5, 16, 25, 18, 12, 21, 26, 23],
                      [6, 18, 25, 19, 14, 23, 26, 24],
                      [7, 19, 25, 17, 15, 24, 26, 22]], dtype=nm.int32),
}

def compare_mesh(geo_name, coors, conn):
    _coors = expected_coors[geo_name]
    _conn = expected_conn[geo_name]

    ok = nm.allclose(coors, _coors, rtol=0.0, atol=1e-14)
    ok = ok and (conn == _conn).all()

    return ok

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        mesh = Mesh('mesh_tetra',
                    data_dir + '/meshes/various_formats/small3d.mesh')
        domain = Domain('domain', mesh)

        return Test(conf=conf, options=options, domain=domain)

    def test_facets(self):
        ok = True

        ed = self.domain.ed
        _ok = ed.n_unique == 26
        self.report('unique edges: %s' % _ok)
        ok = ok and _ok

        fa = self.domain.fa
        _ok = fa.n_unique == 30
        self.report('unique faces: %s' % _ok)
        ok = ok and _ok

        return ok

    def test_refine_tetra(self):
        refine(self.domain, self.options.out_dir)

        return True

    def test_refine_hexa(self):
        mesh = Mesh('mesh_hexa',
                    data_dir + '/meshes/various_formats/abaqus_hex.inp')
        domain = Domain('domain', mesh)

        refine(domain, self.options.out_dir)

        return True

    def test_refine_2_3(self):
        mesh = Mesh('2_3', data_dir + '/meshes/elements/2_3_1.mesh')
        domain = refine(Domain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('2_3', domain.mesh.coors, domain.mesh.conns[0])

        return ok

    def test_refine_2_4(self):
        mesh = Mesh('2_4', data_dir + '/meshes/elements/2_4_1.mesh')
        domain = refine(Domain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('2_4', domain.mesh.coors, domain.mesh.conns[0])

        return ok

    def test_refine_3_4(self):
        mesh = Mesh('3_4', data_dir + '/meshes/elements/3_4_1.mesh')
        domain = refine(Domain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('3_4', domain.mesh.coors, domain.mesh.conns[0])

        return ok

    def test_refine_3_8(self):
        mesh = Mesh('3_8', data_dir + '/meshes/elements/3_8_1.mesh')
        domain = refine(Domain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('3_8', domain.mesh.coors, domain.mesh.conns[0])

        return ok
