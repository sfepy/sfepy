from __future__ import print_function
from __future__ import absolute_import
import os.path as op

import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy import data_dir
from sfepy.discrete.fem import Mesh, FEDomain
from six.moves import range

def refine(domain, out_dir, level=3):
    for ii in range(level):
        domain = domain.refine()
        filename = op.join(out_dir,
                           'refine_' + op.basename(domain.mesh.name) + '.mesh')
        domain.mesh.write(filename, io='auto')

    return domain

expected_coors = {
    '2_3' : nm.array([[0. , 0. ],
                      [1. , 0. ],
                      [0. , 1. ],
                      [0.5, 0. ],
                      [0.5, 0.5],
                      [0. , 0.5]], dtype=nm.float64),
    '2_4' : nm.array([[0. , 0. ],
                      [1. , 0. ],
                      [1. , 1. ],
                      [0. , 1. ],
                      [0.5, 0. ],
                      [1. , 0.5],
                      [0.5, 1. ],
                      [0. , 0.5],
                      [0.5, 0.5]], dtype=nm.float64),
    '3_4' : nm.array([[0. , 0. , 0. ],
                      [1. , 0. , 0. ],
                      [0. , 1. , 0. ],
                      [0. , 0. , 1. ],
                      [0.5, 0. , 0. ],
                      [0.5, 0.5, 0. ],
                      [0. , 0.5, 0. ],
                      [0. , 0. , 0.5],
                      [0.5, 0. , 0.5],
                      [0. , 0.5, 0.5]], dtype=nm.float64),
    '3_8' : nm.array([[0. , 0. , 0. ],
                      [1. , 0. , 0. ],
                      [1. , 1. , 0. ],
                      [0. , 1. , 0. ],
                      [0. , 0. , 1. ],
                      [1. , 0. , 1. ],
                      [1. , 1. , 1. ],
                      [0. , 1. , 1. ],
                      [0.5, 0. , 0. ],
                      [1. , 0.5, 0. ],
                      [0.5, 1. , 0. ],
                      [0. , 0.5, 0. ],
                      [0.5, 0. , 1. ],
                      [1. , 0.5, 1. ],
                      [0.5, 1. , 1. ],
                      [0. , 0.5, 1. ],
                      [0. , 0. , 0.5],
                      [1. , 0. , 0.5],
                      [1. , 1. , 0.5],
                      [0. , 1. , 0.5],
                      [0.5, 0.5, 0. ],
                      [0. , 0.5, 0.5],
                      [0.5, 0. , 0.5],
                      [0.5, 0.5, 1. ],
                      [1. , 0.5, 0.5],
                      [0.5, 1. , 0.5],
                      [0.5, 0.5, 0.5]], dtype=nm.float64),
}

expected_conn = {
    '2_3' : nm.array([[0, 3, 5],
                      [3, 4, 5],
                      [1, 4, 3],
                      [2, 5, 4]], dtype=nm.int32),
    '2_4' : nm.array([[0, 4, 8, 7],
                      [1, 5, 8, 4],
                      [2, 6, 8, 5],
                      [3, 7, 8, 6]], dtype=nm.int32),
    '3_4' : nm.array([[0, 4, 6, 7],
                      [4, 1, 5, 8],
                      [6, 5, 2, 9],
                      [7, 8, 9, 3],
                      [4, 6, 7, 8],
                      [4, 6, 8, 5],
                      [6, 7, 8, 9],
                      [6, 5, 9, 8]], dtype=nm.int32),
    '3_8' : nm.array([[0,  8, 20, 11, 16, 22, 26, 21],
                      [1,  9, 20,  8, 17, 24, 26, 22],
                      [2, 10, 20,  9, 18, 25, 26, 24],
                      [3, 11, 20, 10, 19, 21, 26, 25],
                      [4, 15, 23, 12, 16, 21, 26, 22],
                      [5, 12, 23, 13, 17, 22, 26, 24],
                      [6, 13, 23, 14, 18, 24, 26, 25],
                      [7, 14, 23, 15, 19, 25, 26, 21]], dtype=nm.int32),
}

def compare_mesh(geo_name, coors, conn):
    _coors = expected_coors[geo_name]
    _conn = expected_conn[geo_name]

    print(coors.__repr__())
    print(conn.__repr__())

    ok = nm.allclose(coors, _coors, rtol=0.0, atol=1e-14)
    ok = ok and (conn == _conn).all()

    return ok

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        mesh = Mesh.from_file(data_dir + '/meshes/various_formats/small3d.mesh')
        domain = FEDomain('domain', mesh)

        return Test(conf=conf, options=options, domain=domain)

    def test_facets(self):
        ok = True
        cmesh = self.domain.cmesh

        _ok = cmesh.num[1] == 26
        self.report('unique edges: %s' % _ok)
        ok = ok and _ok

        _ok = cmesh.num[2] == 30
        self.report('unique faces: %s' % _ok)
        ok = ok and _ok

        return ok

    def test_refine_tetra(self):
        refine(self.domain, self.options.out_dir)

        return True

    def test_refine_hexa(self):
        filename = data_dir + '/meshes/various_formats/abaqus_hex.inp'
        mesh = Mesh.from_file(filename)
        domain = FEDomain('domain', mesh)

        refine(domain, self.options.out_dir)

        return True

    def test_refine_2_3(self):
        mesh = Mesh.from_file(data_dir + '/meshes/elements/2_3_1.mesh')
        domain = refine(FEDomain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('2_3', domain.mesh.coors, domain.mesh.get_conn('2_3'))

        return ok

    def test_refine_2_4(self):
        mesh = Mesh.from_file(data_dir + '/meshes/elements/2_4_1.mesh')
        domain = refine(FEDomain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('2_4', domain.mesh.coors, domain.mesh.get_conn('2_4'))

        return ok

    def test_refine_3_4(self):
        mesh = Mesh.from_file(data_dir + '/meshes/elements/3_4_1.mesh')
        domain = refine(FEDomain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('3_4', domain.mesh.coors, domain.mesh.get_conn('3_4'))

        return ok

    def test_refine_3_8(self):
        mesh = Mesh.from_file(data_dir + '/meshes/elements/3_8_1.mesh')
        domain = refine(FEDomain('domain', mesh), self.options.out_dir, 1)

        ok = compare_mesh('3_8', domain.mesh.coors, domain.mesh.get_conn('3_8'))

        return ok
