import os

import numpy as nm

from sfepy.base.testing import TestCommon
from sfepy import data_dir

# n_vertex, n_edge, n_face, n_cell
# d1 -> d2 : num, n_incident
expected = {
    '1_2_2.mesh' : ([3, 2, 0, 0], {
        (0, 0) : (3, 4),
        (0, 1) : (3, 4),
        (1, 0) : (2, 4),
        (1, 1) : (2, 2),
        }),
    '2_3_2.mesh' : ([4, 5, 2, 0], {
        (0, 0) : (4, 10),
        (0, 1) : (4, 10),
        (0, 2) : (4, 6),
        (1, 0) : (5, 10),
        (1, 1) : (5, 16),
        (1, 2) : (5, 6),
        (2, 0) : (2, 6),
        (2, 1) : (2, 6),
        (2, 2) : (2, 2),
        }),
    '2_4_2.mesh' : ([6, 7, 2, 0], {
        (0, 0) : (6, 22),
        (0, 1) : (6, 14),
        (0, 2) : (6, 8),
        (1, 0) : (7, 14),
        (1, 1) : (7, 20),
        (1, 2) : (7, 8),
        (2, 0) : (2, 8),
        (2, 1) : (2, 8),
        (2, 2) : (2, 2),
        }),
    '3_4_2.mesh' : ([5, 9, 7, 2], {
        (0, 0) : (5, 18),
        (0, 1) : (5, 18),
        (0, 2) : (5, 21),
        (0, 3) : (5, 8),
        (1, 0) : (9, 18),
        (1, 1) : (9, 48),
        (1, 2) : (9, 21),
        (1, 3) : (9, 12),
        (2, 0) : (7, 21),
        (2, 1) : (7, 21),
        (2, 2) : (7, 42),
        (2, 3) : (7, 8),
        (3, 0) : (2, 8),
        (3, 1) : (2, 12),
        (3, 2) : (2, 8),
        (3, 3) : (2, 2),
        }),
    '3_8_2.mesh' : ([12, 20, 11, 2], {
        (0, 0) : (12, 100),
        (0, 1) : (12, 40),
        (0, 2) : (12, 44),
        (0, 3) : (12, 16),
        (1, 0) : (20, 40),
        (1, 1) : (20, 96),
        (1, 2) : (20, 44),
        (1, 3) : (20, 24),
        (2, 0) : (11, 44),
        (2, 1) : (11, 44),
        (2, 2) : (11, 72),
        (2, 3) : (11, 12),
        (3, 0) : (2, 16),
        (3, 1) : (2, 24),
        (3, 2) : (2, 12),
        (3, 3) : (2, 2),
        }),
    'square_triquad.mesh' : ([470, 1127, 658, 0], {
        (0, 0) : (470, 3054),
        (0, 1) : (470, 2254),
        (0, 2) : (470, 2174),
        (1, 0) : (1127, 2254),
        (1, 1) : (1127, 9174),
        (1, 2) : (1127, 2174),
        (2, 0) : (658, 2174),
        (2, 1) : (658, 2174),
        (2, 2) : (658, 6686),
        }),
}

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        filename_meshes = [data_dir + '/meshes/elements/%s_2.mesh' % geom
                           for geom in ['1_2', '2_3', '2_4', '3_4', '3_8']]
        filename_meshes.append(data_dir
                               + '/meshes/2d/special/square_triquad.mesh')

        test = Test(filename_meshes=filename_meshes,
                    conf=conf, options=options)
        return test

    def test_cmesh_counts(self):
        from sfepy.discrete.fem import Mesh
        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        from sfepy.discrete.fem.extmods.cmesh import CMesh, get_cmem_usage

        gels = create_geometry_elements()

        ok = True

        for filename in self.filename_meshes:
            basename = os.path.basename(filename)
            enum, esizes = expected[basename]

            self.report('mesh: %s' % basename)

            mesh = Mesh.from_file(filename)
            cmesh = CMesh.from_mesh(mesh)
            cmesh.set_local_entities(gels)

            cmesh.setup_entities()

            self.report('dim:', cmesh.dim)
            self.report('n_vertex: %d, n_edge: %d, n_face: %d, n_cell: %d' %
                        tuple(cmesh.num))

            _ok = (enum == cmesh.num).all()
            if not _ok:
                self.report('%s == %s failed!' % (enum, cmesh.num))
            ok = ok and _ok

            dim = cmesh.dim
            for ir in range(dim + 1):
                for ic in range(dim + 1):
                    cmesh.setup_connectivity(ir, ic)
                    mem_usage1 = get_cmem_usage()[0]

                    if (ir == dim) and (ic == 0):
                        continue

                    cmesh.free_connectivity(ir, ic)
                    mem_usage2 = get_cmem_usage()[0]

                    cmesh.setup_connectivity(ir, ic)
                    mem_usage3 = get_cmem_usage()[0]

                    conn = cmesh.get_conn(ir, ic)

                    self.report('(%d, %d) : (%d, %d)'
                                % (ir, ic, conn.num, conn.n_incident))
                    sizes = nm.array([conn.num, conn.n_incident])

                    _ok = (esizes[ir, ic] == sizes).all()
                    if not _ok:
                        self.report('%s == %s failed!' % (esizes, sizes))
                    ok = ok and _ok

                    _ok1 = mem_usage3 == mem_usage1
                    _ok2 = mem_usage3 > mem_usage2
                    if not (_ok1 and _ok2):
                        self.report('unexpected memory usage! (%s)'
                                    % (mem_usage1, mem_usage2, mem_usage3))
                    ok = ok and (_ok1 and _ok2)

        return ok
