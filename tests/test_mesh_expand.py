from __future__ import absolute_import
from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_mesh_expand(self):
        import numpy as nm
        from sfepy.mesh.mesh_tools import expand2d
        from sfepy.discrete.fem.mesh import Mesh
        from sfepy import data_dir

        mesh2d = Mesh.from_file(data_dir + '/meshes/2d/square_quad.mesh')
        mesh3d_ref = Mesh.from_file(data_dir +\
                                    '/meshes/3d/cube_medium_hexa.mesh')
        mesh3d_gen = expand2d(mesh2d, 0.1, 10)
        mesh3d_gen.coors[:,2] += -0.5

        d0 = nm.sort(nm.linalg.norm(mesh3d_ref.coors, axis=1))
        d1 = nm.sort(nm.linalg.norm(mesh3d_gen.coors, axis=1))
        if nm.linalg.norm(d0 - d1) < 1e-6:
            return True

        else:
            return False

        return True
