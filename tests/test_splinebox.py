from sfepy.base.testing import TestCommon
import numpy as nm
from sfepy import data_dir

def tetravolume(cells, vertices):
    vol = 0.0
    c1 = nm.ones((4,4), dtype=nm.float64)
    mul = 1.0 / 6.0
    for ic in cells:
        c1[:,:3] = vertices[ic,:]
        vol += mul * nm.linalg.det(c1)

    return -vol

expected_volumes = (1.22460186e-4, 1.46950423e-4, 1.22460186e-4)

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_spbox(self):
        """
        Check volume change of the mesh which is deformed using
        the SplineBox functions.
        """
        from sfepy.discrete.fem import Mesh
        from sfepy.mesh.splinebox import SplineBox
        mesh = Mesh.from_file(data_dir + '/meshes/3d/cylinder.vtk')
        vol0 = tetravolume(mesh.conns[0], mesh.coors)

        bbox = nm.array(mesh.get_bounding_box()).T
        spbox = SplineBox(bbox, mesh.coors)
        cpoints0 = spbox.get_control_points(init=True)

        for ii in range(4):
            for jj in range(4):
                spbox.change_shape((0, ii, jj), [-0.02, 0, 0])
        coors = spbox.evaluate()
        vol1 = tetravolume(mesh.conns[0], coors)
        mesh.coors = coors

        spbox.set_control_points(cpoints0)
        coors = spbox.evaluate()
        vol2 = tetravolume(mesh.conns[0], coors)

        ok = True
        actual_volumes = (vol0, vol1, vol2)

        for ii in range(3):
            relerr = abs(actual_volumes[ii] - expected_volumes[ii])\
                     / expected_volumes[ii]
            ok = ok and (relerr < 1e-6)

        if not ok:
            self.report('expected volumes:')
            self.report(expected_volumes)
            self.report('actual volumes:')
            self.report(actual_volumes)

        return ok
