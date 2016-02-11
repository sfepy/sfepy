from sfepy.base.testing import TestCommon
import numpy as nm
from sfepy import data_dir
from sfepy.mesh.splinebox import SplineBox, SplineRegion2D
from sfepy.mesh.bspline import BSpline
from sfepy.discrete.fem import Mesh

def tetravolume(cells, vertices):
    vol = 0.0
    c1 = nm.ones((4,4), dtype=nm.float64)
    mul = 1.0 / 6.0
    for ic in cells:
        c1[:,:3] = vertices[ic,:]
        vol += mul * nm.linalg.det(c1)

    return -vol

tolerance = 1e-6

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_spbox_3d(self):
        """
        Check volume change of the mesh which is deformed using
        the SplineBox functions.
        """
        from sfepy.discrete.fem import Mesh
        from sfepy.mesh.splinebox import SplineBox
        mesh = Mesh.from_file(data_dir + '/meshes/3d/cylinder.vtk')
        conn = mesh.get_conn('3_4')
        vol0 = tetravolume(conn, mesh.coors)

        bbox = nm.array(mesh.get_bounding_box()).T
        spbox = SplineBox(bbox, mesh.coors)
        cpoints0 = spbox.get_control_points(init=True)

        for ii in range(4):
            for jj in range(4):
                spbox.move_control_point((0, ii, jj), [-0.02, 0, 0])
        coors = spbox.evaluate()
        vol1 = tetravolume(conn, coors)
        mesh.coors[:] = coors

        spbox.set_control_points(cpoints0)
        coors = spbox.evaluate()
        vol2 = tetravolume(conn, coors)

        ok = True
        actual_volumes = (vol0, vol1, vol2)
        expected_volumes = (1.22460186e-4, 1.46950423e-4, 1.22460186e-4)

        for ii in range(3):
            relerr = abs(actual_volumes[ii] - expected_volumes[ii])\
                     / expected_volumes[ii]
            ok = ok and (relerr < tolerance)

        if not ok:
            self.report('expected volumes:')
            self.report(expected_volumes)
            self.report('actual volumes:')
            self.report(actual_volumes)

        return ok

    def test_spbox_2d(self):
        """
        Check position of a given vertex in the deformed mesh.
        """
        mesh = Mesh.from_file('meshes/2d/square_tri1.mesh')
        spb = SplineBox([[-1, 1], [-1, 0.6]], mesh.coors, nsg=[2,1])
        spb.move_control_point(1, [0.1, -0.2])
        spb.move_control_point(2, [0.2, -0.3])
        spb.move_control_point(3, [0.0, -0.1])

        pt0 = mesh.coors[175,:].copy()
        mesh.cmesh.coors[:] = spb.evaluate()
        pt1 = mesh.coors[175,:]

        expected_distance = 0.165892726387
        actual_distance = nm.linalg.norm(pt0 - pt1)
        ok = nm.fabs(actual_distance - expected_distance)\
            / expected_distance < tolerance

        if not ok:
            self.report('expected distance:')
            self.report(expected_distance)
            self.report('actual distance:')
            self.report(actual_distance)

        return ok

    def test_spregion2d(self):
        """
        Check position of a given vertex in the deformed mesh.
        """
        line_l = nm.array([[-1, 1], [-1, .5], [-1, 0], [-1, -.5]])
        line_r = nm.array([[0, -.2], [.1, .2], [.3, .6], [.4, 1]])
        sp_l = BSpline(3, is_cyclic=False)
        sp_l.approximate(line_l, ncp=4)
        kn_lr = sp_l.get_knot_vector()
        sp_r = BSpline(3, is_cyclic=False)
        sp_r.approximate(line_r, knots=kn_lr)

        line_b = nm.array([[-1, -.5], [-.8, -.6], [-.5, -.4], [-.2, -.2],
                           [0, -.2]])
        line_t = nm.array([[.4, 1], [0, 1], [-.2, 1], [-.6, 1], [-1, 1]])
        sp_b = BSpline(3, is_cyclic=False)
        sp_b.approximate(line_b, ncp=5)
        kn_bt = sp_b.get_knot_vector()
        sp_t = BSpline(3, is_cyclic=False)
        sp_t.approximate(line_t, knots=kn_bt)

        mesh = Mesh.from_file('meshes/2d/square_tri1.mesh')
        spb = SplineRegion2D([sp_b, sp_r, sp_t, sp_l], mesh.coors)
        spb.move_control_point(5, [-.2, .1])
        spb.move_control_point(10, [-.3, .2])
        spb.move_control_point(15, [-.1, .2])

        pt0 = mesh.coors[145,:].copy()
        mesh.cmesh.coors[:] = spb.evaluate()
        pt1 = mesh.coors[145,:]

        expected_distance = 0.0908306614584
        actual_distance = nm.linalg.norm(pt0 - pt1)
        ok = nm.fabs(actual_distance - expected_distance)\
            / expected_distance < tolerance

        if not ok:
            self.report('expected distance:')
            self.report(expected_distance)
            self.report('actual distance:')
            self.report(actual_distance)

        return ok
