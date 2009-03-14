from sfepy.base.base import *
from sfepy.fem.mesh import make_inverse_connectivity

class Probe(Struct):

    def __call__(self, variable):
        return self.probe(variable)

    def probe(self, variable):
        tt = time.clock()
        points = self.get_points()
        vals = variable.interp_to_points(points, self.mesh, self.iconn)
        print time.clock() - tt

        return vals

class LineProbe(Probe):
    """Probe variables along a line."""

    def __init__(self, p0, p1, n_point, mesh):
        p0 = nm.array(p0, dtype=nm.float64)
        p1 = nm.array(p1, dtype=nm.float64)
        name = '[%s, %s] / %d' % (p0, p1, n_point)

        Probe.__init__(self, name=name, p0=p0, p1=p1, n_point=n_point, mesh=mesh)

        iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                          combine_groups=True)
        self.iconn = iconn

        dirvec = self.p1 - self.p0
        self.length = nm.linalg.norm(dirvec)
        self.dirvec = dirvec / self.length

    def get_points(self):
        points = []
        for ii, step in enumerate(nm.linspace(0, self.length, self.n_point)):
            point = nm.array(self.p0 + self.dirvec * step, ndmin=2)
            points.append( point )
        return points

