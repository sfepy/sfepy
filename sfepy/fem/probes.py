try:
    from scipy.spatial import KDTree
except ImportError:
    KDTree = None

from sfepy.base.base import *
from sfepy.fem.mesh import make_inverse_connectivity, TreeItem

class Probe(Struct):
    iconn = None
    ctree = None
    cttype = None
    
    def __call__(self, variable):
        return self.probe(variable)

    def probe(self, variable):
        tt = time.clock()
        points = self.get_points()
        vals = variable.interp_to_points(points, self.mesh,
                                         ctree=self.ctree, iconn=self.iconn)
        print time.clock() - tt

        return vals

class LineProbe(Probe):
    """Probe variables along a line."""

    def __init__(self, p0, p1, n_point, mesh, share_mesh=True, use_tree=0):
        p0 = nm.array(p0, dtype=nm.float64)
        p1 = nm.array(p1, dtype=nm.float64)
        name = '[%s, %s] / %d' % (p0, p1, n_point)

        Probe.__init__(self, name=name, p0=p0, p1=p1, n_point=n_point,
                       mesh=mesh, use_tree=use_tree)

        tt = time.clock()
        if share_mesh and Probe.iconn:
            iconn = Probe.iconn
        else:
            iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                              combine_groups=True)
            Probe.iconn = iconn
        self.iconn = iconn
        output('iconn: %f s' % (time.clock()-tt))

        if use_tree:
            
            tt = time.clock()
            if share_mesh and Probe.ctree:
                self.ctree = Probe.ctree
                self.cttype = Probe.cttype
            else:
                if (use_tree == 1) and (KDTree is not None):
                    tt = time.clock()                    
                    ctree = KDTree(mesh.coors)
                    cttype = 'scipy.kdtree'
                else:
                    tt = time.clock()                    
                    ctree = TreeItem.build_tree(mesh.coors, use_tree, 2)
                    cttype = 'builtin (slow)'
                Probe.ctree = self.ctree = ctree
                Probe.cttype = self.cttype = cttype
            output('ctree (%s): %f s' % (self.cttype, time.clock()-tt))
                
        else:
            self.ctree = None
            
        dirvec = self.p1 - self.p0
        self.length = nm.linalg.norm(dirvec)
        self.dirvec = dirvec / self.length

    def get_points(self):
        points = []
        for ii, step in enumerate(nm.linspace(0, self.length, self.n_point)):
            point = nm.array(self.p0 + self.dirvec * step, ndmin=2)
            points.append( point )
        return points

