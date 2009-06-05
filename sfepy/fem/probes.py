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

    def __init__(self, name, mesh, share_mesh=True, use_tree=0, **kwargs):
        Struct.__init__(self, name=name, mesh=mesh, use_tree=use_tree, **kwargs)

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
    
    def __call__(self, variable):
        return self.probe(variable)

    def probe(self, variable):
        tt = time.clock()
        pars, points = self.get_points()
        vals = variable.interp_to_points(points, self.mesh,
                                         ctree=self.ctree, iconn=self.iconn)
        print time.clock() - tt

        return pars, vals

class LineProbe(Probe):
    """Probe variables along a line."""

    def __init__(self, p0, p1, n_point, mesh, share_mesh=True, use_tree=0):
        p0 = nm.array(p0, dtype=nm.float64)
        p1 = nm.array(p1, dtype=nm.float64)
        name = '[%s, %s] / %d' % (p0, p1, n_point)

        Probe.__init__(self, name=name, mesh=mesh, use_tree=use_tree,
                       p0=p0, p1=p1, n_point=n_point)
            
        dirvec = self.p1 - self.p0
        self.length = nm.linalg.norm(dirvec)
        self.dirvec = dirvec / self.length

    def get_points(self):
        pars = []
        points = []
        for ii, step in enumerate(nm.linspace(0, self.length, self.n_point)):
            point = nm.array(self.p0 + self.dirvec * step, ndmin=2)
            pars.append(step)
            points.append(point)
        return nm.array(pars), nm.array(points)

class RayProbe(Probe):
    """Probe variables along a ray. The points are parametrized by a function of
    radial coordinates from a given point in a given direction."""
    
    def __init__(self, p0, dirvec, p_fun, n_point, both_dirs, mesh,
                 share_mesh=True, use_tree=0):
        p0 = nm.array(p0, dtype=nm.float64)
        dirvec = nm.array(dirvec, dtype=nm.float64)
        dirvec /= nla.norm(dirvec)
        name = '%s [%s, %s] / %d' % (p_fun.__name__, p0, dirvec, n_point)

        Probe.__init__(self, name=name, mesh=mesh, use_tree=use_tree,
                       p0=p0, dirvec=dirvec, p_fun=p_fun, n_point=n_point,
                       both_dirs=both_dirs)

    def gen_points(self, sign):
        pars = []
        points = []
        for ii in xrange(self.n_point):
            step = self.p_fun(ii)
            point = nm.array(self.p0 + sign * self.dirvec * step, ndmin=2)
            pars.append(step)
            points.append( point )
        return nm.array(pars), nm.array(points)

    def get_points(self):
        pars, points = self.gen_points(1.0)
        if self.both_dirs:
            pars0, points0 = self.gen_points(-1.0)
            pars = nm.concatenate((-pars0[::-1], pars))
            points = nm.concatenate((points0[::-1], points))
        return pars, points

class IntegralProbe(Struct):
    """Evaluate integral expressions."""
    def __init__(self, name, problem, expressions, labels):
        Struct.__init__(self, name=name, problem=problem,
                        expressions=expressions, labels=labels)

    def __call__(self, ip, state=None, **kwargs):
        return self.problem.evaluate(self.expressions[ip], state, **kwargs)
