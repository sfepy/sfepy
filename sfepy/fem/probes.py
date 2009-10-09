try:
    from scipy.spatial import KDTree
except ImportError:
    KDTree = None

from sfepy.base.base import *
from sfepy.base.la import make_axis_rotation_matrix
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
        """
        Get the probe points.

        Returns
        -------
        pars : array_like
           The independent coordinate of the probe.
        points : array_like
           The probe points, parametrized by pars.
        """
        pars = nm.linspace(0, self.length, self.n_point)
        points = self.p0 + self.dirvec * pars[:,None]
        return pars, points

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
        pars = self.p_fun(nm.arange(self.n_point, dtype=nm.float64))
        points = self.p0 + sign * self.dirvec * pars[:,None]
        return pars, points

    def get_points(self):
        """
        Get the probe points.

        Returns
        -------
        pars : array_like
           The independent coordinate of the probe.
        points : array_like
           The probe points, parametrized by pars.
        """
        pars, points = self.gen_points(1.0)
        if self.both_dirs:
            pars0, points0 = self.gen_points(-1.0)
            pars = nm.concatenate((-pars0[::-1], pars))
            points = nm.concatenate((points0[::-1], points))
        return pars, points

class CircleProbe(Probe):
    """Probe variables along a circle."""

    def __init__(self, centre, normal, radius, n_point,
                 mesh, share_mesh=True, use_tree=0):
        """
        Parameters
        ----------
        centre : array_like
            The coordinates of the circle centre.
        normal : array_like
            The normal vector perpendicular to the circle plane.
        radius : float
            The radius of the circle.
        n_point : int
            The number of the probe points.
        """
        centre = nm.array(centre, dtype=nm.float64)
        normal = nm.array(normal, dtype=nm.float64)
        normal /= nla.norm(normal)

        name = '[%s, %s, %s] / %d' % (centre, normal, radius, n_point)

        Probe.__init__(self, name=name, mesh=mesh, use_tree=use_tree,
                       centre=centre, normal=normal, radius=radius,
                       n_point=n_point)

    def get_points(self):
        """
        Get the probe points.

        Returns
        -------
        pars : array_like
           The independent coordinate of the probe.
        points : array_like
           The probe points, parametrized by pars.
        """
        # Vector of angles.
        pars = nm.linspace(0.0, 2.0*nm.pi, self.n_point + 1)[:-1]

        # Create the points in xy plane, centered at the origin.
        x = self.radius * nm.cos(pars[:,None])
        y = self.radius * nm.sin(pars[:,None])

        if self.mesh.dim == 3:
            z = nm.zeros((self.n_point, 1), dtype=nm.float64)
            points = nm.c_[x, y, z]

            # Rotate to satisfy the normal, shift to the centre.
            n1 = nm.array([0.0, 0.0, 1.0], dtype=nm.float64)
            axis = nm.cross(n1, self.normal)
            angle = nm.arccos(nm.dot(n1, self.normal))

            if nla.norm(axis) < 0.1:
                # n1 == self.normal
                rot_mtx = nm.eye(3, dtype=nm.float64)
            else:
                rot_mtx = make_axis_rotation_matrix(axis, angle)

            points = nm.dot(points, rot_mtx)

        else:
            points = nm.c_[x, y]
    
        points += self.centre

        return pars, points


class IntegralProbe(Struct):
    """Evaluate integral expressions."""
    def __init__(self, name, problem, expressions, labels):
        Struct.__init__(self, name=name, problem=problem,
                        expressions=expressions, labels=labels)

    def __call__(self, ip, state=None, **kwargs):
        return self.problem.evaluate(self.expressions[ip], state, **kwargs)
