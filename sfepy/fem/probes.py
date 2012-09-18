"""Classes for probing values of Variables, for example, along a line."""
import time

import numpy as nm
import numpy.linalg as nla
try:
    from scipy.spatial import cKDTree as KDTree
except ImportError:
    from scipy.spatial import KDTree

from sfepy.base.base import output, Struct
from sfepy.linalg import make_axis_rotation_matrix, norm_l2_along_axis
from sfepy.fem.mesh import make_inverse_connectivity

class Probe(Struct):
    """
    Base class for all point probes. Enforces two points minimum.
    """
    cache = Struct(name='probe_shared_cache',
                   offsets=None, iconn=None, kdtree=None)
    is_cyclic = False

    def __init__(self, name, mesh, share_mesh=True, n_point=None, **kwargs):
        """
        Parameters
        ----------
        name : str
            The probe name, set automatically by the subclasses.
        mesh : Mesh instance
            The FE mesh where the variables to be probed are defined.
        share_mesh : bool
            Set to True to indicate that all the probes will work on the same
            mesh. Certain data are then computed only for the first probe and
            cached.
        n_point : int
           The (fixed) number of probe points, when positive. When non-positive,
           the number of points is adaptively increased starting from -n_point,
           until the neighboring point distance is less than the diameter of the
           elements enclosing the points. When None, it is set to -10.

        For additional parameters see the __init__() docstrings of the
        subclasses.

        Notes
        -----
        If the mesh contains vertices that are not contained in any element, we
        shift coordinates of such vertices so that they never match in the
        nearest node search.
        """
        Struct.__init__(self, name=name, mesh=mesh, **kwargs)

        self.set_n_point(n_point)

        self.is_refined = False

        tt = time.clock()
        if share_mesh:
            if Probe.cache.iconn is None:
                offsets, iconn = make_inverse_connectivity(mesh.conns,
                                                           mesh.n_nod,
                                                           ret_offsets=True)
                Probe.cache.iconn = iconn
                Probe.cache.offsets = offsets
            self.cache = Probe.cache

        else:
            offsets, iconn = make_inverse_connectivity(mesh.conns,
                                                       mesh.n_nod,
                                                       ret_offsets=True)
            self.cache = Struct(name='probe_cache',
                                offsets=offsets, iconn=iconn, kdtree=None)
        output('iconn: %f s' % (time.clock()-tt))

        i_bad = nm.where(nm.diff(self.cache.offsets) == 0)[0]
        if len(i_bad):
            bbox = mesh.get_bounding_box()
            mesh.coors[i_bad] = bbox[1] + bbox[1] - bbox[0]
            output('warning: some vertices are not in any element!')
            output('warning: vertex-based results will be wrong!')

        tt = time.clock()
        if share_mesh:
            if Probe.cache.kdtree is None:
                self.cache.kdtree = KDTree(mesh.coors)

        else:
            self.cache.kdtree = KDTree(mesh.coors)

        output('kdtree: %f s' % (time.clock()-tt))

    def set_n_point(self, n_point):
        """
        Set the number of probe points.

        Parameters
        ----------
        n_point : int
           The (fixed) number of probe points, when positive. When non-positive,
           the number of points is adaptively increased starting from -n_point,
           until the neighboring point distance is less than the diameter of the
           elements enclosing the points. When None, it is set to -10.
        """
        if n_point is None:
            n_point = -10

        if n_point <= 0:
            n_point = max(-n_point, 2)
            self.n_point_required = -1

        else:
            n_point = max(n_point, 2)
            self.n_point_required = n_point

        self.n_point = n_point

    def report(self):
        """Report the probe parameters."""
        out = [self.__class__.__name__]

        if self.n_point_required == -1:
            aux = 'adaptive'
        else:
            aux = 'fixed'
        out.append('number of points: %s (%s)' % (self.n_point, aux))

        return out

    def __call__(self, variable):
        """
        Probe the given variable. The actual implementation is in self.probe(),
        so that it can be overridden in subclasses.

        Parameters
        ----------
        variable : Variable instance
            The variable to be sampled along the probe.
        """
        return self.probe(variable)

    def probe(self, variable):
        """
        Probe the given variable.

        Parameters
        ----------
        variable : Variable instance
            The variable to be sampled along the probe.
        """
        refine_flag = None
        while True:
            pars, points = self.get_points(refine_flag)

            vals, cells, status = variable.evaluate_at(points,
                                                       strategy='kdtree',
                                                       cache=self.cache,
                                                       ret_status=True)
            ii = nm.where(status > 0)[0]
            vals[ii] = nm.nan

            if self.is_refined:
                break

            else:
                refine_flag = self.refine_points(variable, points, cells)
                if (refine_flag == False).all():
                    break

        self.is_refined = True

        return pars, vals

    def refine_points(self, variable, points, cells):
        """
        Mark intervals between points for a refinement, based on element
        sizes at those points. Assumes the points to be ordered.

        Returns
        -------
        refine_flag : bool array
            True at places corresponding to intervals between subsequent points
            that need to be refined.
        """
        if self.n_point_required == self.n_point:
            refine_flag = nm.array([False])

        else:
            ed = variable.get_element_diameters(cells, 0)
            pd = 0.5 * (ed[1:] + ed[:-1])

            dist = norm_l2_along_axis(points[1:] - points[:-1])

            refine_flag = dist > pd
            if self.is_cyclic:
                pd1 = 0.5 * (ed[0] + ed[-1])
                dist1 = nla.norm(points[0] - points[-1])

                refine_flag = nm.r_[refine_flag, dist1 > pd1]

        return refine_flag

    @staticmethod
    def refine_pars(pars, refine_flag, cyclic_val=None):
        """
        Refine the probe parametrization based on the refine_flag.
        """
        ii = nm.where(refine_flag)[0]
        ip = ii + 1

        if cyclic_val is not None:
            cpars = nm.r_[pars, cyclic_val]
            pp = 0.5 * (cpars[ip] + cpars[ii])

        else:
            pp = 0.5 * (pars[ip] + pars[ii])

        pars = nm.insert(pars, ip, pp)

        return pars

class PointsProbe(Probe):
    """
    Probe variables in given points.
    """
    
    def __init__(self, points, mesh, share_mesh=True):
        """
        Parameters
        ----------
        points : array_like
            The coordinates of the points.
        """
        points = nm.array(points, dtype=nm.float64)
        if points.ndim == 1:
            points.shape = points.shape + (1,)
        n_point = points.shape[0]
        name = 'points %d' % n_point

        Probe.__init__(self, name=name, mesh=mesh,
                       points=points, n_point=n_point)

        self.n_point_single = n_point

    def report(self):
        """Report the probe parameters."""
        out = Probe.report(self)
        for ii, point in enumerate(self.points):
            out.append('point %d: %s' % (ii, point))
        out.append('-----')
        return out

    def refine_points(self, variable, points, cache):
        """No refinement for this probe."""
        refine_flag = nm.array([False])
        return refine_flag
    
    def get_points(self, refine_flag=None):
        """
        Get the probe points.

        Returns
        -------
        pars : array_like
           The independent coordinate of the probe.
        points : array_like
           The probe points, parametrized by pars.
        """
        pars = nm.arange(self.n_point, dtype=nm.float64)
        return pars, self.points

class LineProbe(Probe):
    """
    Probe variables along a line.

    If n_point is positive, that number of evenly spaced points is used. If
    n_point is None or non-positive, an adaptive refinement based on element
    diameters is used and the number of points and their spacing are determined
    automatically. If it is negative, -n_point is used as an initial guess.
    """

    def __init__(self, p0, p1, n_point, mesh, share_mesh=True):
        """
        Parameters
        ----------
        p0 : array_like
            The coordinates of the start point.
        p1 : array_like
            The coordinates of the end point.
        """
        p0 = nm.array(p0, dtype=nm.float64)
        p1 = nm.array(p1, dtype=nm.float64)
        name = 'line [%s, %s]' % (p0, p1)

        Probe.__init__(self, name=name, mesh=mesh, p0=p0, p1=p1, n_point=n_point)
            
        dirvec = self.p1 - self.p0
        self.length = nm.linalg.norm(dirvec)
        self.dirvec = dirvec / self.length

    def report(self):
        """Report the probe parameters."""
        out = Probe.report(self)
        out.append('point 0: %s' % self.p0)
        out.append('point 1: %s' % self.p1)
        out.append('-----')
        return out

    def get_points(self, refine_flag=None):
        """
        Get the probe points.

        Returns
        -------
        pars : array_like
           The independent coordinate of the probe.
        points : array_like
           The probe points, parametrized by pars.
        """
        if self.is_refined:
            return self.pars, self.points

        if refine_flag is None:
            pars = nm.linspace(0, self.length, self.n_point)
            
        else:
            pars = Probe.refine_pars(self.pars, refine_flag)

            self.n_point = pars.shape[0]

        self.pars = pars

        self.points = self.p0 + self.dirvec * pars[:,None]

        return pars, self.points

class RayProbe(Probe):
    """
    Probe variables along a ray. The points are parametrized by a function of
    radial coordinates from a given point in a given direction.
    """
    
    def __init__(self, p0, dirvec, p_fun, n_point, both_dirs, mesh,
                 share_mesh=True):
        """
        Parameters
        ----------
        p0 : array_like
            The coordinates of the start point.
        dirvec : array_like
            The probe direction vector.
        p_fun : function
            The function returning the probe parametrization along the dirvec
            direction.
        both_dirs : bool
            If True, the probe works, starting at p0, symmetrically in both
            dirvec and -dirvec directions.
        """
        p0 = nm.array(p0, dtype=nm.float64)
        dirvec = nm.array(dirvec, dtype=nm.float64)
        dirvec /= nla.norm(dirvec)
        name = 'ray %s [%s, %s]' % (p_fun.__name__, p0, dirvec)

        if both_dirs:
            n_point_true  = 2 * n_point
        else:
            n_point_true = n_point

        Probe.__init__(self, name=name, mesh=mesh,
                       p0=p0, dirvec=dirvec, p_fun=p_fun, n_point=n_point_true,
                       both_dirs=both_dirs)

        self.n_point_single = n_point

    def report(self):
        """Report the probe parameters."""
        out = Probe.report(self)
        out.append('point 0: %s' % self.p0)
        out.append('direction vector: %s' % self.dirvec)
        out.append('both directions: %s' % self.both_dirs)
        out.append('distribution function: %s' % self.p_fun.__name__)
        out.append('-----')
        return out

    def refine_points(self, variable, points, cache):
        """No refinement for this probe."""
        refine_flag = nm.array([False])
        return refine_flag
    
    def gen_points(self, sign):
        """Generate the probe points and their parametrization."""
        pars = self.p_fun(nm.arange(self.n_point_single, dtype=nm.float64))
        points = self.p0 + sign * self.dirvec * pars[:,None]
        return pars, points

    def get_points(self, refine_flag=None):
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
    """
    Probe variables along a circle.

    If n_point is positive, that number of evenly spaced points is used. If
    n_point is None or non-positive, an adaptive refinement based on element
    diameters is used and the number of points and their spacing are determined
    automatically. If it is negative, -n_point is used as an initial guess.
    """
    is_cyclic = True

    def __init__(self, centre, normal, radius, n_point,
                 mesh, share_mesh=True):
        """
        Parameters
        ----------
        centre : array_like
            The coordinates of the circle centre.
        normal : array_like
            The normal vector perpendicular to the circle plane.
        radius : float
            The radius of the circle.
        """
        centre = nm.array(centre, dtype=nm.float64)
        normal = nm.array(normal, dtype=nm.float64)
        normal /= nla.norm(normal)

        name = 'circle [%s, %s, %s]' % (centre, normal, radius)

        Probe.__init__(self, name=name, mesh=mesh,
                       centre=centre, normal=normal, radius=radius,
                       n_point=n_point)

    def report(self):
        """Report the probe parameters."""
        out = Probe.report(self)
        out.append('centre: %s' % self.centre)
        out.append('normal: %s' % self.normal)
        out.append('radius: %s' % self.radius)
        out.append('-----')
        return out

    def get_points(self, refine_flag=None):
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
        if self.is_refined:
            return self.pars, self.points

        if refine_flag is None:
            pars = nm.linspace(0.0, 2.0*nm.pi, self.n_point + 1)[:-1]

        else:
            pars = Probe.refine_pars(self.pars, refine_flag,
                                     cyclic_val=2.0 * nm.pi)

            self.n_point = pars.shape[0]

        self.pars = pars

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

        self.points = points

        return pars, points


class IntegralProbe(Struct):
    """Evaluate integral expressions."""
    def __init__(self, name, problem, expressions, labels):
        Struct.__init__(self, name=name, problem=problem,
                        expressions=expressions, labels=labels)

    def __call__(self, ip, state=None, **kwargs):
        return self.problem.evaluate(self.expressions[ip], state, **kwargs)
