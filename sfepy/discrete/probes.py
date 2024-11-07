"""Classes for probing values of Variables, for example, along a line."""
from __future__ import absolute_import
import hashlib

from ast import literal_eval

import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import assert_, basestr, Struct
from sfepy.linalg import make_axis_rotation_matrix, norm_l2_along_axis
import six

def write_results(filename, probe, results):
    """
    Write probing results into a file.

    Parameters
    ----------
    filename : str or file object
        The output file name.
    probe : Probe subclass instance
        The probe used to obtain the results.
    results : dict
        The dictionary of probing results. Keys are data names, values are
        the probed values.
    """
    fd = open(filename, 'w') if isinstance(filename, basestr) else filename

    fd.write('\n'.join(probe.report()) + '\n')
    for key, result in six.iteritems(results):
        pars, vals = result
        vals = vals.reshape(pars.shape[0], -1)
        fd.write('\n# %s %d\n' % (key, vals.shape[-1]))

        aux = nm.concatenate((pars[:, None], vals), axis=1)

        nm.savetxt(fd, aux)

    if isinstance(filename, basestr):
        fd.close()

def read_results(filename, only_names=None):
    """
    Read probing results from a file.

    Parameters
    ----------
    filename : str or file object
        The probe results file name.

    Returns
    -------
    header : Struct instance
        The probe data header.
    results : dict
        The dictionary of probing results. Keys are data names, values are
        the probed values.
    """
    from sfepy.base.ioutils import read_array

    is_only_names = only_names is not None

    fd = open(filename, 'r') if isinstance(filename, basestr) else filename

    header = read_header(fd)
    results = {}
    for name, nc in get_data_name(fd):
        if is_only_names and (name not in only_names): continue

        result = read_array(fd, header.n_point, nc + 1, nm.float64)
        results[name] = result

    return header, results

def read_header(fd):
    """
    Read the probe data header from file descriptor fd.

    Returns
    -------
    header : Struct instance
        The probe data header.
    """
    header = Struct(name='probe_data_header')
    header.probe_class = fd.readline().strip()

    aux = fd.readline().strip().split(':')[1]
    header.n_point = int(aux.strip().split()[0])

    details = []
    while 1:
        line = fd.readline().strip()

        if line == '-----':
            break
        else:
            details.append(line)

    cls = globals()[header.probe_class]
    header.details = cls.parse_report(details)

    return header

def get_data_name(fd):
    """
    Try to read next data name in file fd.

    Returns
    -------
    name : str
        The data name.
    nc : int
        The number of data columns.
    """
    name = None
    while 1:
        try:
            line = fd.readline()
            if (len(line) == 0): break
            if len(line) == 1: continue
        except:
            return

        line = line.strip().split()
        if (len(line) == 3) and (line[0] == '#'):
            name = line[1]
            nc = int(line[2])

            yield name, nc

def parse_vector(line):
    svec = line.split(':')[1].strip('[] ')
    vec = nm.fromstring(svec, dtype=nm.float64, sep=' ')
    return vec

def parse_scalar(line):
    sval = line.split(':')[1].strip('[] ')
    val = literal_eval(sval)
    return val

class Probe(Struct):
    """
    Base class for all point probes. Enforces two points minimum.
    """
    cache = Struct(name='probe_shared_evaluate_cache')
    is_cyclic = False

    def __init__(self, name, share_geometry=True, n_point=None, **kwargs):
        """
        Parameters
        ----------
        name : str
            The probe name, set automatically by the subclasses.
        share_geometry : bool
            Set to True to indicate that all the probes will work on the same
            domain. Certain data are then computed only for the first probe and
            cached.
        n_point : int
           The (fixed) number of probe points, when positive. When non-positive,
           the number of points is adaptively increased starting from -n_point,
           until the neighboring point distance is less than the diameter of the
           elements enclosing the points. When None, it is set to -10.

        For additional parameters see the __init__() docstrings of the
        subclasses.
        """
        Struct.__init__(self, name=name, share_geometry=share_geometry,
                        **kwargs)

        self.set_n_point(n_point)

        self.options = Struct(close_limit=0.1, size_hint=None)
        self.cache = Struct(name='probe_local_evaluate_cache')
        self.acache = Struct(name='probe_actual_evaluate_cache',
                             pars_digest='')

        self.is_refined = False

    def get_evaluate_cache(self):
        """
        Return the evaluate cache for domain-related data given by
        `self.share_geometry`.
        """
        return Probe.cache if self.share_geometry else self.cache

    def get_actual_cache(self, pars, cache, hash_chunk_size=100000):
        """
        Return the actual evaluate cache, which is a combination of the
        (mesh-based) evaluate cache and probe-specific data, like the reference
        element coordinates. The reference element coordinates are reused, if
        the sha1 hash of the probe parameter vector does not change.
        """
        self.acache += cache

        def _gen_array_chunks(arr):
            ii = 0
            while len(arr[ii:]):
                yield arr[ii:ii+hash_chunk_size].tobytes()
                ii += hash_chunk_size

        sha1 = hashlib.sha1()
        for chunk in _gen_array_chunks(pars):
            sha1.update(chunk)

        digest = sha1.hexdigest()
        if digest != self.acache.pars_digest:
            self.acache.pars_digest = digest
            self.acache.ref_coors = None
            self.acache.cells = None
            self.acache.status = None

        return self.acache

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

        self.n_point0 = self.n_point = n_point

    def set_options(self, close_limit=None, size_hint=None):
        """
        Set the probe options.

        Parameters
        ----------
        close_limit : float
            The maximum limit distance of a point from the closest
            element allowed for extrapolation.
        size_hint : float
            Element size hint for the refinement of probe parametrization.
        """
        if close_limit is not None:
            self.options.close_limit = close_limit

        if size_hint is not None:
            self.options.size_hint = size_hint

    def report(self):
        """Report the probe parameters."""
        out = [self.__class__.__name__]

        if self.n_point_required == -1:
            aux = 'adaptive'
        else:
            aux = 'fixed'
        out.append('number of points: %s (%s)' % (self.n_point, aux))

        return out

    def __call__(self, variable, **kwargs):
        """
        Probe the given variable. The actual implementation is in self.probe(),
        so that it can be overridden in subclasses.

        Parameters
        ----------
        variable : Variable instance
            The variable to be sampled along the probe.
        **kwargs : additional arguments
            See :func:`Probe.probe()`.
        """
        return self.probe(variable, **kwargs)

    def probe(self, variable, mode='val', ret_points=False):
        """
        Probe the given variable.

        Parameters
        ----------
        variable : Variable instance
            The variable to be sampled along the probe.
        mode : {'val', 'grad'}, optional
            The evaluation mode: the variable value (default) or the
            variable value gradient.
        ret_points : bool
            If True, return also the probe points.

        Returns
        -------
        pars : array
            The parametrization of the probe points.
        points : array, optional
            If `ret_points` is True, the coordinates of points corresponding to
            `pars`, where the `variable` is evaluated.
        vals : array
            The probed values.
        """
        refine_flag = None

        ev = variable.evaluate_at
        field = variable.field

        cache = field.get_evaluate_cache(cache=self.get_evaluate_cache(),
                                         share_geometry=self.share_geometry)
        self.reset_refinement()

        while True:
            pars, points = self.get_points(refine_flag)
            if not nm.isfinite(points).all():
                raise ValueError('Inf/nan in probe points!')

            acache = self.get_actual_cache(pars, cache)

            vals, ref_coors, cells, status = ev(
                points, mode=mode, strategy='general',
                close_limit=self.options.close_limit, cache=acache,
                ret_ref_coors=True, ret_status=True, ret_cells=True)

            acache.ref_coors = ref_coors
            acache.cells = cells
            acache.status = status

            if self.is_refined:
                break

            else:
                refine_flag = self.refine_points(variable, points, cells)
                if (refine_flag == False).all():
                    break

        self.is_refined = True

        if ret_points:
            return pars, points, vals

        else:
            return pars, vals

    def reset_refinement(self):
        """
        Reset the probe refinement state.
        """
        self.is_refined = False
        self.n_point = self.n_point0

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
            if self.options.size_hint is None:
                ed = variable.get_element_diameters(cells, 0)
                pd = 0.5 * (ed[1:] + ed[:-1])

            else:
                pd = self.options.size_hint

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

    def __init__(self, points, share_geometry=True):
        """
        Parameters
        ----------
        points : array_like
            The coordinates of the points.
        """
        points = nm.array(points, dtype=nm.float64, order='C')
        if points.ndim == 1:
            points.shape = points.shape + (1,)
        n_point = points.shape[0]
        name = 'points %d' % n_point

        Probe.__init__(self, name=name, share_geometry=share_geometry,
                       points=points, n_point=n_point)

        self.n_point_single = n_point

    def report(self):
        """Report the probe parameters."""
        out = Probe.report(self)
        for ii, point in enumerate(self.points):
            out.append('point %d: %s' % (ii, point))
        out.append('-----')
        return out

    @staticmethod
    def parse_report(lines):
        """
        Parse report lines to get the probe parameters.
        """
        out = [parse_vector(line) for line in lines]

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

    def __init__(self, p0, p1, n_point, share_geometry=True):
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

        Probe.__init__(self, name=name, share_geometry=share_geometry,
                       p0=p0, p1=p1, n_point=n_point)

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

    @staticmethod
    def parse_report(lines):
        """
        Parse report lines to get the probe parameters.
        """
        out = [parse_vector(line) for line in lines]
        assert_(len(out) == 2)

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

    def __init__(self, p0, dirvec, p_fun, n_point, both_dirs,
                 share_geometry=True):
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

        Probe.__init__(self, name=name, share_geometry=share_geometry,
                       p0=p0, dirvec=dirvec, p_fun=p_fun,
                       n_point=n_point_true, both_dirs=both_dirs)

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

    @staticmethod
    def parse_report(lines):
        """
        Parse report lines to get the probe parameters.
        """
        out = [parse_vector(lines[0]), parse_vector(lines[1]),
               parse_scalar(lines[2]), lines[3]]

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

    def __init__(self, centre, normal, radius, n_point, share_geometry=True):
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

        Probe.__init__(self, name=name, share_geometry=share_geometry,
                       centre=centre, normal=normal,
                       radius=radius, n_point=n_point)

    def report(self):
        """Report the probe parameters."""
        out = Probe.report(self)
        out.append('centre: %s' % self.centre)
        out.append('normal: %s' % self.normal)
        out.append('radius: %s' % self.radius)
        out.append('-----')
        return out

    @staticmethod
    def parse_report(lines):
        """
        Parse report lines to get the probe parameters.
        """
        out = [parse_vector(lines[0]), parse_vector(lines[1]),
               parse_scalar(lines[2])]

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

        if len(self.centre) == 3:
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
