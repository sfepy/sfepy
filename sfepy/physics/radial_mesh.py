import numpy as np
import math
from scipy.integrate import simps, cumtrapz
from sfepy.linalg import norm_l2_along_axis
import scipy.interpolate as si

class RadialVector(object):

    """ Read array of radial vectors from file """
    @staticmethod
    def from_file(file):
        array = np.genfromtxt(file)
        mesh = ExplicitRadialMesh(array[:, 0])
        return [RadialVector(mesh, array[:, r])
                for r in xrange(1, array.shape[1])]

    def __repr__(self):
        return self.pretty(10)

    def pretty(self, values):
        """ Pretty print <values> values of vector 
        @param int number of values printed 
        """
        size = self.values.size
        if values > size:
            n = np.arange(size, dtype=np.int)
        else:
            n = np.array(np.linspace(0, self.mesh.size - 1, values),
                         dtype=np.int)
        out = 'RadialVector: '
        for x in n[0:-1:1]:
            out += '%f: %f,  ' % (self.mesh.get_r(x), self.values[x])
        x = n[-1]
        out += '%f: %f' % (self.mesh.get_r(x), self.values[x])
        return out

    def __init__(self, mesh, values=None):
        if not isinstance(mesh, RadialMesh):
            mesh = ExplicitRadialMesh(mesh)
        self.mesh = mesh
        self.values = (values if values
                       is not None else np.zeros(self.mesh.size))
        self.extrapolated = None
        self.precision = None
        self.result_precision = None
        self.running = None
        self.derivate_cache = {}
        self.interpolator = {}
        
    def running_mean(self):
        """ Smooths vector """
        if self.running is None:
            weights = np.array((1.0, 3, 4, 3, 1))
            weights = weights / weights.sum()
            wsize = int(weights.size - 1)
            data = np.hstack([np.ones(wsize) * self.values[0], self.values,
                              np.ones(wsize / 2) * self.values[-1]])
            self.running = np.convolve(data, weights)[wsize:-(wsize / 2)
                                                      - wsize]
        return self.running

    def integrate(self, factor = 'spherical'):
        """ Compute integral of vector with given integral factor. 
            See RadialMesh.integrate for possible factors.
        """
        return self.mesh.integrate(self.values, factor)

    def linear_integrate(self):
        """ Compute (not-spherical) integral of vector """
        return self.mesh.linear_integrate(self.values)

    def linear_integral(self, from_zero=False):
        """ Compute running integral of vector. See mesh.linear_integral """
        return self.mesh.linear_integral(self.values, from_zero)

    def get_extrapolated(self, precision=0.0001, grade=10, attempts=10):
        """ Try to smooth and extrapolate curve trough vector points with given precision,
        if can't be done, it raise the precision
        @param precision int max error in points
        @param grade int grade of spline
        @param attempts max number of tries to gues the precision 
        """
        
        if precision is None:
            return si.InterpolatedUnivariateSpline(self.mesh.get_coors(),
                                                   self.values, k=5)

        if ((self.extrapolated is None)
            or (self.precision == self.result_precision)
            and (precision < self.precision)):

            self.precision = precision
            data = self.running_mean()
            for attempt in xrange(attempts):
                self.extrapolated = si.UnivariateSpline(self.mesh.get_coors(),
                                                        data, k=5, s=precision)
                der = self.extrapolated(self.mesh.get_coors(), 1)
                sig = np.sign(der)
                if np.abs(sig[1:] - sig[:-1]).sum() <= grade:
                    break
                precision = precision * 2
        self.result_precision = precision
        a = self.extrapolated
        return a

    def extrapolated_values(self, at=None, precision=0.0001, grade=10,
                            attempts=10):
        """ Smoth the vector by interposing spline curve and return smoothed values """
        if at is None:
            at = self.mesh
        elif not isinstance(at, RadialMesh):
            at = ExplicitRadialMesh(at)
        val = self.get_extrapolated(precision, grade, attempts)(at.get_coors())
        return RadialVector(at, val)

    def extrapolated_derivatives(self, at=None, precision=0.0001, attempts=10):
         """ Smoth the vector by interposing spline curve and return derivatives """
        if at is None:
            at = self.mesh
        elif not isinstance(at, RadialMesh):
            at = ExplicitRadialMesh(at)
        val = self.get_extrapolated(precision=0.0001, grade=10,
                                    attempts=10)(at.get_coors(), 1)
        return RadialVector(at, val)

    def derivation(self, at = None, factor = 'spherical', force=False):
        """ Return radial vector of derivatives and if 'at' is not none return derivatives 
            in given points 
            Can derivate with respect to given integral factor, see mesh.integrate 
        """
        radial = bool(radial) 
        if not self.derivate_cache.has_key(factor) or force:
           self.derivate_cache[factor] = RadialVector(self.mesh, self.mesh.derivation(self.values, factor))
        if at is not None:
           return self.derivate_cache[factor].interpolate(at)
        return self.derivate_cache[factor]
        
    def linear_derivation(self, at=None, force = False):
        """ Return radial vector of derivatives and if 'at' is not none return derivatives 
            in given points.
            Linear derivation (no integration factor).
        """
        return self.derivation(at = at, radial = False, force = force)

    def slice(self, x, y):
        """ Return slice of vector, given by two indexes (start, stop+1) or
            float bounds <a, b) """
        if isinstance(x, float):
            x = self.get_index(x)
        if isinstance(y, float):
            y = self.get_index(y)
        return RadialVector(self.mesh.slice(x, y), self.values[x:y])
    
   def interpolate(self, x, kind = None, centre = None):
        """ Return values interpolated in given coordinates x.
        Kind of interpolation can be None for linear interpolation or 
        some of scipy.interpolate.interp1d kinds.
        
        Caches computed interpolator for non-linear cases.
        If centre is not None, coordinates are 3d and the resulting coordinates 
        are the distances from the given centre.
        """
        if kind is None:
           return self.mesh.interpolate(self.values, x, centre, kind) 
        if not self.interpolator.has_key(kind):
           self.interpolator[kind]=self.mesh.interpolator(self.values, kind)
        if(centre is not None):
           if centre.any():
              x = x - centre
           x = norm_l2_along_axis(x, axis=1)
        return self.interpolator[kind](x)

    def interpolate_3d(self, x, kind = None, centre = None):
        """ Return interpolated values at points given by 3d coordinates and centre """
        if centre is None:
          centre = np.zeros((3,))
        return self.interpolate(x, kind, centre)    def output_vector(self, filename=None):
        """ """
        return self.mesh.output_vector(self, filename)

    @staticmethod
    def sparse_merge(vectors):
        mesh = RadialMesh.merge([v.mesh for v in vectors])
        return [mesh.sparse_vector(v) for v in vectors]

    def _get_values_from_object(self, object):
        return (object.values if isinstance(object, RadialVector) else object)

    def __add__(self, vector):
        values = self._get_values_from_object(vector)
        return RadialVector(self.mesh, self.values + values)

    def __sub__(self, vector):
        values = self._get_values_from_object(vector)
        return RadialVector(self.mesh, self.values - values)

    def __mul__(self, vector):
        values = self._get_values_from_object(vector)
        return RadialVector(self.mesh, self.values * values)

    def __div__(self, vector):
        values = self._get_values_from_object(vector)
        return RadialVector(self.mesh, self.values / values)

    def __getitem__(self, r):
        return self.values[r]

    def plot(self):
        """ Plot vector using given shell script """
        return self.mesh.plot(self)

    def get_coors(self):
        """ Return mesh coors """
        return self.mesh.get_coors()

    def __call__(self, name, *args, **kwargs):
        """ Call numpy functions on vector """
        return getattr(self.values, name)(*args, **kwargs)

    def __getattr__(self, name):
        """ Read numpy array attributes """
        return getattr(self.values, name)

    def partition(self, block=None):
        """
        Split again vector joined from more radial vectors.
        The mesh.coors should look like a saw, it's divided to single tooths.  
        """
        meshes = self.mesh.partition()
        out = []
        start = 0
        for m in meshes:
            size = m.size
            out.append(RadialVector(m, self.values[start:start + size]))
            start += size
        if block is not None:
           out=out[block]
        return out
    
    def dot(self, vector, factor='spherical'):
        """ Return dot product with given vector using given integral factor """
        return self.mesh.dot(self.values, vector, factor)
            
    def v_dot(self, a, b, factor='spherical'):
        """ Return dot product of two vectors in self-norm using given integral factor 
            .. math::
               (a, b)_{vector} = \int a * (b * vector)' * factor dx 
        """
        return self.mesh.dot(a, self.values*b, factor)

class RadialMesh(object):
    """
    Radial mesh.
    """
    @staticmethod
    def integral_factor(vector, factor):
        """
        return values multiplied by given integral factor
        It can be 
        'linear' - no integral factor
        'spehrical' - spherical integral factor: 4*math.pi * r**2
        'r2' - physics spherical integral factor (used e.g. with spherical harmonics): r**2
        """
        if isinstance(vector, RadialVector):
           vector = vector.values
        if factor == 'r2':
           return vector * self.coors ** 2
        elif factor == 'linear':
           return vector
        elif factor == 'spherical' or True:
           return vector * self.coors ** 2 * 4 * math.pi
        raise Exception('Unknown integral factor: %s' % norm)        
 

    
    def derivation(self, vector, factor = 'spherical'):
        """
        Return vector from derivations of given vector that was integrated using given integral factor
        .. math::
          f_n = \frac{v(x_n - \epsilon) - v(x_n + \epsilon)}{2 * \epsilon * factor} 
        """
        spaces = self.integral_factor(coors, factor)
        spaces = np.convolve(spaces, [1, -1])[1:-1]        
        diffs = np.convolve(vector, [1, -1])[1:-1]
        diffs /= spaces
        out = np.convolve(diffs, [0.5, 0.5])
        out[0] = diffs[0]
        out[-1] = diffs[-1]
        return RadialVector(self, out)

    def linear_derivation(self, vector, at = None):
        """
        Return vector from derivations of given vector (with no integral factor)
        """
        return self.derivation(vector, radial = False)

    def interpolate_3d(self, potential, coors, centre=None, kind=None):
        """
        Interpolate values in points giving by 3d coordinates with respect to given centre.
        """        
        if centre is None:
           centre = np.zeros((3,))
        return self.interpolate(potential, coors, centre, kind)        

    def integrate(self, vector, factor = 'spherical'):
        """
        .. math::
           \int f(r) r^2 dr
        """
        v = RadialMesh.integral_factor(vector, factor)
        return simps(v, self.coors)

    def linear_integrate(self, vector):
        """
        .. math::
           \int f(r) dr
        """
        return simps(vector, self.coors)

    def integral(self, vector, factor = 'spherical', from_zero=False):
        """
        .. math::
          a_n = \int_{r_0}^{r_n} f(r) * integral_factor dr

        from_zero starts to integrate from zero, instead of starting between
        the first two points
        """        
        v = RadialMesh.integral_factor(vector, factor)
        r = self.get_coors()
        if from_zero:
            v = cumtrapz(vector, r, initial=vector[0] / 2 * max(0.0, r[0]))
            return RadialVector(self, v)
        else:
            v = cumtrapz(vector, r)
            r = r[1:]
            return RadialVector(r, v)
        
    def linear_integral(self, vector, from_zero=False):
        """
        .. math::
          a_n = \int_{r_0}^{r_n} f(r) dr

        from_zero starts to integrate from zero, instead of starting between
        the first two points
        """
        return self.integral(vector, 'linear', from_zero)
        

    def dot(self, vector_a, vector_b, factor='spherical'):
        """
        Dot product with respect to given integral factor (see RadialMesh.integral_factor)
        .. math::
           \int f(r) g(r) factor r^2 dr
        """
        
        if isinstance(vector_a, RadialVector):
           vector_a=vector_a.values
        if isinstance(vector_b, RadialVector):
           vector_b=vector_b.values
        return self.integrate(vector_a * vector_b, factor)

    def norm(self, vector, factor = 'spherical'):
        """
        L2 norm with respect to given integral factor (see RadialMesh.integral_factor)
        
        .. math::
           \sqrt(\int f(r) f(r) factor r^2 dr)
        """
        return np.sqrt(self.dot(vector, vector, factor))

    def output_vector(self, vector, filename=None):
        """
        Output vector or vectors in columns prepended by mesh coordinates to file or stdout if no filename given.
        Output format
        corrs[0] v[0]
        corrs[1] v[1]
        etc... or 
        corrs[0] v0[0] v1[0]] ...
        corrs[1] v0[1] v1[1]] ...        
        """
        if filename is None:
            import sys
            filename = sys.stdout
        if isinstance(vector, RadialVector):
            vector = [vector]
        if isinstance(vector[0], RadialVector):
            vector = [v.values for v in vector]
        if isinstance(vector, list):
          vector = np.vstack([self.coors] + vector)
        else:
          vector = np.vstack([self.coors, vector])
        np.savetxt(filename, vector.T)

    def plot(self, vector, cmd='plot'):
        """ Plot given vectors using given shell script"""
        import tempfile
        import os
        fhandle, fname = tempfile.mkstemp()
        fil = os.fdopen(fhandle, 'w')
        self.output_vector(vector, fil)
        fil.close()
        os.system('%s %s' % (cmd, fname))
        os.remove(fname)

    @staticmethod
    def merge(meshes):
        """ Merge more radial meshes to one """
        merged = np.concatenate(tuple(m.get_coors() for m in meshes))
        return ExplicitRadialMesh(np.unique(merged))

    def intervals(self):
        """ Return distances between coors """
        return np.convolve(self.coors, [1, -1])[1:-1]

class ExplicitRadialMesh(RadialMesh):
    """   Radial mesh given by explicit mesh point coordinates   """
    def __init__(self, coors):
        self.coors = coors
        self.midpointMesh = {}
        self.parentMesh = None

    @property
    def shape(self):
        return self.coors.shape

    @property
    def size(self):
        return self.coors.size

    def get_coors(self):
        return self.coors

    def last_point(self):
        return self.coors[self.coors.size - 1]

    def get_r(self, index):
        return self.coors[index]

    def get_index(self, r):
        pos = self.coors.searchsorted(r)
        return (pos if pos < self.coors.size else self.coors.size - 1)

    def get_mixing(self, r):
        pos = self.get_index(r)
        if pos == self.coors.size - 1 and self.coors[pos] < r:
            out = [(pos, 1.0)]

        elif pos == 0 or self.coors[pos] == r:
            out = [(pos, 1.0)]

        else:
            pos_c = (r - self.coors[pos - 1]) / (self.coors[pos]
                                                 - self.coors[pos - 1])
            out = [(pos - 1, 1.0 - pos_c), (pos, pos_c)]

        return out

    def interpolator(self, potential, kind = 'cubic'):
        return si.interp1d(self.coors, potential, kind, copy = False, bounds_error = False, fill_value = potential[-1])

    def interpolate(self, potential, r , centre=None, kind = None):
        if centre is not None:
           if centre.any():
              coors = coors - centre
           r = norm_l2_along_axis(coors, axis=1)
        if isinstance(potential, RadialVector):
           potential = potential.values
        if kind is None:
           return np.interp(r, self.coors, potential, right=0)
        return self.interpolator(potential, kind)(r)
        #return np.interp(r, self.coors, potential, right=0)

    def get_midpoint_mesh(self, to=None):
        if to is None:
            to = len(self.coors)
        else:
            if not isinstance(to, int):
                to = self.get_r(to)
        if self.midpointMesh.has_key(to):
            return self.midpointMesh[to]

        if to is None:
            coors = self.coors
        else:
            coors = self.coors[0:to]

        midpoints = np.convolve(coors, [0.5, 0.5], 'valid')
        midmesh = ExplicitRadialMesh(midpoints)
        self.midpointMesh[to] = midmesh
        midmesh.parentMesh = self
        return midmesh

    def get_parent_mesh(self):
        return self.parentMesh

    def slice(self, x, y):
        if isinstance(x, float):
            x = self.get_index(x)
        if isinstance(y, float):
            y = self.get_index(y)
        return ExplicitRadialMesh(self.coors[x:y])

    def sparse_vector(self, vector):
        values = np.tile(float('NAN'), self.size)
        ii = self.coors.searchsorted(vector.mesh.get_coors())
        values[ii] = vector.values
        return RadialVector(self, values)

    def partition(self):
        at = np.where(self.coors[1:]<self.coors[:-1])[0]+1
        return [ExplicitRadialMesh(v) for v in np.split(self.coors, at)]


class RadialHyperbolicMesh(ExplicitRadialMesh):
    """   Radial hyperbolic mesh given by params, see __init__"""

    def __init__(self, jm, ap=None, size=None, from_zero=False):
        """
        If three (jp, ap, size) params is given, mesh has size points at  
         a(n) = ap * n  / ( jm - n)
        If two params given, mesh has ap points hyperbolicaly distributed over interval
           (0 if from_zero else 1, jm)
        If one params given, mesh has jm points hyperbolicaly distributed over interval
           (0 if from_zero else 1, jm)
        """
        if size is None:
            # range, number of points
            self.size = (ap if not ap is None else jm)
            self.ap = 1.0
            self.jm = self.size / jm + self.size
        else:
            # clasical
            self.size = size
            self.jm = jm
            self.ap = ap

        coors = np.arange((0.0 if from_zero else 1.0), self.size + 1,
                          dtype=np.float64)
        coors = self.ap * coors / (self.jm - coors)

        super(RadialHyperbolicMesh, self).__init__(np.asfortranarray(coors))
