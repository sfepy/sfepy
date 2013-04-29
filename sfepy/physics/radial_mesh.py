import numpy as np
import math
from scipy.integrate import simps, cumtrapz
from sfepy.linalg import norm_l2_along_axis
import scipy.interpolate as si

class RadialVector(object):

    @staticmethod
    def from_file(file):
        array = np.genfromtxt(file)
        mesh = ExplicitRadialMesh(array[:, 0])
        return [RadialVector(mesh, array[:, r])
                for r in xrange(1, array.shape[1])]

    def __repr__(self):
        return self.pretty(10)

    def pretty(self, values):
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

    def to_file(self, filename=None):
        if filename is None:
            import sys
            filename = sys.stdout
        vector = np.vstack([self.mesh.get_coors(), self.values])
        np.savetxt(filename, vector.T)

    def running_mean(self):
        if self.running is None:
            weights = np.array((1.0, 3, 4, 3, 1))
            weights = weights / weights.sum()
            wsize = int(weights.size - 1)
            data = np.hstack([np.ones(wsize) * self.values[0], self.values,
                              np.ones(wsize / 2) * self.values[-1]])
            self.running = np.convolve(data, weights)[wsize:-(wsize / 2)
                                                      - wsize]
        return self.running

    def integrate(self, precision=0.0001):
        return self.mesh.integrate(self.values)

    def linear_integrate(self):
        return self.mesh.linear_integrate(self.values)

    def linear_integral(self, from_zero=False):
        return self.mesh.linear_integral(self.values, from_zero)

    def get_extrapolated(self, precision=0.0001, grade=10, attempts=10):
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
        if at is None:
            at = self.mesh
        elif not isinstance(at, RadialMesh):
            at = ExplicitRadialMesh(at)
        val = self.get_extrapolated(precision, grade, attempts)(at.get_coors())
        return RadialVector(at, val)

    def extrapolated_derivatives(self, at=None, precision=0.0001, attempts=10):
        if at is None:
            at = self.mesh
        elif not isinstance(at, RadialMesh):
            at = ExplicitRadialMesh(at)
        val = self.get_extrapolated(precision=0.0001, grade=10,
                                    attempts=10)(at.get_coors(), 1)
        return RadialVector(at, val)

    def derivatives(self, radial=True):
        if radial:
            factor = np.convolve(self.mesh.coors ** 2 * math.pi, [1, -1])[1:-1]
        else:
            factor = self.mesh.intervals()
        diffs = np.convolve(self.values, [1, -1])[1:-1]
        diffs /= factor
        out = np.convolve(diffs, [0.5, 0.5])
        out[0] = diffs[0]
        out[-1] = diffs[-1]
        return RadialVector(self.mesh, out)

    def linear_derivatives(self):
        return self.derivatives(False)

    def slice(self, x, y):
        if isinstance(x, float):
            x = self.get_index(x)
        if isinstance(y, float):
            y = self.get_index(y)
        return RadialVector(self.mesh.slice(x, y), self.values[x:y])

    def interpolate(self, x):
        return self.mesh.interpolate(self.values, x)

    def interpolate_3d(self, coors, centre=(0, 0, 0)):
        return self.mesh.interpolate_3d(self.values, coors, centre)

    def output_vector(self, filename=None):
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
        return self.mesh.plot(self)

    def get_coors(self):
        return self.mesh.get_coors()

    def __call__(self, name, *args, **kwargs):
        return getattr(self.values, name)(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self.values, name)

class RadialMesh(object):
    """
    Radial mesh.
    """

    def interpolate_3d(self, potential, coors, centre=None):
        if not centre is None:
            coors = coors - centre
        r = norm_l2_along_axis(coors, axis=1)
        return self.interpolate(potential, r)

    def integrate(self, vector, norm):
        """
        .. math::
           \int f(r) r^2 dr
        """
        if norm == 'r2':
           v = vector * self.coors ** 2
        elif norm == 'linear':
           v = vector
        elif norm == 'spherical' or True:
           v = vector * self.coors ** 2 * 4 * math.pi
        return simps(v, self.coors)

    def linear_integrate(self, vector):
        """
        .. math::
           \int f(r) dr
        """
        return simps(vector, self.coors)

    def linear_integral(self, vector, from_zero=False):
        """
        .. math::
          a_n = \int_{r_0}^{r_n} f(r) dr

        from_zero starts to integrate from zero, instead of starting between
        the first two points
        """
        r = self.get_coors()
        if from_zero:
            v = cumtrapz(vector, r, initial=vector[0] / 2 * max(0.0, r[0]))
        else:
            v = cumtrapz(vector, r)
            r = r[1:]
        return RadialVector(r, v)

    def dot(self, vector_a, vector_b, norm='spherical'):
        """
        .. math::
           \int f(r) g(r) r^2 dr
        """
        return self.integrate(vector_a * vector_b, norm)

    def norm(self, vector, norm = 'spherical'):
        return np.sqrt(self.dot(vector, vector, norm))

    def output_vector(self, vector, filename=None):
        if filename is None:
            import sys
            filename = sys.stdout
        if isinstance(vector, RadialVector):
            vector = [vector]
        if isinstance(vector[0], RadialVector):
            vector = [v.values for v in vector]
        vector = np.vstack([self.coors, vector])
        np.savetxt(filename, vector.T)

    def plot(self, vector, cmd='plot'):
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
        merged = np.concatenate(tuple(m.get_coors() for m in meshes))
        return ExplicitRadialMesh(np.unique(merged))

    def intervals(self):
        return np.convolve(self.coors, [1, -1])[1:-1]

class ExplicitRadialMesh(RadialMesh):

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

    def interpolate(self, potential, r):
        return np.interp(r, self.coors, potential, right=0)

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


class RadialHyperbolicMesh(ExplicitRadialMesh):
    size = None

    def __init__(self, jm, ap=None, size=None, from_zero=False):
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
