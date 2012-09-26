import numpy as np
import math
from scipy.integrate import simps
from sfepy.linalg import norm_l2_along_axis
import scipy.interpolate

class RadialVector(object):
     @staticmethod
     def ByXY(x,y):
        return RadialVector(ExplicitMesh(x), y)

     def __init__(self, mesh, values):
        if isinstance(mesh, np.ndarray):
           mesh=ExplicitRadialMesh(mesh)
        self.mesh=mesh
        self.values=values
        self.interpolated=None
        self.precision=None;
        self.resultPrecision=None
        self.running=None

     def running_mean(self):
         if self.running is None:
           window=5
           weights=np.array((1.0,3,4,3,1))
           weights=weights/weights.sum()
           wsize=int(weights.size-1)
           data = np.hstack([ np.ones(wsize)*self.values[0], self.values, np.ones(wsize/2)*self.values[-1] ])
           self.running = np.convolve(data, weights)[wsize:-(wsize/2)-wsize]
         return self.running

     def integrate(self, precision=0.0001):
         return self.mesh.integrate(self)

     def get_interpolated(self, precision=0.0001, grade=10):
        if precision is None:
           return scipy.interpolate.InterpolatedUnivariateSpline(self.mesh.get_coors(), self.values, k=5)

        if self.interpolated is None or (self.precision == self.resultPrecision and precision < self.precision):
           self.precision=precision
           data=self.runningMean()
           zpresnuji=None
           while True:
             prec=(data.mean()*precision)**2 * self.values.size;
             self.interpolated=scipy.interpolate.UnivariateSpline(self.mesh.get_coors(), data, k=5, s=precision)
             der=self.interpolated(self.mesh.get_coors(),1)
             sig=np.sign(der)
             if np.abs(sig[1:]-sig[:-1]).sum()<=grade:
                break
             precision=precision*2
        self.resultPrecision=precision
        a=self.interpolated
        return a

     def interpolated_values(self,at=None, precision=0.0001, grade=10):
         if at is None:
            at=self.mesh.get_coors()
         return self.getInterpolated(precision, grade)(at)

     def interpolated_derivatives(self,at=None, precision=0.0001):
         if at is None:
            at=self.mesh.get_coors()
         return self.get_interpolated(precision)(at,1)

     """ Treat vector as integral of spheres, return corresponding radial function"""
     def radial_derivatives(self):
         difference=np.convolve(self.values, [-1,1])
         factor=np.convolve(self.mesh.coors**2*math.pi, [-1,1])
         parent=self.mesh.getParentMesh()
         if parent is None:
           return RadialVector.ByXY(self.mesh.get_midpoint_mesh(), difference/factor)
         return RadialVector(parent, self.interpolated_values(parent.get_coors(), None))

     def slice(self, x, y):
        if isinstance(x, float):
           x=self.get_index(x)
        if isinstance(y, float):
           y=self.get_index(y)
        return RadialVector(self.mesh.slice(x,y), self.values[x:y])

     @staticmethod
     def fromfile(file):
         array=np.genfromtxt(file)
         mesh = ExplicitRadialMesh(array[:,0])
         return [ RadialVector(mesh, array[:,r]) for r in xrange(1, array.shape[1]) ]

     def tofile(self, filename=None):
         if filename is None:
            import sys
            filename=sys.stdout
         vector=np.vstack([self.mesh.get_coors(), self.values])
         np.savetxt(filename, vector.T)

     def extrapolate(self, x):
         return self.mesh.extrapolate(self.values, x)

     def extrapolate_3d(self, coors, centre = (0,0,0)):
         return self.mesh.extrapolate_3d(self.values, coors, centre)

     def output_vector(self, filename=None):
         return self.mesh.output_vector(self, filename)

     @staticmethod
     def SparseMerge(vectors):
        mesh = RadialMesh.Merge([v.mesh for v in vectors])
        return [mesh.sparse_vector(v) for v in vectors]

class RadialMesh(object):
    """
    Radial mesh.
    """

    def extrapolate_3d(self, potential, coors, centre = None):
        if not centre is None:
           coors = coors - centre
        r = norm_l2_along_axis(coors, axis=1)
        return self.extrapolate(potential, r)

    def integrate(self, vector):
        """
        .. math::
           \int f(r) r^2 dr
        """
        return simps(vector * self.coors**2, self.coors)

    def dot(self, vector_a, vector_b):
        """
        .. math::
           \int f(r) g(r) r^2 dr
        """
        return self.integrate(vector_a * vector_b)

    def norm(self, vector):
        return np.sqrt(self.vdot(vector, vector))

    def output_vector(self, vector, filename=None):
        if filename is None:
           import sys
           filename=sys.stdout
        if isinstance(vector, RadialVector):
           vector = [vector]
        if isinstance(vector[0], RadialVector):
           vector = [v.values for v in vector]
        vector=np.vstack([self.coors, vector])
        np.savetxt(filename, vector.T)

    @staticmethod
    def Merge(meshes):
        merged=np.concatenate(tuple(m.get_coors() for m in meshes))
        return ExplicitRadialMesh(np.unique(merged))

class ExplicitRadialMesh(RadialMesh):
    def __init__(self, coors):
        self.coors=coors
        self.midpointMesh={}
        self.parentMesh=None

    @property
    def shape(self):
        return self.coors.shape

    @property
    def size(self):
        return self.coors.size

    def get_coors(self):
        return self.coors

    def last_point(self):
        return self.coors[self.coors.size-1]

    def get_r(self, index):
        return self.coors[index]

    def get_index(self, r):
        pos = self.coors.searchsorted(r)
        return pos if pos < self.coors.size else self.coors.size -1

    def get_mixing(self, r):
        pos = self.get_index(r)
        if (pos == self.coors.size - 1) and (self.coors[pos] < r):
            out = [(pos, 1.0)]

        elif (pos == 0) or (self.coors[pos] == r):
            out = [(pos, 1.0)]

        else:
            pos_c = (r - self.coors[pos-1]) \
                    / (self.coors[pos] - self.coors[pos-1])
            out = [(pos - 1, 1.0 - pos_c), (pos, pos_c)]

        return out

    def extrapolate(self, potential, r):
        return np.interp(r, self.coors, potential, right = 0)

    def get_midpoint_mesh(self,to=None):
        if to is None:
          to=len(self.coors)
        else:
          if not isinstance(to, int):
             to=self.get_r(to);
        if self.midpointMesh.has_key(to):
          return self.midpointMesh[to]
        if to is None:
          coors=self.coors
        else:
          coors=self.coors[0:to]
        midpoints=np.convolve(coors, [0.5,0.5], "valid")
        midmesh=ExplicitRadialMesh(midpoints)
        self.midpointMesh[to]=midmesh
        midmesh.parentMesh=self
        return midmesh;

    def get_parent_mesh(self):
        return self.parentMesh

    def slice(self, x, y):
        if isinstance(x, float):
           x=self.get_index(x)
        if isinstance(y, float):
           y=self.get_index(y)
        return ExplicitRadialMesh(self.coors[x:y])

    def sparse_vector(self, vector):
         values=np.tile(float('NAN'), (self.size))
         values[self.coors.searchsorted(vector.mesh.get_coors())]=vector.values
         return RadialVector(self, values)



class RadialHyperbolicMesh(ExplicitRadialMesh):
    size = None

    def __init__(self, jm, ap = None, size = None, from_zero=False):
        if size is None:
          #range, number of points
          self.size = ap if not ap is None else jm
          self.ap = 1.0
          self.jm = self.size / jm + self.size
        else:
          #clasical
          self.size = size
          self.jm = jm
          self.ap = ap

        coors = np.arange(0.0 if from_zero else 1.0, self.size + 1, dtype=np.float64)
        coors = self.ap * coors / (self.jm - coors)
        super(RadialHyperbolicMesh, self).__init__( np.asfortranarray(coors))
