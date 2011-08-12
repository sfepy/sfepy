import numpy as nm
from scipy.integrate import simps

from sfepy.base.base import Struct
from sfepy.linalg import norm_l2_along_axis

class RadialMesh(Struct):
    """
    Radial mesh.
    """

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

    def extrapolate(self, r, potential):
        return nm.interp(r, self.coors, potential, right = 0)

    def extrapolate3D(self, coors, potential):
        r = norm_l2_along_axis(coors, axis=1)
        return self.extrapolate(r, potential)

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
        return nm.sqrt(self.vdot(vector, vector))

class RadialHyperbolicMesh(RadialMesh):

    def __init__(self, jm, ap, size):
        self.size = size
        self.jm = jm
        self.ap = ap

        self.coors = nm.arange(1.0, self.size + 1, dtype=nm.float64)
        self.coors = self.ap * self.coors / (self.jm - self.coors)
        self.coors = nm.asfortranarray(self.coors)

    def last_point(self):
        return self.coors[self.coors.size-1]
