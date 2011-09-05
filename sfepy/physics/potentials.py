"""
Classes for constructing potentials of atoms and molecules.
"""
import numpy as nm

from sfepy.base.base import as_float_or_complex, Container, Struct
from sfepy.linalg import norm_l2_along_axis

class CompoundPotential(Container):
    """
    Sum of several potentials.
    """

    def __init__(self, objs=None):
        Container.__init__(self, objs=objs)

        self.update_expression()

    def insert(self, ii, obj):
        Container.insert(self, ii, obj)
        self.update_expression()

    def append(self, obj):
        Container.append(self, obj)
        self.update_expression()

    def update_expression(self):
        self.expression = []
        for pot in self:
            aux = [pot.sign, pot.name, pot.centre]
            self.expression.append(aux)

    def __mul__(self, other):
        out = CompoundPotential()
        for name, pot in self.iteritems():
            out.append(pot * other)

        return out

    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        if isinstance(other, PotentialBase):
            out = self.copy()
            out.append(other)

        elif isinstance(other, CompoundPotential):
            out = CompoundPotential(self._objs + other._objs)

        else:
            raise ValueError('cannot add CompoundPotential with %s!' % other)

        return out

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, PotentialBase):
            out = self + (-other)

        elif isinstance(other, CompoundPotential):
            out = self + (-other)

        else:
            raise ValueError('cannot subtract CompoundPotential with %s!' \
                             % other)

        return out

    def __rsub__(self, other):
        return -self + other

    def __pos__(self):
        return self

    def __neg__(self):
        return -1.0 * self

    def __call__(self, coors):
        val = 0.0
        for pot in self:
            val += pot(coors)

        return val

class PotentialBase(Struct):
    """
    Base class for potentials.
    """

    def __mul__(self, other):
        try:
            mul = as_float_or_complex(other)

        except ValueError:
            raise ValueError('cannot multiply PotentialBase with %s!' % other)

        out = self.copy(name=self.name)

        out.sign = mul * self.sign

        return out

    def __rmul__(self, other):
        return self * other

    def __add__(self, other):
        if isinstance(other, PotentialBase):
            out = CompoundPotential([self, other])

        elif nm.isscalar(other):
            if other == 0:
                out = self

            else:
                out = NotImplemented

        else:
            out = NotImplemented

        return out

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, PotentialBase):
            out = CompoundPotential([self, -1.0 * other])

        elif nm.isscalar(other):
            if other == 0:
                out = self

            else:
                out = NotImplemented

        else:
            out = NotImplemented

        return out

    def __rsub__(self, other):
        return -self + other

    def __pos__(self):
        return self

    def __neg__(self):
        out = -1.0 * self
        return out

class Potential(PotentialBase):
    """
    Single spherically symmetric potential.
    """

    def __init__(self, name, function, centre=None, dim=3, args=None):
        self.name = name
        self.function = function
        self.args = args if args is not None else ()

        if centre is None:
            centre = nm.array([0.0] * dim, dtype=nm.float64)

        self.centre = nm.asarray(centre, dtype=nm.float64)

        self.sign = 1.0

    def __call__(self, coors):
        r = self.get_distance(coors)

        pot = self.sign * self.function(r, *self.args)

        return pot

    def __iter__( self ):
        """
        Allow iteration even over a single potential.
        """
        yield self

    def __len__(self):
        """
        Allow length even of a single potential.
        """
        return 1

    def get_distance(self, coors):
        """
        Get the distance of points with coordinates `coors` of the
        potential centre.
        """
        return norm_l2_along_axis(coors - self.centre)

    def get_charge(self, coors, eps=1e-6):
        """
        Get charge corresponding to the potential by numerically
        applying Laplacian in spherical coordinates.
        """
        r = self.get_distance(coors)

        args = self.args

        f0 = self.function(r, *args)
        fp1 = self.function(r + eps, *args)
        fp2 = self.function(r + 2.0 * eps, *args)
        fm1 = self.function(r - eps, *args)
        fm2 = self.function(r - 2.0 * eps, *args)

        # Second derivative w.r.t. r.
        d2 = (fp2 - 2.0 * f0 + fm2) / (4.0 * eps * eps)
        # First derivative w.r.t. r.
        d1 = (fp1 - fm1) / (2.0 * eps)

        charge = - self.sign / (4.0 * nm.pi) * (d2 + 2.0 * d1 / r)

        return charge
