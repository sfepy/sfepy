"""
Classes for constructing potentials of atoms and molecules.
"""
from sfepy.base.base import *
from sfepy.base.la import norm_l2_along_axis

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

        else:
            out = NotImplemented

        return out

    def __sub__(self, other):
        if isinstance(other, PotentialBase):
            out = CompoundPotential([self, -1.0 * other])

        else:
            out = NotImplemented

        return out

    def __pos__(self):
        return self

    def __neg__(self):
        out = -1.0 * self
        return out

class Potential(PotentialBase):
    """
    Single potential.
    """

    def __init__(self, name, function, centre=None, dim=3):
        self.name = name
        self.function = function

        if centre is None:
            centre = nm.array([0.0] * dim, dtype=nm.float64)

        self.centre = nm.asarray(centre, dtype=nm.float64)

        self.sign = 1.0

    def __call__(self, coors):
        r = norm_l2_along_axis(coors - self.centre)

        pot = self.sign * self.function(r)

        return pot
