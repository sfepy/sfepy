import numpy as nm

from sfepy.base.base import OneTypeList, Container, Struct, basestr
from quadratures import QuadraturePoints

class Integrals(Container):
    """
    Container for instances of :class:`Integral`.
    """
    
    @staticmethod
    def from_conf(conf):
        objs = OneTypeList(Integral)

        for desc in conf.itervalues():
            if hasattr(desc, 'vals'):
                aux = Integral(desc.name,
                               kind=desc.kind,
                               coors=desc.vals,
                               weights=desc.weights)

            else:
                aux = Integral(desc.name,
                               kind=desc.kind,
                               order=desc.order)

            objs.append(aux)

        obj = Integrals(objs)
        return obj

    def get(self, name, kind='v'):
        """
        Return existing or new integral.

        Parameters
        ----------
        name : str
            The name can either be a non-negative integer, a string
            representation of a non-negative integer (the integral
            order) or 'a' (automatic order) or a string beginning with
            'i' (existing custom integral name).
        """
        if name == 'a':
            raise NotImplementedError

        elif isinstance(name, basestr) and (name[0] == 'i'):
            try:
                obj = self[name]

            except IndexError:
                raise ValueError('integral %s is not defined!' % name)

        else:
            try:
                order = int(name)

            except:
                raise ValueError('unsupported integral reference! (%s)' % name)

            name = '__%s_o%d' % (kind, order)
            if self.has_key(name):
                obj = self[name]

            else:
                # Create new integral, and add it to self.
                obj = Integral(name, kind, order=order)

                self.append(obj)

        return obj
    

class Integral(Struct):
    """
    Wrapper class around quadratures.
    """

    def __init__(self, name, kind='v', order=1, coors=None, weights=None):
        self.name = name
        self.kind = kind
        self.qps = {}

        if coors is None:
            self.mode = 'builtin'

        else:
            self.mode = 'custom'
            self.coors = coors
            self.weights = weights

        self.order = 0

        if order in ('auto', 'custom', 'a', 'c'):
            self.order = -1

        else:
            self.order = int(order)

    def get_key(self):
        """
        Get the key string corresponding to the integral kind and order,
        that can be used to distinguish various cached data evaluated
        using the integral.
        """
        return '__%s_o%d' % (self.kind, self.order)

    def get_qp(self, geometry):
        """
        Get quadrature point coordinates and corresponding weights for given
        geometry. For built-in quadratures, the integration order is given by
        `self.order`.

        Parameters
        ----------
        geometry : str
            The geometry key describing the integration domain,
            see the keys of `sfepy.fem.quadratures.quadrature_tables`.

        Returns
        -------
        coors : array
           The coordinates of quadrature points.
        weights: array
           The quadrature weights.
        """

        if geometry in self.qps:
            qp = self.qps[geometry]

        else:
            if self.mode == 'builtin':
                qp = QuadraturePoints.from_table(geometry, self.order)

            else:
                qp = QuadraturePoints(None,
                                      coors=self.coors, weights=self.weights)

            self.qps[geometry] = qp

        return qp.coors, qp.weights

    def integrate(self, function, order=1, geometry='1_2'):
        """
        Integrate numerically a given scalar function.

        Parameters
        ----------
        function : callable(coors)
            The function of space coordinates to integrate.
        order : int, optional
            The integration order. For tensor product geometries, this is the
            1D (line) order.
        geometry : str
            The geometry key describing the integration domain. Default
            is `'1_2'`, i.e. a line integral in [0, 1]. For other values
            see the keys of `sfepy.fem.quadratures.quadrature_tables`.

        Returns
        -------
        val : float
            The value of the integral.
        """
        qp = QuadraturePoints.from_table(geometry, order)

        fvals = function(qp.coors)
        val = nm.sum(fvals * qp.weights)

        return val
