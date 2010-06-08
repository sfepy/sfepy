from sfepy.base.base import *
from quadratures import QuadraturePoints, quadrature_tables

import re

_match_order_dim = re.compile( '.*_o([0-9]+)_d([0-9]+)$' ).match

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
                               quad_name=desc.quadrature,
                               coors=desc.vals,
                               weights=desc.weights)

            else:
                aux = Integral(desc.name,
                               kind=desc.kind,
                               quad_name=desc.quadrature)
                
            objs.append(aux)

        obj = Integrals(objs)
        return obj

    def get(self, name, dim, kind='v'):
        """
        Return existing or new integral.

        Parameters
        ----------
        name : str
            The name can either be a string representation of a
            non-negative integer (the integral order) or 'a' (automatic
            order) or a string beginning with 'i' (existing custom
            integral name).
        """
        if name == 'a':
            raise NotImplementedError

        elif name[0] == 'i':
            try:
                obj = self[name]

            except IndexError:
                raise ValueError('integral %s is not defined!' % name)

        else:
            try:
                order = int(name)

            except:
                raise ValueError('unsupported integral reference! (%s)' % name)

            name = '__%s_o%s_d%d' % (kind, name, dim)
            if self.has_key(name):
                obj = self[name]

            else:
                # Create new integral, and add it to self.
                quad_name = '_o%d_d%d' % (order, dim)
                obj = Integral(name, kind, quad_name=quad_name)

                self.append(obj)

        return obj
    

class Integral(Struct):
    """
    Wrapper class around quadratures.
    """
    _msg1 = 'WARNING: quadrature order %d is not available for geometry %s!'
    _msg2 = 'WARNING: using %d instead!'
    
    def __init__(self, name, kind='v', quad_name='auto',
                 coors=None, weights=None):
        self.name = name
        self.kind = kind
        self.quad_name = quad_name
        self.qps = {}

        if coors is None:
            self.mode = 'builtin'

        else:
            self.mode = 'custom'
            self.coors = coors
            self.weights = weights

        if self.quad_name in ('auto', 'custom'):
            self.order = -1

        else:
            match = _match_order_dim(self.quad_name)
            self.order, self.dim = [int( ii ) for ii in match.groups()]

    def get_actual_order(self, geometry):
        """
        Return the actual integration order for given geometry.

        Parameters
        ----------
        geometry : str
            The geometry key describing the integration domain,
            see the keys of `sfepy.fem.quadratures.quadrature_tables`.

        Returns
        -------
        order : int
            If `self.order` is in quadrature tables it is this
            value. Otherwise it is the closest higher order. If no
            higher order is available, a warning is printed and the
            highest available order is used.
        """
        table = quadrature_tables[geometry]
        if self.order in table:
            order = self.order

        else:
            orders = table.keys()
            ii = nm.searchsorted(orders, self.order)
            if ii >= len(orders):
                order = max(orders)
                output(self._msg1 % (self.order, geometry))
                output(self._msg2 % order)

            else:
                order = orders[ii]

        return order

    def get_qp(self, geometry):
        """
        Get quadrature point coordinates and corresponding weights for given
        geometry.  For built-in quadratures, the integration order is
        given by `self.get_actual_order()`.

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
                order = self.get_actual_order(geometry)
                qp = quadrature_tables[geometry][order]

            else:
                qp = QuadraturePoints(None,
                                      coors=self.coors, weights=self.weights)

            self.qps[geometry] = qp

        return qp.coors, qp.weights

    def integrate(self, function, order=None, geometry=None):
        """
        Integrate numerically a given scalar function.

        Parameters
        ----------
        function : callable(coors)
            The function of space coordinates to integrate.
        order : int, optional
            The integration order. Default is given by
            `self.get_actual_order()`.
        geometry : str
            The geometry key describing the integration domain. Default
            is `'1_2'`, i.e. a line integral in [0, 1]. For other values
            see the keys of `sfepy.fem.quadratures.quadrature_tables`.

        Returns
        -------
        val : float
            The value of the integral.
        """
        if geometry is None:
            geometry = '1_2'

        table = quadrature_tables[geometry]

        if order is None:
            order = sorted(table.keys())[0]

        qp = table[order]
        
        fvals = function(qp.coors)
        val = nm.sum(fvals * qp.weights)

        return val
