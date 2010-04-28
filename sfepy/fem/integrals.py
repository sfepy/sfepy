from sfepy.base.base import *
from quadratures import QuadraturePoints, quadrature_tables

import re

_match_order_dim = re.compile( '.*_o([0-9]+)_d([0-9]+)$' ).match

class Integrals(Container):
    """
    Container for instances of :class:`Integral`.
    """
    
    @staticmethod
    def from_conf(conf, names):
        objs = OneTypeList(Integral)

        name_map = {}
        for desc in conf.itervalues():
            name_map[desc.name] = desc

        for name in names:
            if not name_map.has_key(name): continue

            int_conf = name_map[name]

            if hasattr(int_conf, 'vals'):
                aux = Integral(int_conf.name,
                               kind=int_conf.kind,
                               quad_name=int_conf.quadrature,
                               coors=int_conf.vals,
                               weights=int_conf.weights)

            else:
                aux = Integral(int_conf.name,
                               kind=int_conf.kind,
                               quad_name=int_conf.quadrature)
                
            objs.append(aux)

        obj = Integrals(objs)
        return obj

class Integral(Struct):
    """
    Wrapper class around quadratures.
    """
    _msg1 = 'WARNING: quadrature order %d is not available for geometry %s!'
    _msg2 = 'WARNING: using %d instead!'
    
    def __init__(self, name, kind='v', quad_name='auto',
                 coors=None, weights=None, term=None):
        self.name = name
        self.kind = kind
        self.quad_name = quad_name
        self.term = term
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

    def set_term(self, term):
        self.term = term

    def get_qp(self, geometry):

        if geometry in self.qps:
            qp = self.qps[geometry]

        else:
            if self.mode == 'builtin':
                table = quadrature_tables[geometry]
                if self.order in table:
                    qp = table[self.order]

                else:
                    orders = table.keys()
                    ii = nm.searchsorted(orders, self.order)
                    if ii >= len(orders):
                        order = max(orders)
                        output(self._msg1 % (self.order, geometry))
                        output(self._msg2 % order)

                    else:
                        order = orders[ii]

                    qp = table[order]

            else:
                qp = QuadraturePoints(None,
                                      coors=self.coors, weights=self.weights)

            self.qps[geometry] = qp

        return qp.coors, qp.weights
