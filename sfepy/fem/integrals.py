from sfepy.base.base import *
from quadratures import CustomQuadrature, quadratures

import re

_match_order_dim = re.compile( '.*_o([0-9]+)_d([0-9]+)$' ).match
_match_order_geom = re.compile( '.*_o([0-9]+)_g([0-9]+_[0-9]+)$' ).match

##
# 16.11.2007, c
class Integrals( Container ):

    ##
    # 16.11.2007, c
    def from_conf( conf, names ):
        objs = OneTypeList( Integral )

        name_map = {}
        for desc in conf.itervalues():
            name_map[desc.name] = desc

        for name in names:
            if not name_map.has_key( name ): continue

            int_conf = name_map[name]
            aux = Integral( name = int_conf.name,
                            kind = int_conf.kind,
                            quad_name = int_conf.quadrature,
                            mode = 'builtin' )
            if hasattr( int_conf, 'vals' ):
                aux.vals = nm.array( int_conf.vals, nm.float64 )
                aux.weights = nm.array( int_conf.weights, nm.float64 )
                aux.mode = 'custom'
                
            objs.append( aux )

        obj = Integrals( objs )
        return obj
    from_conf = staticmethod( from_conf )

    def set_quadratures( self, quadratures ):
        for intl in self:
            intl.qcs = quadratures[intl.quad_name]
            intl.setup()

##
# 16.11.2007, c
class Integral( Struct ):

    def setup( self ):
        if self.mode == 'builtin':
            match = _match_order_dim( self.quad_name )
            self.order, self.dim = [int( ii ) for ii in match.groups()]
            qcs = {}
            for key, quad_contructor in self.qcs.iteritems():
#                print key
                match = _match_order_geom( key )
                order = int( match.group( 1 ) )
                geom = match.group( 2 )
                qcs[geom] = quad_contructor
                assert_( order == self.order )
            self.qcs = qcs
        else:
            self.order = None
            self.dim = self.vals.shape[1]

    def create_qp( self ):
        self.qp = {}
        if self.mode == 'builtin':
            for geom, quad_contructor in self.qcs.iteritems():
                self.qp[geom] = quad_contructor()
        else:
            self.qp['custom'] = CustomQuadrature.from_conf( self )

    def get_qp( self, geometry ):
        if self.mode == 'builtin':
            return self.qp[geometry]
        else:
            return self.qp['custom']
