from sfe.base.base import *
from quadratures import CustomQuadrature, quadratures

import re

_matchOrderDim = re.compile( '.*_o([0-9]+)_d([0-9]+)$' ).match
_matchOrderGeom = re.compile( '.*_o([0-9]+)_g([0-9]+_[0-9]+)$' ).match

##
# 16.11.2007, c
class Integrals( Container ):

    ##
    # 16.11.2007, c
    def fromConf( conf, names ):
        objs = OneTypeList( Integral )

        nameMap = {}
        for desc in conf.itervalues():
            nameMap[desc.name] = desc

        for name in names:
            if not nameMap.has_key( name ): continue

            intConf = nameMap[name]
            aux = Integral( name = intConf.name,
                            kind = intConf.kind,
                            quadName = intConf.quadrature,
                            mode = 'builtin' )
            if hasattr( intConf, 'vals' ):
                aux.vals = nm.array( intConf.vals, nm.float64 )
                aux.weights = nm.array( intConf.weights, nm.float64 )
                aux.mode = 'custom'
                
            objs.append( aux )

        obj = Integrals( objs )
        return obj
    fromConf = staticmethod( fromConf )

    def setQuadratures( self, quadratures ):
        for intl in self:
            intl.qcs = quadratures[intl.quadName]
            intl.setup()

##
# 16.11.2007, c
class Integral( Struct ):

    def setup( self ):
        if self.mode == 'builtin':
            match = _matchOrderDim( self.quadName )
            self.order, self.dim = [int( ii ) for ii in match.groups()]
            qcs = {}
            for key, quadContructor in self.qcs.iteritems():
#                print key
                match = _matchOrderGeom( key )
                order = int( match.group( 1 ) )
                geom = match.group( 2 )
                qcs[geom] = quadContructor
                assert order == self.order
            self.qcs = qcs
        else:
            self.order = None
            self.dim = self.vals.shape[1]

    def createQP( self ):
        self.qp = {}
        if self.mode == 'builtin':
            for geom, quadContructor in self.qcs.iteritems():
                self.qp[geom] = quadContructor()
        else:
            self.qp['custom'] = CustomQuadrature.fromConf( self )

    def getQP( self, geometry ):
        if self.mode == 'builtin':
            return self.qp[geometry]
        else:
            return self.qp['custom']
