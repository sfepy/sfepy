import extmods.terms as terms
from cache import DataCache

##
# 13.03.2007, c
class StateInVolumeQPDataCache( DataCache ):
    name = 'state_in_volume_qp'
    argTypes = ('state',)

    ##
    # 13.03.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['state'], historySizes,
                            terms.dq_state_in_qp )
        
    ##
    # created:       13.03.2007
    # last revision: 13.12.2007
    def initData( self, key, ckey, **kwargs ):
        state, = self.getArgs( **kwargs )

        nEl, nQP = state.getDataShapes( ckey )[:2]
        shape = (nEl, nQP, state.dpn, 1)

#        print self.name, key, ckey, shape
        DataCache.initData( self, key, ckey, shape )

    ##
    # created:       13.03.2007
    # last revision: 13.12.2007
    def update( self, key, groupIndx, ih, **kwargs ):
        state, = self.getArgs( **kwargs )
        ap, vg = state.getApproximation( groupIndx, 'Volume' )
        ckey = self.gToC( groupIndx )

        vec, indx = state()
        bf = ap.getBase( 'v', 0, groupIndx[0] )
        self.function( self.data[key][ckey][ih], vec, indx.start,
                       bf, ap.econn )


##
# 24.04.2007, c
class StateInSurfaceQPDataCache( DataCache ):
    name = 'state_in_surface_qp'
    argTypes = ('state', 'region')

    ##
    # 24.04.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['state'], historySizes,
                            terms.dq_state_in_qp )
        
    ##
    # 24.04.2007, c
    # 03.05.2007
    def initData( self, key, ig, **kwargs ):
        state, region = self.getArgs( **kwargs )
        sh = state.getDataShapes( surface = True, key = region.name, ig = ig )
        nFa, nQP = sh[ig][:2]
        shape = (nFa, nQP, state.dpn, 1)
##         print self.name, key, ig, shape, region.name

        DataCache.initData( self, key, ig, shape )

    ##
    # 24.04.2007, c
    def update( self, key, ig, ih, **kwargs ):
        state, region = self.getArgs( **kwargs )

        vec, indx = state()
        ap, sg = state.getCurrentApproximation( surface = True,
                                                key = region.name )
        sd = ap.surfaceData[region.name]
        self.function( self.data[key][ig][ih], vec, indx.start,
                       ap.bf[sd.faceType], sd.econn )

##
# 27.02.2007, c
class CauchyStrainDataCache( DataCache ):
    name = 'cauchy_strain'
    argTypes = ('state',)

    ##
    # 27.02.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['strain', 'dstrain'],
                            historySizes, terms.dq_cauchy_strain )
        
    ##
    # 27.02.2007, c
    def initData( self, key, ig, **kwargs ):
        state, = self.getArgs( **kwargs )

        sh = state.getDataShapes()
        nEl, nQP, dim = sh[ig][:3]
        sym = dim * (dim + 1) / 2
        shape = (nEl, nQP, sym, 1)

#        print self.name, key, ig, shape
        DataCache.initData( self, key, ig, shape )

    ##
    # 27.02.2007, c
    # 08.06.2007
    # 27.08.2007
    def update( self, key, ig, ih, **kwargs ):
        if not self.valid['strain'][ig]:
            state, = self.getArgs( **kwargs )

            vec, indx = state()
            ap, vg = state.getCurrentApproximation()
            self.function( self.data[key][ig][ih], vec, indx.start,
                           vg, ap.econn )
            self.valid['strain'][ig] = True
        if key == 'dstrain':
            if self.step > 0:
                self.data[key][ig][ih] = self.data['strain'][ig][ih] \
                                         - self.data['strain'][ig][ih+1]
                if ih > 0:
                    print 'history update!'
                    print kwargs['history']
                    raise NotImplementedError
            else:
                self.data[key][ig][ih].fill( 0.0 )
                
##
# 12.03.2007, c
class GradScalarDataCache( DataCache ):
    name = 'grad_scalar'
    argTypes = ('state',)

    ##
    # 12.03.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['grad'], historySizes,
                            terms.dq_grad_scalar )
        
    ##
    # 12.03.2007, c
    def initData( self, key, ig, **kwargs ):
        state, = self.getArgs( **kwargs )

        sh = state.getDataShapes()
        nEl, nQP, dim = sh[ig][:3]
        shape = (nEl, nQP, dim, 1)

#        print self.name, key, ig, shape
        DataCache.initData( self, key, ig, shape )

    ##
    # 12.03.2007, c
    # 27.08.2007
    def update( self, key, ig, ih, **kwargs ):
        state, = self.getArgs( **kwargs )

        vec, indx = state()
        ap, vg = state.getCurrentApproximation()
        self.function( self.data[key][ig][ih], vec, indx.start, vg, ap.econn )

##
# 13.03.2007, c
class DivVectorDataCache( DataCache ):
    name = 'div_vector'
    argTypes = ('state',)

    ##
    # 13.03.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['div'], historySizes,
                            terms.dq_div_vector )
        
    ##
    # 13.03.2007, c
    def initData( self, key, ig, **kwargs ):
        state, = self.getArgs( **kwargs )

        sh = state.getDataShapes()
        nEl, nQP = sh[ig][:2]
        shape = (nEl, nQP, 1, 1)

#        print self.name, key, ig, shape
        DataCache.initData( self, key, ig, shape )

    ##
    # 13.03.2007, c
    # 27.08.2007
    def update( self, key, ig, ih, **kwargs ):
        state, = self.getArgs( **kwargs )

        vec, indx = state()
        ap, vg = state.getCurrentApproximation()
        self.function( self.data[key][ig][ih], vec, indx.start, vg, ap.econn )

##
# 23.04.2007, c
class VolumeDataCache( DataCache ):
    name = 'volume'
    argTypes = ('region','field')

    ##
    # 23.04.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['volume'], historySizes )
        
    ##
    # 23.04.2007, c
    def initData( self, key, ckey, **kwargs ):
        shape = (1, 1, 1, 1)

        DataCache.initData( self, key, ckey, shape )

    ##
    # created:       23.04.2007
    # last revision: 13.12.2007
    def update( self, key, groupIndx, ih, **kwargs ):
        region, field = self.getArgs( **kwargs )
        ckey = self.gToC( groupIndx )
        self.data[key][ckey][ih] = region.getVolume( field, ckey, update = True )
