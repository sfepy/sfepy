import numpy as nm
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
    argTypes = ('state',)
    regionMatters = True
    
    ##
    # 24.04.2007, c
    # 08.06.2007
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['state'], historySizes,
                            terms.dq_state_in_qp )
        
    ##
    # c: 24.04.2007, r: 15.01.2008
    def initData( self, key, ckey, **kwargs ):
        state, = self.getArgs( **kwargs )
        nFa, nQP = state.getDataShapes( ckey, kind = 'Surface' )[:2]
        shape = (nFa, nQP, state.dpn, 1)

        DataCache.initData( self, key, ckey, shape )

    ##
    # c: 24.04.2007, r: 15.01.2008
    def update( self, key, groupIndx, ih, **kwargs ):
        ckey = self.gToC( groupIndx )
        state, = self.getArgs( **kwargs )

        vec, indx = state()
        ap, sg = state.getApproximation( groupIndx, 'Surface' )
        sd = ap.surfaceData[groupIndx[1]]
        bf = ap.getBase( sd.faceType, 0, groupIndx[0] )
        self.function( self.data[key][ckey][ih], vec, indx.start,
                       bf, sd.econn )

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
    # c: 27.02.2007, r: 15.01.2008
    def initData( self, key, ckey, **kwargs ):
        state, = self.getArgs( **kwargs )

        nEl, nQP, dim = state.getDataShapes( ckey )[:3]
        sym = dim * (dim + 1) / 2
        shape = (nEl, nQP, sym, 1)

#        print self.name, key, ckey, shape
        DataCache.initData( self, key, ckey, shape )

    ##
    # c: 27.02.2007, r: 07.03.2008
    def update( self, key, groupIndx, ih, **kwargs ):
        ckey = self.gToC( groupIndx )
        if not self.valid['strain'][ckey]:
            state, = self.getArgs( **kwargs )

            vec, indx = state()
            ap, vg = state.getApproximation( groupIndx, 'Volume' )
            self.function( self.data[key][ckey][ih], vec, indx.start,
                           vg, ap.econn )
            isFinite = nm.isfinite( self.data[key][ckey][ih] )
            if not nm.alltrue( isFinite ):
                ii = nm.where( isFinite == False )
                print ii
                print self.data[key][ckey][ih][ii]
                print 'infinite strains in', ckey
#                from sfe.base.base import debug; debug()
                raise ValueError
            self.valid['strain'][ckey] = True
        if key == 'dstrain':
            if self.step > 0:
                self.data[key][ckey][ih] = self.data['strain'][ckey][ih] \
                                         - self.data['strain'][ckey][ih+1]
                if ih > 0:
                    print 'history update!'
                    print kwargs['history']
                    raise NotImplementedError
            else:
                self.data[key][ckey][ih].fill( 0.0 )
                
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
    # c: 12.03.2007, r: 14.01.2008
    def initData( self, key, ckey, **kwargs ):
        state, = self.getArgs( **kwargs )

        nEl, nQP, dim = state.getDataShapes( ckey )[:3]
        shape = (nEl, nQP, dim, 1)

#        print self.name, key, ckey, shape
        DataCache.initData( self, key, ckey, shape )

    ##
    # c: 12.03.2007, r: 14.01.2008
    def update( self, key, groupIndx, ih, **kwargs ):
        state, = self.getArgs( **kwargs )
        ap, vg = state.getApproximation( groupIndx, 'Volume' )
        ckey = self.gToC( groupIndx )

        vec, indx = state()
        self.function( self.data[key][ckey][ih], vec, indx.start, vg, ap.econn )

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
    # c: 13.03.2007, r: 17.01.2008
    def initData( self, key, ckey, **kwargs ):
        state, = self.getArgs( **kwargs )

        nEl, nQP = state.getDataShapes( ckey )[:2]
        shape = (nEl, nQP, 1, 1)

#        print self.name, key, ig, shape
        DataCache.initData( self, key, ckey, shape )

    ##
    # c: 13.03.2007, r: 17.01.2008
    def update( self, key, groupIndx, ih, **kwargs ):
        state, = self.getArgs( **kwargs )
        ap, vg = state.getApproximation( groupIndx, 'Volume' )
        ckey = self.gToC( groupIndx )

        vec, indx = state()
        self.function( self.data[key][ckey][ih], vec, indx.start, vg, ap.econn )

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


##
# c: 23.01.2008, r: 23.01.2008
class MatInQPDataCache( DataCache ):
    name = 'mat_in_qp'
    argTypes = ('mat', 'ap', 'assumedShapes', 'modeIn')

    ##
    # c: 23.01.2008, r: 23.01.2008
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames, ['matqp'], historySizes,
                            terms.dq_state_in_qp )
        self.shape = {}
        
    ##
    # c: 23.01.2008, r: 01.02.2008
    def initData( self, key, ckey, **kwargs ):
        mat, ap, assumedShapes, modeIn = self.getArgs( **kwargs )
        if modeIn is None:
            if mat.ndim == 3:
                modeIn = 'element_avg'
            elif mat.ndim == 2:
                ashape = assumedShapes[0]
                if ashape[2:] != mat.shape:
                    modeIn = 'vertex'
                else:
                    modeIn = 'const'
            else:
                raise ValueError

        shape = None
        for ashape in assumedShapes:
            if ashape[0] == 1:
                if ashape[1] == 1:
                    modeOut = 'const'
                else:
                    modeOut = 'const_in_qp'
            else:
                if ashape[1] == 1:
                    modeOut = 'variable'
                else:
                    modeOut = 'variable_in_qp'

            if modeIn == 'const':
                shape = ashape
                break
            elif modeIn == 'element_avg':
                if modeOut in ['variable_in_qp', 'variable']:
                    shape = ashape
                    break
            elif modeIn == 'vertex':
                if modeOut in ['variable_in_qp']:
                    shape = ashape
                    break

        if shape is None:
            raise ValueError

        self.modeIn = modeIn
        self.modeOut = modeOut
        self.shape[ckey] = shape
        DataCache.initData( self, key, ckey, shape )

    ##
    # c: 23.01.2008, r: 01.02.2008
    def update( self, key, groupIndx, ih, **kwargs ):
        import numpy as nm
        mat, ap, assumedShapes, modeIn = self.getArgs( **kwargs )
        ckey = self.gToC( groupIndx )

        shape = self.shape[ckey]
        if self.modeIn == 'const':
            mat2 = nm.reshape( mat.copy(), (1, 1) + mat.shape )
            mat2 = mat2.repeat( shape[1], 1 )
            matQP = mat2.repeat( shape[0], 0 )
            self.data[key][ckey][ih][:] = matQP

        elif self.modeIn == 'vertex':
            iname, ig = ckey[0], ckey[-1]
            
            gbf = ap.getBase( 'v', 0, iname, fromGeometry = True )
            group = ap.region.domain.groups[ig]
            conn = group.conn

            # dq_state_in_qp() works for vectors -> make a view of
            # shape (nEl, nQP, nRow * nCol, 1).
            vshape = shape[0:2] + (mat.shape[1], 1)
##             print self
##             print self.shape, ckey
##             print vshape
##             print self.data[key][ckey][ih].shape
##             from sfe.base.base import debug
##             debug()
            matQP = self.data[key][ckey][ih].reshape( vshape )
            self.function( matQP, mat, 0, gbf, conn )

        else:
            raise NotImplementedError

