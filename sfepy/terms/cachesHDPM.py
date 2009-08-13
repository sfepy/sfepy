from sfepy.terms.extmods import terms
from sfepy.terms.cache import DataCache

##
# 11.06.2007, c
class HDPMDeltaStateDataCache( DataCache ):
    name = 'hdpm_dstate'
    argTypes = ('state1', 'state2', 'state10', 'state20',)

    ##
    # 11.06.2007, c
    def __init__( self, name, argNames, historySizes = None ):
        DataCache.__init__( self, name, argNames,
                            ['dstate_sub', 'dstate_add'],
                            historySizes, terms.state_in_qp )
        
    ##
    # 11.06.2007, c
    def initData( self, key, ig, **kwargs ):
        args = self.getArgs( **kwargs )
        state = args[0]
        
        sh = state.getDataShapes()
        nEl, nQP = sh[ig][:2]
        shape = (nEl, nQP, state.dpn, 1)

#        print self.name, key, ig, shape
        DataCache.initData( self, key, ig, shape )

    ##
    # 11.06.2007, c
    def update( self, key, ig, ih, **kwargs ):
        state1, state2, state10, state20 = self.getArgs( **kwargs )

        vec1, i1 = state1()
        vec2, i2 = state2()
        vec10, i10 = state10()
        vec20, i20 = state20()

        dvec1 = vec1[i1] - vec2[i2]
        dvec0 = vec10[i10] - vec20[i20]
        if key == 'dstate_sub':
            vec = dvec1 - dvec0
        else:
            vec = dvec1 + dvec0
    
        ap, vg = state1.getCurrentApproximation()
        self.function( self.data[key][ig][ih], vec, 0, ap.bf['v'], ap.econn )
