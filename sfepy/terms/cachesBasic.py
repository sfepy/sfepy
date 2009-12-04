import numpy as nm
from sfepy.terms.extmods import terms
from sfepy.terms.cache import DataCache
from sfepy.base.base import pause, debug

class StateInVolumeQPDataCache( DataCache ):
    name = 'state_in_volume_qp'
    arg_types = ('state', 'get_vector')

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['state'], history_sizes,
                            terms.dq_state_in_qp )
        
    def init_data( self, key, ckey, **kwargs ):
        state, aux = self.get_args( **kwargs )

        n_el, n_qp = state.get_data_shapes( ckey )[:2]
        shape = (n_el, n_qp, state.dpn, 1)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        state, get_vector = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )
        ckey = self.g_to_c( group_indx )

        if ih == 0:
            bf = ap.get_base( 'v', 0, group_indx[0] )
            vec = get_vector( state )
            self.function( self.data[key][ckey][ih], vec, 0, bf, ap.econn )
        else:
            print 'history update!'
            print kwargs['history']
            raise NotImplementedError

class StateInSurfaceQPDataCache( DataCache ):
    name = 'state_in_surface_qp'
    arg_types = ('state',)
    region_matters = True
    
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['state'], history_sizes,
                            terms.dq_state_in_qp )
        
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )
        n_fa, n_qp = state.get_data_shapes( ckey, kind = 'Surface' )[:2]
        shape = (n_fa, n_qp, state.dpn, 1)

        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        ckey = self.g_to_c( group_indx )
        state, = self.get_args( **kwargs )

        ap, sg = state.get_approximation( group_indx, 'Surface' )
        sd = ap.surface_data[group_indx[1]]
        bf = ap.get_base( sd.face_type, 0, group_indx[0] )
        self.function( self.data[key][ckey][ih], state(), 0, bf, sd.econn )

class CauchyStrainDataCache( DataCache ):
    name = 'cauchy_strain'
    arg_types = ('state', 'get_vector')

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['strain'],
                            history_sizes, terms.dq_cauchy_strain )
        
    def init_data( self, key, ckey, **kwargs ):
        state, aux = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shapes( ckey )[:3]
        sym = dim * (dim + 1) / 2
        shape = (n_el, n_qp, sym, 1)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        ckey = self.g_to_c( group_indx )
        if ih > 0:
            print 'history update!'
            print kwargs['history']
            raise NotImplementedError
        state, get_vector = self.get_args( **kwargs )

        ap, vg = state.get_approximation( group_indx, 'Volume' )
        vec = get_vector( state )
        self.function( self.data['strain'][ckey][ih], vec, 0, vg, ap.econn )
        is_finite = nm.isfinite( self.data[key][ckey][ih] )
        if not nm.alltrue( is_finite ):
            ii = nm.where( is_finite == False )
            print ii
            print self.data[key][ckey][ih][ii]
            print 'infinite strains in', ckey
            raise ValueError
        self.valid['strain'][ckey] = True
                
class GradScalarDataCache( DataCache ):
    name = 'grad_scalar'
    arg_types = ('state',)

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['grad'], history_sizes,
                            terms.dq_grad )
        
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shapes( ckey )[:3]
        shape = (n_el, n_qp, dim, 1)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )
        ckey = self.g_to_c( group_indx )

        self.function( self.data[key][ckey][ih], state(), 0, vg, ap.econn )

class GradVectorDataCache( GradScalarDataCache ):
    name = 'grad_vector'
    arg_types = ('state',)

    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shapes( ckey )[:3]
        shape = (n_el, n_qp, dim, dim)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )

class DivVectorDataCache( DataCache ):
    name = 'div_vector'
    arg_types = ('state',)

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['div'], history_sizes,
                            terms.dq_div_vector )
        
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp = state.get_data_shapes( ckey )[:2]
        shape = (n_el, n_qp, 1, 1)

#        print self.name, key, ig, shape
        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )
        ckey = self.g_to_c( group_indx )

        self.function( self.data[key][ckey][ih], state(), 0, vg, ap.econn )

class VolumeDataCache( DataCache ):
    name = 'volume'
    arg_types = ('region','field')

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['volume'], history_sizes )
        
    def init_data( self, key, ckey, **kwargs ):
        shape = (1, 1, 1, 1)

        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        region, field = self.get_args( **kwargs )
        ckey = self.g_to_c( group_indx )
        self.data[key][ckey][ih] = region.get_volume( field, ckey,
                                                      update = True )

class SurfaceDataCache( DataCache ):
    name = 'surface'
    arg_types = ('region','field')

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['surface'], history_sizes )
         
    def init_data( self, key, ckey, **kwargs ):
        shape = (1, 1, 1, 1)

        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        region, field = self.get_args( **kwargs )
        ckey = self.g_to_c( group_indx )
        self.data[key][ckey][ih] = region.get_volume( field, ckey,
                                                      update = True,
                                                      mode = 'Surface' )
