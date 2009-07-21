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


class MatInQPDataCache( DataCache ):
    name = 'mat_in_qp'
    arg_types = ('mat', 'ap', 'assumed_shapes', 'mode_in')

    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['matqp'], history_sizes,
                            terms.dq_state_in_qp )
        self.shape = {}
        self.mode_in = {}
        self.mode_out = {}
        
    def init_data( self, key, ckey, **kwargs ):
        mat, ap, assumed_shapes, mode_in = self.get_args( **kwargs )
        if mode_in is None:
            if mat.ndim == 3:
                ig = ckey[1]
                rshape = ap.region.shape[ig]
                if rshape.n_vertex == rshape.n_cell:
                    msg = ('cannot determine mode_in! (%d nodes, %d cells ' +
                           'material data shape: %s)') \
                           % (rshape.n_vertex, rshape.n_cell, mat.shape)
                    raise ValueError( msg )
                if mat.shape[0] == rshape.n_vertex:
                    mode_in = 'vertex'
                elif mat.shape[0] == rshape.n_cell:
                    mode_in = 'element_avg'
                else:
                    msg = ('cannot determine mode_in! (%d nodes, %d cells ' +
                           'material data shape: %s)') \
                           % (rshape.n_vertex, rshape.n_cell, mat.shape)
                    raise ValueError( msg )
            elif mat.ndim == 2:
                ashape = assumed_shapes[0]
                if ashape[2:] != mat.shape:
                    mode_in = 'vertex'
                else:
                    mode_in = 'const'
            else:
                raise ValueError

        if 'mode_out' in kwargs:
            mode_out = kwargs['mode_out']
            shape = assumed_shapes[0]

        else:
            shape = None
            for ashape in assumed_shapes:
                if ashape[0] == 1:
                    if ashape[1] == 1:
                        mode_out = 'const'
                    else:
                        mode_out = 'const_in_qp'
                else:
                    if ashape[1] == 1:
                        mode_out = 'variable'
                    else:
                        mode_out = 'variable_in_qp'

                if mode_in == 'const':
                    shape = ashape
                    break
                elif mode_in == 'element_avg':
                    if mode_out in ['variable_in_qp', 'variable']:
                        shape = ashape
                        break
                elif mode_in == 'vertex':
                    if mode_out in ['variable_in_qp']:
                        shape = ashape
                        break

        if shape is None:
            raise ValueError

        self.dtype = mat.dtype
        self.geometry = kwargs.get('geometry', 'volume')
        self.region_name = kwargs.get('region_name', None)
        
        self.mode_in[ckey] = mode_in
        self.mode_out[ckey] = mode_out
        self.shape[ckey] = shape
        DataCache.init_data( self, key, ckey, shape )

    def update( self, key, group_indx, ih, **kwargs ):
        import numpy as nm
        mat, ap, assumed_shapes, mode_in = self.get_args( **kwargs )
        ckey = self.g_to_c( group_indx )
        shape = self.shape[ckey]
        if self.mode_in[ckey] == 'const':
            mat2 = nm.reshape( mat.copy(), (1, 1) + mat.shape )
            mat2 = mat2.repeat( shape[1], 1 )
            mat_qp = mat2.repeat( shape[0], 0 )
            self.data[key][ckey][ih][:] = mat_qp

        elif self.mode_in[ckey] == 'vertex':
            """no group.lconn, so it is built here..."""
            iname, ig = ckey[0], ckey[-1]


            if self.geometry == 'volume':
                gbf = ap.get_base( 'v', 0, iname, from_geometry = True )
                group = ap.region.domain.groups[ig]
                conn = group.conn
            else:
                sd = ap.surface_data[self.region_name]
                gbf = ap.get_base(sd.face_type, 0, iname, from_geometry=True)
                conn = sd.leconn
                
            # dq_state_in_qp() works for vectors -> make a view of
            # shape (n_el, n_qp, n_row * n_col, 1).
            vshape = shape[0:2] + (nm.prod( mat.shape[1:] ), 1)
##             print self
##             print self.shape, ckey
##             print vshape
##             print self.data[key][ckey][ih].shape
##             debug()
            mat_qp = self.data[key][ckey][ih].reshape( vshape )
            if ((self.geometry == 'volume')
                and (ap.region.n_v_max > group.shape.n_vertex)):
                remap = nm.zeros( (ap.region.n_v_max,), dtype = nm.int32 )
                remap[group.vertices] = nm.arange( group.shape.n_vertex,
                                                   dtype = nm.int32 )
                lconn = remap[conn]
            else:
                lconn = conn

            if mat.dtype == nm.float64:
                self.function( mat_qp, mat, 0, gbf, lconn )
            else: # complex128
                ac = nm.ascontiguousarray
                val_r = nm.zeros_like(mat_qp.real)
                val_i = nm.zeros_like(mat_qp.imag)
                self.function(val_r, ac(mat.real), 0, gbf, lconn)
                self.function(val_i, ac(mat.imag), 0, gbf, lconn)
                mat_qp.flat[:] = val_r + 1j * val_i
        else:
            raise NotImplementedError
