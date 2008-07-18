import numpy as nm
import extmods.terms as terms
from cache import DataCache
from sfepy.base.base import pause, debug

##
# 13.03.2007, c
class StateInVolumeQPDataCache( DataCache ):
    name = 'state_in_volume_qp'
    arg_types = ('state',)

    ##
    # 13.03.2007, c
    # 08.06.2007
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['state'], history_sizes,
                            terms.dq_state_in_qp )
        
    ##
    # created:       13.03.2007
    # last revision: 13.12.2007
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp = state.get_data_shapes( ckey )[:2]
        shape = (n_el, n_qp, state.dpn, 1)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )

    ##
    # c: 13.03.2007, r: 02.04.2008
    def update( self, key, group_indx, ih, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )
        ckey = self.g_to_c( group_indx )

        if ih == 0:
            vec, indx = state()
            bf = ap.get_base( 'v', 0, group_indx[0] )
            self.function( self.data[key][ckey][ih], vec, indx.start,
                           bf, ap.econn )
        else:
            print 'history update!'
            print kwargs['history']
            raise NotImplementedError

##
# 24.04.2007, c
class StateInSurfaceQPDataCache( DataCache ):
    name = 'state_in_surface_qp'
    arg_types = ('state',)
    region_matters = True
    
    ##
    # 24.04.2007, c
    # 08.06.2007
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['state'], history_sizes,
                            terms.dq_state_in_qp )
        
    ##
    # c: 24.04.2007, r: 15.01.2008
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )
        n_fa, n_qp = state.get_data_shapes( ckey, kind = 'Surface' )[:2]
        shape = (n_fa, n_qp, state.dpn, 1)

        DataCache.init_data( self, key, ckey, shape )

    ##
    # c: 24.04.2007, r: 15.01.2008
    def update( self, key, group_indx, ih, **kwargs ):
        ckey = self.g_to_c( group_indx )
        state, = self.get_args( **kwargs )

        vec, indx = state()
        ap, sg = state.get_approximation( group_indx, 'Surface' )
        sd = ap.surface_data[group_indx[1]]
        bf = ap.get_base( sd.face_type, 0, group_indx[0] )
        self.function( self.data[key][ckey][ih], vec, indx.start,
                       bf, sd.econn )

##
# 27.02.2007, c
class CauchyStrainDataCache( DataCache ):
    name = 'cauchy_strain'
    arg_types = ('state',)

    ##
    # 27.02.2007, c
    # 08.06.2007
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['strain', 'dstrain'],
                            history_sizes, terms.dq_cauchy_strain )
        
    ##
    # c: 27.02.2007, r: 15.04.2008
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shapes( ckey )[:3]
        sym = dim * (dim + 1) / 2
        shape = (n_el, n_qp, sym, 1)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )
        if key == 'dstrain': # dstrain uses strain
            DataCache.init_data( self, 'strain', ckey, shape )

    ##
    # c: 27.02.2007, r: 15.04.2008
    def update( self, key, group_indx, ih, **kwargs ):
        ckey = self.g_to_c( group_indx )
        if not self.valid['strain'][ckey]:
            if ih > 0:
                print 'history update!'
                print kwargs['history']
                raise NotImplementedError
            state, = self.get_args( **kwargs )

            vec, indx = state()
            ap, vg = state.get_approximation( group_indx, 'Volume' )
            self.function( self.data['strain'][ckey][ih], vec, indx.start,
                           vg, ap.econn )
            is_finite = nm.isfinite( self.data[key][ckey][ih] )
            if not nm.alltrue( is_finite ):
                ii = nm.where( is_finite == False )
                print ii
                print self.data[key][ckey][ih][ii]
                print 'infinite strains in', ckey
#                from sfepy.base.base import debug; debug()
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
    arg_types = ('state',)

    ##
    # 12.03.2007, c
    # 08.06.2007
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['grad'], history_sizes,
                            terms.dq_grad_scalar )
        
    ##
    # c: 12.03.2007, r: 14.01.2008
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shapes( ckey )[:3]
        shape = (n_el, n_qp, dim, 1)

#        print self.name, key, ckey, shape
        DataCache.init_data( self, key, ckey, shape )

    ##
    # c: 12.03.2007, r: 14.01.2008
    def update( self, key, group_indx, ih, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )
        ckey = self.g_to_c( group_indx )

        vec, indx = state()
        self.function( self.data[key][ckey][ih], vec, indx.start, vg, ap.econn )

##
# 13.03.2007, c
class DivVectorDataCache( DataCache ):
    name = 'div_vector'
    arg_types = ('state',)

    ##
    # 13.03.2007, c
    # 08.06.2007
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['div'], history_sizes,
                            terms.dq_div_vector )
        
    ##
    # c: 13.03.2007, r: 17.01.2008
    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp = state.get_data_shapes( ckey )[:2]
        shape = (n_el, n_qp, 1, 1)

#        print self.name, key, ig, shape
        DataCache.init_data( self, key, ckey, shape )

    ##
    # c: 13.03.2007, r: 17.01.2008
    def update( self, key, group_indx, ih, **kwargs ):
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )
        ckey = self.g_to_c( group_indx )

        vec, indx = state()
        self.function( self.data[key][ckey][ih], vec, indx.start, vg, ap.econn )

##
# 23.04.2007, c
class VolumeDataCache( DataCache ):
    name = 'volume'
    arg_types = ('region','field')

    ##
    # 23.04.2007, c
    # 08.06.2007
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['volume'], history_sizes )
        
    ##
    # 23.04.2007, c
    def init_data( self, key, ckey, **kwargs ):
        shape = (1, 1, 1, 1)

        DataCache.init_data( self, key, ckey, shape )

    ##
    # created:       23.04.2007
    # last revision: 13.12.2007
    def update( self, key, group_indx, ih, **kwargs ):
        region, field = self.get_args( **kwargs )
        ckey = self.g_to_c( group_indx )
        self.data[key][ckey][ih] = region.get_volume( field, ckey, update = True )


##
# c: 23.01.2008, r: 23.01.2008
class MatInQPDataCache( DataCache ):
    name = 'mat_in_qp'
    arg_types = ('mat', 'ap', 'assumed_shapes', 'mode_in')

    ##
    # c: 23.01.2008, r: 06.05.2008
    def __init__( self, name, arg_names, history_sizes = None ):
        DataCache.__init__( self, name, arg_names, ['matqp'], history_sizes,
                            terms.dq_state_in_qp )
        self.shape = {}
        self.mode_in = {}
        self.mode_out = {}

    ##
    # c: 23.01.2008, r: 06.05.2008
    def init_data( self, key, ckey, **kwargs ):
        mat, ap, assumed_shapes, mode_in = self.get_args( **kwargs )
        if mode_in is None:
            if mat.ndim == 3:
                ig = ckey[1]
                rshape = ap.region.shape[ig]
                if rshape.n_vertex == rshape.n_cell:
                    output( 'cannot determine mode_in! (%d nodes, %d cells '
                            'material data shape: %s)'\
                            % (rshape.n_vertex, rshape.n_cell, mat.shape) )
                    raise ValueError
                if mat.shape[0] == rshape.n_vertex:
                    mode_in = 'vertex'
                elif mat.shape[0] == rshape.n_cell:
                    mode_in = 'element_avg'
                else:
                    output( 'cannot determine mode_in! (%d nodes, %d cells '
                            'material data shape: %s)'\
                            % (rshape.n_vertex, rshape.n_cell, mat.shape) )
                    raise ValueError
            elif mat.ndim == 2:
                ashape = assumed_shapes[0]
                if ashape[2:] != mat.shape:
                    mode_in = 'vertex'
                else:
                    mode_in = 'const'
            else:
                raise ValueError

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

        self.mode_in[ckey] = mode_in
        self.mode_out[ckey] = mode_out
        self.shape[ckey] = shape
        DataCache.init_data( self, key, ckey, shape )

    ##
    # c: 23.01.2008, r: 06.05.2008
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
            
            gbf = ap.get_base( 'v', 0, iname, from_geometry = True )
            group = ap.region.domain.groups[ig]
            conn = group.conn

            # dq_state_in_qp() works for vectors -> make a view of
            # shape (n_el, n_qp, n_row * n_col, 1).
            vshape = shape[0:2] + (nm.prod( mat.shape[1:] ), 1)
##             print self
##             print self.shape, ckey
##             print vshape
##             print self.data[key][ckey][ih].shape
##             debug()
            mat_qp = self.data[key][ckey][ih].reshape( vshape )
            if ap.region.n_v_max > group.shape.n_vertex:
                remap = nm.zeros( (ap.region.n_v_max,), dtype = nm.int32 )
                remap[group.vertices] = nm.arange( group.shape.n_vertex,
                                                   dtype = nm.int32 )
                lconn = remap[conn]
            else:
                lconn = conn
            self.function( mat_qp, mat, 0, gbf, lconn )

        else:
            raise NotImplementedError

