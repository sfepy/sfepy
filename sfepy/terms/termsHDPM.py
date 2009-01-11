from terms import *
from termsLaplace import DiffusionVelocityTerm

class HDPMGDPTerm( Term ):
    name = 'dw_hdpm_gdp'
    arg_types = ('ts', 'material', 'virtual', 'state_1', 'state_2',
                'parameter_1', 'parameter_2')
    geometry = [(Volume, 'virtual')]
    use_caches = {'hdpm_dstate' : [['state_1', 'state_2',
                                   'parameter_1', 'parameter_2']]}
##     use_caches = {'state_in_volume_qp' : [['state_1', {'state' : (2,2)}],
##                                          ['state_2', {'state' : (2,2)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_hdpm_g )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat, virtual, state1, state2, state10, state20\
            = self.get_args( **kwargs )
        ap, vg = virtual.get_current_approximation()
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state_1' ):
            shape = (chunk_size, 1, n_ep, n_ep)
            mode = 1
        elif diff_var == self.get_arg_name( 'state_2' ):
            shape = (chunk_size, 1, n_ep, n_ep)
            mode = 2
        else:
            raise StopIteration

        cache = self.get_cache( 'hdpm_dstate', 0 )
        ddp12 = cache( 'dstate_sub', self.char_fun.ig, 0,
                       state1 = state1, state2 = state2,
                       state10 = state10, state20 = state20 )

        if ts.step == 0:
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                yield out, chunk, 0
            return

        # Which is faster??
        
##         cache1 = self.get_cache( 'state_in_volume_qp', 0 )
##         cache2 = self.get_cache( 'state_in_volume_qp', 1 )
##         p1 = cache1( 'state', self.char_fun.ig, 0, state = state1 )
##         p2 = cache2( 'state', self.char_fun.ig, 0, state = state2 )
##         p10 = cache1( 'state', self.char_fun.ig, 1, state = state10 )
##         p20 = cache2( 'state', self.char_fun.ig, 1, state = state20 )
##         ddp12 = (p1 - p2) - (p10 - p20)

        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
#        print mat_qp
        for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
            status = self.function( out, 1.0 / ts.dt, ddp12, ap.bf['v'],
                                    mat_qp, vg, chunk, mode )
            yield out, chunk, status

class HDPMGPTerm( Term ):
    name = 'dw_hdpm_gp'
    arg_types = ('material', 'virtual', 'state_1', 'state_2')
    geometry = [(Volume, 'virtual')]
    use_caches = {'state_in_volume_qp' : [['state_1'], ['state_2']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_hdpm_g )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state1, state2 = self.get_args( **kwargs )
        ap, vg = virtual.get_current_approximation()
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state_1' ):
            shape = (chunk_size, 1, n_ep, n_ep)
            mode = 1
        elif diff_var == self.get_arg_name( 'state_2' ):
            shape = (chunk_size, 1, n_ep, n_ep)
            mode = 2
        else:
            raise StopIteration

        cache1 = self.get_cache( 'state_in_volume_qp', 0 )
        cache2 = self.get_cache( 'state_in_volume_qp', 1 )
        p1 = cache1( 'state', self.char_fun.ig, 0, state = state1 )
        p2 = cache2( 'state', self.char_fun.ig, 0, state = state2 )
        dp = p1 - p2

        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
        for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
            status = self.function( out, 1.0, dp, ap.bf['v'],
                                    mat_qp, vg, chunk, mode )
            yield out, chunk, status

class HDPMDiffusionVelocitySIntegratedTerm( Term ):
    name = 'd_hdpm_surfdvel'
    arg_types = ('material', 'parameter')
    geometry = [(SurfaceExtra, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_hdpm_surfdvel )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """? move surface pressure grad part into a cache ?"""
        mat, par = self.get_args( **kwargs )
        ap, sg = par.get_approximation( self.get_current_group(),
                                        'SurfaceExtra' )
        n_fa, n_qp = ap.get_s_data_shape( self.integral_name,
                                          self.region.name )[:2]
        shape = (chunk_size, 1, 1, 1)

        sd = ap.surface_data[self.region.name]

        vec = par()
        mat = nm.asarray( mat, dtype = nm.float64  )
        mat_qp = mat[nm.newaxis,:,:].repeat( n_qp, 0 )
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            status = self.function( out, vec, 0,
                                    mat_qp, sg, sd.fis, lchunk, ap.econn, chunk )
            out1 = nm.sum( out )
            yield out1, chunk, status

class HDPMPressureGradientTerm( DiffusionVelocityTerm ):
    """Pressure gradient."""
    name = 'de_hdpm_pgrad'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        parameter, = self.get_args( **kwargs )
        ap, vg = parameter.get_current_approximation()
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape()

        if diff_var is None:
            shape = (chunk_size, 1, dim, 1)
        else:
            raise StopIteration

        vec = parameter()
        mat_qp = nm.eye( dim )[nm.newaxis,:,:].repeat( n_qp, 0 )
#        print mat_qp
        for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
            status = self.function( out, vec, 0,
                                    mat_qp, vg, ap.econn, chunk )
            yield out, chunk, status

class HDPMConvRTerm( Term ):
    name = 'dw_hdpm_bh'
    arg_types = ('ts', 'material', 'virtual', 'state_1', 'state_2', 'history')
    geometry = [(Volume, 'virtual')]
    use_caches = {'state_in_volume_qp' : [['state_1', {'state' : (-1,-1)}],
                                         ['state_2', {'state' : (-1,-1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_hdpm_b1 )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state1, state2, history\
            = self.get_args( **kwargs )
        apr, vgr = virtual.get_current_approximation()
        apc, vgc = state1.get_current_approximation()
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, dim * n_epr, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state_1' ):
            n_epc = apc.get_v_data_shape()[3]
            shape = (chunk_size, 1, dim * n_epr, n_epc)
            mode = 1
        elif diff_var == self.get_arg_name( 'state_2' ):
            n_epc = apc.get_v_data_shape()[3]
            shape = (chunk_size, 1, dim * n_epr, n_epc)
            mode = 2
        else:
            raise StopIteration

        if ts.step == 0:
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                yield out, chunk, 0
            return

        if mode == 1:
            mat_qp = mats[0][nm.newaxis,:,:].repeat( n_qp, 0 )
##             print mat_qp
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
                status = self.function( out, 0.5, nm.empty( 0 ),
                                        apc.bf['v'], mat_qp, vgr, chunk, 1 )
                yield out, chunk, status
        elif mode == 2:
            mat_qp = mats[0][nm.newaxis,:,:].repeat( n_qp, 0 )
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
                status = self.function( out, -0.5, nm.empty( 0 ),
                                        apc.bf['v'], mat_qp, vgr, chunk, 1 )
                yield out, chunk, status
        else:
            cache1 = self.get_cache( 'state_in_volume_qp', 0 )
            cache2 = self.get_cache( 'state_in_volume_qp', 1 )
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,:].repeat( n_qp, 0 )
                    p1 = cache1( 'state', self.char_fun.ig, ii,
                                 state = state1, history = history )
                    p2 = cache2( 'state', self.char_fun.ig, ii,
                                 state = state2, history = history )
                    p10 = cache1( 'state', self.char_fun.ig, ii + 1,
                                  state = state1, history = history )
                    p20 = cache2( 'state', self.char_fun.ig, ii + 1,
                                  state = state2, history = history )
                    adp12 = (p1 - p2) + (p10 - p20)
                    status = self.function( out1, 0.5, adp12,
                                            apc.bf['v'], mat_qp, vgr, chunk, 0 )
                    out += out1
                yield out, chunk, status

class HDPMConvRRTerm( Term ):
    name = 'dw_hdpm_bh_r'
    arg_types = ('ts', 'material', 'virtual', 'parameter_1', 'parameter_2')
    geometry = [(Volume, 'virtual')]
    use_caches = {'state_in_volume_qp' : [['parameter_1'], ['parameter_1']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_hdpm_b1 )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, bstate1, bstate2\
            = self.get_args( **kwargs )
        apr, vgr = virtual.get_current_approximation()
        apc, vgc = bstate1.get_current_approximation()
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, dim * n_epr, 1)
            mode = 0
        else:
            raise StopIteration

        if ts.step == 0:
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                yield out, chunk, 0
            return

        cache1 = self.get_cache( 'state_in_volume_qp', 0 )
        cache2 = self.get_cache( 'state_in_volume_qp', 1 )
        p1 = cache1( 'state', self.char_fun.ig, 0, state = bstate1 )
        p2 = cache2( 'state', self.char_fun.ig, 0, state = bstate2 )
        dvec = p1 - p2

        for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                zero = True ):
            out1 = nm.empty_like( out )
            for ii, mat in enumerate( mats ):
                mat_qp = mat[nm.newaxis,:,:].repeat( n_qp, 0 )
##                    print ii, -ii-1
                status = self.function( out1, 1.0, dvec, apc.bf['v'],
                                        mat_qp, vgr, chunk, 0 )
                out += out1
            yield out, chunk, status

class HDPMConvQTerm( Term ):
    name = 'dw_hdpm_bbh'
    arg_types = ('ts', 'material', 'virtual', 'state', 'history')
    geometry = [(Volume, 'state')]
    use_caches = {'cauchy_strain' : [['state', {'strain' : (2,2),
                                               'dstrain' : (-1,-1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_hdpm_b2 )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state, history = self.get_args( **kwargs )
        apr, vgr = virtual.get_current_approximation()
        apc, vgc = state.get_current_approximation()
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, n_epr, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape()[3]
            shape = (chunk_size, 1, n_epr, dim * n_epc)
            mode = 1
        else:
            raise StopIteration
        
        if ts.step == 0:
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                yield out, chunk, 0
            return

        if mode == 1:
            mat_qp = mats[0][nm.newaxis,:,:].repeat( n_qp, 0 )
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
                status = self.function( out, 1.0 / ts.dt, nm.empty( 0 ),
                                        apr.bf['v'], mat_qp, vgc, chunk, 1 )
                yield out, chunk, status
        else:
            cache = self.get_cache( 'cauchy_strain', 0 )
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,:].repeat( n_qp, 0 )
                    dstrain = cache( 'dstrain', self.char_fun.ig, ii,
                                     state = state, history = history )
                    status = self.function( out1, 1.0 / ts.dt, dstrain,
                                            apr.bf['v'], mat_qp, vgc, chunk, 0 )
                    out += out1
                yield out, chunk, status

class HDPMConvGDPTerm( Term ):
    name = 'dw_hdpm_gdph'
    geometry = [(Volume, 'virtual')]
    arg_types = ('ts', 'material', 'virtual', 'state_1', 'state_2', 'history')
    use_caches = {'state_in_volume_qp' : [['state_1', {'state' : (-1,-1)}],
                                         ['state_2', {'state' : (-1,-1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_hdpm_g )


    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state1, state2, history\
            = self.get_args( **kwargs )
        ap, vg = virtual.get_current_approximation()
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state_1' ):
            shape = (chunk_size, 1, n_ep, n_ep)
            mode = 1
        elif diff_var == self.get_arg_name( 'state_2' ):
            shape = (chunk_size, 1, n_ep, n_ep)
            mode = 2
        else:
            raise StopIteration
        
        if ts.step == 0:
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                yield out, chunk, 0
            return

        if (mode == 1) or (mode == 2):
            mat_qp = mats[0][nm.newaxis,:,:].repeat( n_qp, 0 )
##             print mat_qp
##             pause( 'aa')
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape ):
                status = self.function( out, 1.0 / ts.dt, nm.empty( 0 ),
                                        ap.bf['v'], mat_qp, vg, chunk, mode )
                yield out, chunk, status
        else:
            cache1 = self.get_cache( 'state_in_volume_qp', 0 )
            cache2 = self.get_cache( 'state_in_volume_qp', 1 )
            for out, chunk in vector_chunk_generator( n_el, chunk_size, shape,
                                                    zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,:].repeat( n_qp, 0 )
                    p1 = cache1( 'state', self.char_fun.ig, ii,
                                 state = state1, history = history )
                    p2 = cache2( 'state', self.char_fun.ig, ii,
                                 state = state2, history = history )
                    p10 = cache1( 'state', self.char_fun.ig, ii + 1,
                                  state = state1, history = history )
                    p20 = cache2( 'state', self.char_fun.ig, ii + 1,
                                  state = state2, history = history )
                    ddp12 = (p1 - p2) - (p10 - p20)
                    status = self.function( out1, 1.0 / ts.dt, ddp12,
                                            ap.bf['v'], mat_qp, vg, chunk, 0 )
                    out += out1
##                     print chunk
##                     print out1
##                     pause()
                yield out, chunk, status
