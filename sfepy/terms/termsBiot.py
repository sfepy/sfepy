from terms import *
from terms_base import CouplingVectorScalar
from termsNavierStokes import GradTerm, DivTerm

class BiotGrad( CouplingVectorScalar ):

    def get_fargs_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
        if diff_var is None:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            vec_qp = cache( 'state', self.get_current_group(), 0, state = state )
        else:
            vec_qp = aux

        n_qp = self.data_shape_r[1]
        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
#        print mat_qp
        bf = apc.get_base( 'v', 0, self.integral_name )

        return (1.0, vec_qp, bf, mat_qp, vgr), shape, mode

class BiotDiv( CouplingVectorScalar ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, state, virtual = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
        if diff_var is None:
            cache = self.get_cache( 'cauchy_strain', 0 )
            strain = cache( 'strain', self.get_current_group(), 0,
                            state = state )
        else:
            strain = aux

        n_qp = self.data_shape_r[1]
        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
        bf = apr.get_base( 'v', 0, self.integral_name )

        return (1.0, strain, bf, mat_qp, vgc), shape, mode

class BiotEval( Struct ):

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        raise NotImplementedError

class BiotTerm( BiotDiv, BiotGrad, BiotEval, Term ):
    r""":description: Biot coupling term with $\alpha_{ij}$
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$. Corresponds to weak
    forms of Biot gradient and divergence terms. Can be evaluated.
    :definition: $\int_{\Omega}  p\ \alpha_{ij} e_{ij}(\ul{v})$, $\int_{\Omega}
    q\ \alpha_{ij} e_{ij}(\ul{u})$
    """
    name = 'dw_biot'
    arg_types = ('material', 'virtual|state', 'state|virtual')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    def set_arg_types( self ):
        """Dynamically inherits from either BiotGrad or
        BiotDiv."""
        if self.ats[1] == 'virtual':
            self.mode = 'grad'
            self.function = terms.dw_biot_grad
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['state']]}
        elif self.ats[2] == 'virtual':
            self.mode = 'div'
            self.function = terms.dw_biot_div
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['state']]}
        else:
            self.mode = 'eval'
            use_method_with_name( self, self.call_eval, '_call' )
            raise NotImplementedError

##
# 01.08.2006, c
class BiotGradTerm( GradTerm ):
    r""":description: Biot gradient-like term (weak form) with $\alpha_{ij}$
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega}  p\ \alpha_{ij} e_{ij}(\ul{v})$
    """
    name = 'dw_biot_grad'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'state_in_volume_qp' : [['state']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_grad )

    ##
    # c: 31.03.2008, r: 31.03.2008
    def build_c_fun_args( self, state, apc, vgr, **kwargs ):
        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.get_current_group(), 0, state = state )

        mat, = self.get_args( ['material'], **kwargs )
        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( self.data_shape[1], 0 )
#        print mat_qp
        bf = apc.get_base( 'v', 0, self.integral_name )
        return 1.0, vec_qp, bf, mat_qp, vgr

##
# c: 18.05.2008
class BiotGradDtTerm( GradTerm ):
    r""":description: Biot gradient-like term (weak form) with time-discretized
    $\dot{p}$ and $\alpha_{ij}$
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega}  \frac{p - p_0}{\dt}\ \alpha_{ij} e_{ij}(\ul{v})$
    :arguments: ts.dt : $\dt$, parameter : $p_0$
    """
    name = 'dw_biot_grad_dt'
    arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'state_in_volume_qp' : [['state', {'state' : (2, 2)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_grad )

    ##
    # c: 18.05.2008, r: 18.05.2008
    def build_c_fun_args( self, state, apc, vgr, **kwargs ):
        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.get_current_group(), 0, state = state )
        vec0_qp = cache( 'state', self.get_current_group(), 1,
                         state = parameter )

        ts, mat = self.get_args( ['ts', 'material'], **kwargs )
        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( self.data_shape[1], 0 )
#        print mat_qp
        bf = apc.get_base( 'v', 0, self.integral_name )
        return 1.0 / ts.dt, vec_qp - vec0_qp, bf, mat_qp, vgr

##
# c: 03.04.2008
class BiotGradTHTerm( BiotGradTerm ):
    r""":definition: $\int_{\Omega} \left [\int_0^t
    \alpha_{ij}(t-\tau)\,p(\tau)) \difd{\tau} \right]\,e_{ij}(\ul{v})$
    """
    name = 'dw_biot_grad_th'
    arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'state_in_volume_qp' : [['state', {'state' : (-1,-1)}]]}

    ##
    # c: 03.04.2008, r: 18.06.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """history for now is just state_0, it is not used anyway, as the
        history is held in the dstrain cache"""
        ts, mats, virtual, state, history = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, apr, apc )
        n_el, n_qp, dim, n_ep = self.data_shape

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        bf = apc.get_base( 'v', 0, self.integral_name )

        if mode == 1:
            mat_qp = mats[0][nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = self.function( out, ts.dt, nm.empty( 0 ), bf,
                                        mat_qp, vgr, chunk, 1 )
                yield out, chunk, status
        else:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            for out, chunk in self.char_fun( chunk_size, shape, zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
                    vec_qp = cache( 'state', self.get_current_group(), ii,
                                    state = state, history = history )
                    status = self.function( out1, ts.dt, vec_qp, bf,
                                            mat_qp, vgr, chunk, 0 )
                    out += out1
                yield out, chunk, status

##
# 01.08.2006, c
class BiotDivDtTerm( DivTerm ):
    r""":description: Biot divergence-like rate term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega} q\ \alpha_{ij} \frac{e_{ij}(\ul{u}) -
    e_{ij}(\ul{u_0})}{\dt}$
    """
    name = 'dw_biot_div_dt'
    arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'cauchy_strain' : [['state', {'strain' : (2,2),
                                               'dstrain' : (1,1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_div )

    ##
    # c: 31.03.2008, r: 02.04.2008
    def build_c_fun_args( self, state, apr, apc, vgc, **kwargs ):
        cache = self.get_cache( 'cauchy_strain', 0 )
        dstrain = cache( 'dstrain', self.get_current_group(), 0, state = state )

        ts, mat = self.get_args( ['ts', 'material'], **kwargs )
        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( self.data_shape[1], 0 )
#        print mat_qp
        bf = apr.get_base( 'v', 0, self.integral_name )
        return 1.0 / ts.dt, dstrain, bf, mat_qp, vgc


##
# 20.09.2006, c
class BiotDivTerm( DivTerm ):
    r""":description: Biot divergence-like term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega} q\ \alpha_{ij} e_{ij}(\ul{u})$
    """
    name = 'dw_biot_div'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'cauchy_strain' : [['state']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_biot_div )

    ##
    # c: 31.03.2008, r: 31.03.2008
    def build_c_fun_args( self, state, apr, apc, vgc, **kwargs ):
        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.get_current_group(), 0, state = state )
        mat, = self.get_args( ['material'], **kwargs )
        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( self.data_shape[1], 0 )
        bf = apr.get_base( 'v', 0, self.integral_name )
        return 1.0, strain, bf, mat_qp, vgc

##
# c: 05.03.2008
class BiotDivRIntegratedTerm( Term ):
    r""":description: Integrated Biot divergence-like term (weak form) with
    $\alpha_{ij}$ given in vector form exploiting symmetry: in 3D it has the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it has
    the indices ordered as $[11, 22, 12]$.
    :definition: $\int_{\Omega} r\ \alpha_{ij} e_{ij}(\ul{w})$
    """
    name = 'd_biot_div'
    arg_types = ('material', 'parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    use_caches = {'state_in_volume_qp' : [['parameter_1']],
                 'cauchy_strain' : [['parameter_2']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_biot_div )

    ##
    # c: 05.03.2008, r: 05.03.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par_p, par_u = self.get_args( **kwargs )
        apr, vgr = par_p.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = par_u.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape( self.integral_name )
        
        shape = (chunk_size, 1, 1, 1)
        
        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.get_current_group(), 0, state = par_p )

        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.get_current_group(), 0, state = par_u )

        mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, 1.0, vec_qp, strain, mat_qp, vgc, chunk )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

##
# c: 03.04.2008
class BiotDivTHTerm( BiotDivTerm ):
    r""":definition: $\int_{\Omega} \left [\int_0^t
    \alpha_{ij}(t-\tau) \tdiff{e_{kl}(\ul{u}(\tau))}{\tau}
    \difd{\tau} \right] q$
    """
    name = 'dw_biot_div_th'
    arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'cauchy_strain' : [['state', {'strain' : (2,2),
                                               'dstrain' : (-1,-1)}]]}

    ##
    # c: 03.04.2008, r: 18.06.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """history for now is just state_0, it is not used anyway, as the
        history is held in the dstrain cache"""
        ts, mats, virtual, state, history = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(), 'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, apr, apc )
        n_el, n_qp, dim, n_ep = self.data_shape

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        bf = apr.get_base( 'v', 0, self.integral_name )

        if mode == 1:
            mat_qp = mats[0][nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = self.function( out, 1.0, nm.empty( 0 ), bf,
                                        mat_qp, vgc, chunk, 1 )
                yield out, chunk, status
        else:
            cache = self.get_cache( 'cauchy_strain', 0 )
            for out, chunk in self.char_fun( chunk_size, shape, zero = True ):
                out1 = nm.empty_like( out )
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
                    dstrain = cache( 'dstrain', self.get_current_group(), ii,
                                     state = state, history = history )
                    status = self.function( out1, 1.0, dstrain,
                                            bf, mat_qp, vgc, chunk, 0 )
                    out += out1
                yield out, chunk, status
