from sfepy.terms.terms import *
from sfepy.terms.terms_base import CouplingVectorScalar, CouplingVectorScalarTH

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
            vec_qp = cache( 'state', self.get_current_group(), 0,
                            state = state, get_vector = self.get_vector )
        else:
            vec_qp = aux

        bf = apc.get_base( 'v', 0, self.integral_name )

        return (1.0, vec_qp, bf, mat, vgr), shape, mode

class BiotDiv( CouplingVectorScalar ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, state, virtual = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_div( diff_var, chunk_size )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
        if diff_var is None:
            cache = self.get_cache( 'cauchy_strain', 0 )
            strain = cache( 'strain', self.get_current_group(), 0,
                            state = state, get_vector = self.get_vector )
        else:
            strain = aux

        bf = apr.get_base( 'v', 0, self.integral_name )

        return (1.0, strain, bf, mat, vgc), shape, mode

class BiotEval( CouplingVectorScalar ):

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par_v, par_s = self.get_args( **kwargs )
        aps, vgs = par_s.get_approximation( self.get_current_group(), 'Volume' )
        apv, vgv = par_v.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( aps, apv )
        return (mat, par_v, par_s, vgv), (chunk_size, 1, 1, 1), 0

    def d_eval( self, out, mat, par_v, par_s, vgv, chunk ):
        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.get_current_group(), 0,
                        state = par_s, get_vector = self.get_vector )

        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.get_current_group(), 0,
                        state = par_v, get_vector = self.get_vector )

        function = terms.d_biot_div
        status = function( out, 1.0, vec_qp, strain, mat, vgv, chunk )
        return status

class BiotTerm( BiotGrad, BiotDiv, BiotEval, Term ):
    r"""
    :Description:
    Biot coupling term with :math:`\alpha_{ij}`
    given in vector form exploiting symmetry: in 3D it has the
    indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has
    the indices ordered as :math:`[11, 22, 12]`. Corresponds to weak
    forms of Biot gradient and divergence terms. Can be evaluated. Can
    use derivatives.
    
    :Definition:
    .. math::
        \int_{\Omega}  p\ \alpha_{ij} e_{ij}(\ul{v}) \mbox{ , } \int_{\Omega}
        q\ \alpha_{ij} e_{ij}(\ul{u})
    """
    name = 'dw_biot'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    geometry = ([(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'parameter_v'), (Volume, 'parameter_s')])
    modes = ('grad', 'div', 'eval')

    def set_arg_types( self ):
        """Dynamically inherits from either BiotGrad, BiotDiv or BiotEval."""
        if self.mode == 'grad':
            self.function = terms.dw_biot_grad
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['state']]}
        elif self.mode == 'div':
            self.function = terms.dw_biot_div
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['state']]}
        else:
            self.function = self.d_eval
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['parameter_s']],
                               'cauchy_strain' : [['parameter_v']]}


class BiotGradTH( CouplingVectorScalarTH ):

    def get_fargs_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        bf = apc.get_base( 'v', 0, self.integral_name )
        n_el, n_qp = self.data_shape_r[:2]

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat = mats[0]
            mat = nm.tile(mat, (n_el, n_qp, 1, 1))
            return (ts.dt, aux, bf, mat, vgr), shape, mode

        else:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    vec_qp = cache( 'state', self.get_current_group(), ii,
                                    state = state, get_vector = self.get_vector )
                    mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                    yield ii, (ts.dt, vec_qp, bf, mat, vgr)
            return iter_kernel, shape, mode

class BiotDivTH( CouplingVectorScalarTH ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, state, virtual = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_div( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        bf = apr.get_base( 'v', 0, self.integral_name )
        n_el, n_qp = self.data_shape_r[:2]

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat = mats[0]
            mat = nm.tile(mat, (n_el, n_qp, 1, 1))
            return (ts.dt, aux, bf, mat, vgc), shape, mode

        else:
            cache = self.get_cache( 'cauchy_strain', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    strain = cache( 'strain', self.get_current_group(), ii,
                                    state = state, get_vector = self.get_vector )
                    mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                    yield ii, (ts.dt, strain, bf, mat, vgc)
            return iter_kernel, shape, mode

class BiotTHTerm( BiotGradTH, BiotDivTH, Term ):
    r"""
    :Description:
    Fading memory Biot term. Can use derivatives.

    :Definition:
    .. math::
        \begin{array}{l}
        \int_{\Omega} \left [\int_0^t \alpha_{ij}(t-\tau)\,p(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v}) \mbox{ ,} \\
        \int_{\Omega} \left [\int_0^t
        \alpha_{ij}(t-\tau) e_{kl}(\ul{u}(\tau)) \difd{\tau} \right] q
        \end{array}
    """
    name = 'dw_biot_th'
    arg_types = (('ts', 'material', 'virtual', 'state'),
                 ('ts', 'material', 'state', 'virtual'))
    geometry = ([(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'virtual'), (Volume, 'state')])
    modes = ('grad', 'div')

    def set_arg_types( self ):
        """Dynamically inherits from either BiotGradTH or
        BiotDivTH."""
        if self.mode == 'grad':
            self.function = terms.dw_biot_grad
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['state',
                                                        {'state' : (-1,-1)}]]}
        elif self.mode == 'div':
            self.function = terms.dw_biot_div
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['state',
                                                   {'strain' : (-1,-1)}]]}
        else:
            raise NotImplementedError

class BiotGradETH( CouplingVectorScalar ):

    def get_fargs_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat0, mat1, virtual, state = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        bf = apc.get_base( 'v', 0, self.integral_name )
        if diff_var is None:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            vec_qp = cache( 'state', self.get_current_group(), 0,
                            state = state, get_vector = self.get_vector )

            cache = self.get_cache('exp_history', 0)
            increment = cache('increment', self.get_current_group(), 0,
                              decay=mat1, values=vec_qp)
            history = cache('history', self.get_current_group(), 0)

            fargs = (ts.dt, history + increment, bf, mat0, vgr)
            if ts.step == 0: # Just init the history in step 0.
                raise StopIteration

        else:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            fargs = (ts.dt, aux, bf, mat0, vgr)

        return fargs, shape, mode

class BiotDivETH( CouplingVectorScalar ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat0, mat1, state, virtual = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_div( diff_var, chunk_size )

        bf = apr.get_base( 'v', 0, self.integral_name )
        if diff_var is None:
            cache = self.get_cache( 'cauchy_strain', 0 )
            strain = cache( 'strain', self.get_current_group(), 0,
                            state = state, get_vector = self.get_vector )

            cache = self.get_cache('exp_history', 0)
            increment = cache('increment', self.get_current_group(), 0,
                              decay=mat1, values=strain)
            history = cache('history', self.get_current_group(), 0)

            fargs = (ts.dt, history + increment, bf, mat0, vgc)
            if ts.step == 0: # Just init the history in step 0.
                raise StopIteration

        else:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            fargs = (ts.dt, aux, bf, mat0, vgc)

        return fargs, shape, mode

class BiotETHTerm( BiotGradETH, BiotDivETH, Term ):
    r"""
    :Description:
    This term has the same definition as dw_biot_th, but assumes an
    exponential approximation of the convolution kernel resulting in much
    higher efficiency. Can use derivatives.

    :Definition:
    .. math::
        \begin{array}{l}
        \int_{\Omega} \left [\int_0^t \alpha_{ij}(t-\tau)\,p(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v}) \mbox{ ,} \\
        \int_{\Omega} \left [\int_0^t
        \alpha_{ij}(t-\tau) e_{kl}(\ul{u}(\tau)) \difd{\tau} \right] q
        \end{array}
    
    :Arguments:
    material_0 : :math:`\alpha_{ij}(0)`,
    material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
    """
    name = 'dw_biot_eth'
    arg_types = (('ts', 'material_0', 'material_1', 'virtual', 'state'),
                 ('ts', 'material_0', 'material_1', 'state', 'virtual'))
    geometry = ([(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'virtual'), (Volume, 'state')])
    modes = ('grad', 'div')
    use_caches = {'exp_history' : [['material_0', 'material_1', 'state']]}

    def set_arg_types( self ):
        """Dynamically inherits from either BiotGradETH or
        BiotDivETH."""
        if self.mode == 'grad':
            self.function = terms.dw_biot_grad
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
            use_caches = {'state_in_volume_qp' : [['state']]}
        elif self.mode == 'div':
            self.function = terms.dw_biot_div
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
            use_caches = {'cauchy_strain' : [['state']]}
        else:
            raise NotImplementedError

        self.use_caches.update(use_caches)
