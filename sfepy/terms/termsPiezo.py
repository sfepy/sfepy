from terms import *
from utils import fix_mat_shape

class PiezoCouplingGrad( Struct ):

    def call_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )

        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
##         print diff_var, self.get_arg_name( 'state' ), self.mode, state.name
##        debug()

        data_shape = apr.get_v_data_shape( self.integral_name ) 
        n_el, n_qp, dim, n_epr = data_shape

        if diff_var is None:
            cache = self.get_cache( 'grad_scalar', 0 )
            p_grad = cache( 'grad', self.get_current_group(), 0, state = state )

            shape, mode = (chunk_size, 1, dim * n_epr, 1 ), 0

        elif diff_var == self.get_arg_name( 'state' ):
            p_grad  = aux
            
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            shape, mode = (chunk_size, 1, dim * n_epr, n_epc ), 1
        else:
            raise StopIteration

        mat = nm.asarray( mat, dtype = nm.float64 )
        mat = fix_mat_shape( mat, n_qp )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, aux, p_grad, mat, vgr, chunk, mode )
            yield out, chunk, status

class PiezoCouplingDiv( Struct ):

    def call_div( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, state, virtual = self.get_args( **kwargs )

        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
##         print diff_var, self.get_arg_name( 'state' ), self.mode, state.name
##        debug()

        data_shape = apr.get_v_data_shape( self.integral_name ) 
        n_el, n_qp, dim, n_epr = data_shape

        if diff_var is None:
            cache = self.get_cache( 'cauchy_strain', 0 )
            strain = cache( 'strain', self.get_current_group(), 0,
                            state = state )

            shape, mode =  (chunk_size, 1, n_epr, 1 ), 2

        elif diff_var == self.get_arg_name( 'state' ):
            strain = aux

            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            shape, mode = (chunk_size, 1, n_epr, dim * n_epc ), 3
        else:
            raise StopIteration

        mat = nm.asarray( mat, dtype = nm.float64 )
        mat = fix_mat_shape( mat, n_qp )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, strain, aux, mat, vgc, chunk, mode )
            yield out, chunk, status

class PiezoCouplingTerm( Term, PiezoCouplingDiv, PiezoCouplingGrad ):
    r""":description: Piezoelectric coupling term.
    :definition: $\int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q$,
    $\int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p$
    """
    name = 'dw_piezo_coupling'
    arg_types = ('material', 'virtual|state', 'state|virtual')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    def set_arg_types( self ):
        """Dynamically inherits from either PiezoCouplingGrad or
        PiezoCouplingDiv."""
        if self.ats[1] == 'virtual':
            self.mode = 'grad'
            self.function = terms.dw_piezo_coupling
            use_method_with_name( self, self.call_grad, '_call' )
            self.use_caches = {'grad_scalar' : [['state']]}
        elif self.ats[2] == 'virtual':
            self.mode = 'div'
            self.function = terms.dw_piezo_coupling
            use_method_with_name( self, self.call_div, '_call' )
            self.use_caches = {'cauchy_strain' : [['state']]}
        else:
            self.mode = 'eval'
            raise NotImplementedError
