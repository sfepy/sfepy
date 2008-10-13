from terms import *
from terms_base import CouplingVectorScalar
from utils import fix_mat_shape

class PiezoCouplingGrad( CouplingVectorScalar ):

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
            cache = self.get_cache( 'grad_scalar', 0 )
            p_grad = cache( 'grad', self.get_current_group(), 0, state = state )
        else:
            p_grad  = aux

        n_qp = self.data_shape_r[1]
        mat = nm.asarray( mat, dtype = nm.float64 )
        mat = fix_mat_shape( mat, n_qp )

        return (aux, p_grad, mat, vgr), shape, mode

class PiezoCouplingDiv( CouplingVectorScalar ):

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
                            state = state )
        else:
            strain = aux

        n_qp = self.data_shape_r[1]
        mat = nm.asarray( mat, dtype = nm.float64 )
        mat = fix_mat_shape( mat, n_qp )

        return (strain, aux, mat, vgc), shape, mode + 2

class PiezoCouplingTerm( PiezoCouplingDiv, PiezoCouplingGrad, Term ):
    r""":description: Piezoelectric coupling term.
    :definition: $\int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q$,
    $\int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p$
    """
    name = 'dw_piezo_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'))
    geometry = ([(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'virtual'), (Volume, 'state')])

    def set_arg_types( self ):
        """Dynamically inherits from either PiezoCouplingGrad or
        PiezoCouplingDiv."""
        if self.ats[1] == 'virtual':
            self.mode = 'grad'
            self.function = terms.dw_piezo_coupling
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
            self.use_caches = {'grad_scalar' : [['state']]}
        elif self.ats[2] == 'virtual':
            self.mode = 'div'
            self.function = terms.dw_piezo_coupling
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['state']]}
        else:
            self.mode = 'eval'
            raise NotImplementedError
