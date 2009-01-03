from sfepy.terms.terms import *
from sfepy.terms.terms_base import CouplingVectorScalar
from sfepy.terms.utils import fix_mat_shape

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
                            state = state, get_vector = self.get_vector )
        else:
            strain = aux

        n_qp = self.data_shape_r[1]
        mat = nm.asarray( mat, dtype = nm.float64 )
        mat = fix_mat_shape( mat, n_qp )

        return (strain, aux, mat, vgc), shape, mode + 2

class  PiezoCouplingEval( CouplingVectorScalar ):

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        if diff_var is not None:
            raise StopIteration
        
        mat, par_v, par_s = self.get_args( **kwargs )
        aps, vgs = par_s.get_approximation( self.get_current_group(),
                                            'Volume' )
        apv, vgv = par_v.get_approximation( self.get_current_group(),
                                            'Volume' )
        self.set_data_shape( aps, apv )

        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.get_current_group(), 0,
                        state = par_v, get_vector = self.get_vector )
        cache = self.get_cache( 'grad_scalar', 0 )
        p_grad = cache( 'grad', self.get_current_group(), 0, state = par_s )

        n_qp = self.data_shape_r[1]
        mat = nm.asarray( mat, dtype = nm.float64 )
        mat = fix_mat_shape( mat, n_qp )

        return (strain, p_grad, mat, vgv), (chunk_size, 1, 1, 1), 0

class PiezoCouplingTerm( PiezoCouplingDiv, PiezoCouplingGrad,
                         PiezoCouplingEval, Term ):
    r""":description: Piezoelectric coupling term.
    :definition: $\int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q$,
    $\int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p$
    """
    name = 'dw_piezo_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    geometry = ([(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'virtual'), (Volume, 'state')],
                [(Volume, 'parameter_v'), (Volume, 'parameter_s')])
    modes = ('grad', 'div', 'eval')

    def set_arg_types( self ):
        """Dynamically inherits from either PiezoCouplingGrad or
        PiezoCouplingDiv."""
        if self.mode == 'grad':
            self.function = terms.dw_piezo_coupling
            use_method_with_name( self, self.get_fargs_grad, 'get_fargs' )
            self.use_caches = {'grad_scalar' : [['state']]}
        elif self.mode == 'div':
            self.function = terms.dw_piezo_coupling
            use_method_with_name( self, self.get_fargs_div, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['state']]}
        else:
            self.function = terms.d_piezo_coupling
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'grad_scalar' : [['parameter_s']],
                               'cauchy_strain' : [['parameter_v']]}
