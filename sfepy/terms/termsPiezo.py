from sfepy.terms.terms import *
from sfepy.terms.terms_base import CouplingVectorScalar

class PiezoCouplingGrad( CouplingVectorScalar ):

    def get_fargs_grad( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )

        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(state)

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_grad( diff_var, chunk_size )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
        if diff_var is None:
            cache = self.get_cache( 'grad_scalar', 0 )
            p_grad = cache('grad', self, 0, state=state)
        else:
            p_grad  = aux

        return (aux, p_grad, mat, vgr), shape, mode

class PiezoCouplingDiv( CouplingVectorScalar ):

    def get_fargs_div( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, state, virtual = self.get_args( **kwargs )

        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(state)

        self.set_data_shape( apr, apc )
        shape, mode = self.get_shape_div( diff_var, chunk_size )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
        if diff_var is None:
            cache = self.get_cache( 'cauchy_strain', 0 )
            strain = cache('strain', self, 0,
                           state=state, get_vector=self.get_vector)
        else:
            strain = aux

        return (strain, aux, mat, vgc), shape, mode + 2

class  PiezoCouplingEval( CouplingVectorScalar ):

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        if diff_var is not None:
            raise StopIteration

        mat, par_v, par_s = self.get_args( **kwargs )

        aps, vgs = self.get_approximation(par_s)
        apv, vgv = self.get_approximation(par_v)

        self.set_data_shape( aps, apv )

        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache('strain', self, 0,
                       state=par_v, get_vector=self.get_vector)
        cache = self.get_cache( 'grad_scalar', 0 )
        p_grad = cache('grad', self, 0, state=par_s)

        return (strain, p_grad, mat, vgv), (chunk_size, 1, 1, 1), 0

class PiezoCouplingTerm( PiezoCouplingDiv, PiezoCouplingGrad,
                         PiezoCouplingEval, Term ):
    r"""
    :Description:
    Piezoelectric coupling term. Can be evaluated.

    :Definition:
    .. math::
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p \mbox{ , }
        \int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q

    :Arguments 1:
        material : :math:`g_{kij}`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`p`

    :Arguments 2:
        material : :math:`g_{kij}`,
        state    : :math:`\ul{u}`,
        virtual  : :math:`q`

    :Arguments 3:
        material    : :math:`g_{kij}`,
        parameter_v : :math:`\ul{u}`,
        parameter_s : :math:`p`
    """
    name = 'dw_piezo_coupling'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
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
