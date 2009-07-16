from sfepy.terms.terms import *
from sfepy.terms.terms_base import ScalarScalar
from sfepy.terms.utils import fix_scalar_constant, fix_scalar_in_el

class LaplaceTerm( ScalarScalar, Term ):
    r""":description: Laplace term with $c$ constant or constant per element.
    :definition: $c \int_{\Omega}\nabla s \cdot \nabla r$
    or $\sum_{K \in \Tcal_h}\int_{T_K} c_K\ \nabla s \cdot \nabla r$
    """
    name = 'dw_laplace'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    symbolic = {'expression': 'c * div( grad( u ) )',
                'map' : {'u' : 'state', 'c' : 'material'}}

    def __init__(self, region, name=name, sign=1):
        Term.__init__(self, region, name, sign, terms.dw_laplace)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )

        if state.is_real():
            fargs = vec, 0, mat, vg, ap.econn
        else:
            ac = nm.ascontiguousarray
            fargs = [(ac( vec.real ), 0, mat, vg, ap.econn),
                     (ac( vec.imag ), 0, mat, vg, ap.econn)]
            mode += 1j
            
        return fargs, shape, mode

class DiffusionTerm( ScalarScalar, Term ):
    r""":description: General diffusion term with permeability $K_{ij}$
    constant or  given in mesh vertices. Can be evaluated. Can use derivatives.

    :definition: $\int_{\Omega} K_{ij} \nabla_i q \nabla_j p$, $\int_{\Omega}
    K_{ij} \nabla_i \bar{p} \nabla_j r$
    """
    name = 'dw_diffusion'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter_1'), (Volume, 'parameter_2')])
    modes = ('weak', 'eval')
    symbolic = {'expression': 'div( K * grad( u ) )',
                'map' : {'u' : 'state', 'K' : 'material'}}

    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )

        n_el, n_qp, dim, n_ep = self.data_shape
        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                        mat = nm.asarray( mat ), ap = ap,
                        assumed_shapes = [(1, n_qp, dim, dim),
                                          (n_el, n_qp, dim, dim)],
                        mode_in = None )
        return (1.0, vec, 0, mat_qp, vg, ap.econn), shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par1, par2 = self.get_args( **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        n_el, n_qp, dim, n_ep = self.data_shape

        cache = self.get_cache( 'grad_scalar', 0 )
        gp1 = cache( 'grad', self.get_current_group(), 0,
                     state = par1, get_vector = self.get_vector )
        cache = self.get_cache( 'grad_scalar', 1 )
        gp2 = cache( 'grad', self.get_current_group(), 0,
                     state = par2, get_vector = self.get_vector )

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                        mat = mat, ap = ap,
                        assumed_shapes = [(1, n_qp, dim, dim),
                                          (n_el, n_qp, dim, dim)],
                        mode_in = None )

        return (1.0, gp1, gp2, mat_qp, vg), (chunk_size, 1, 1, 1), 0

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = terms.dw_diffusion
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
            self.use_caches = {'mat_in_qp' : [['material']]}
        else:
            self.function = terms.d_diffusion
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'grad_scalar' : [['parameter_1'],
                                                ['parameter_2']],
                               'mat_in_qp' : [['material']]}

class PermeabilityRTerm( Term ):
    r""":description: Special-purpose diffusion-like term with permeability
    $K_{ij}$ constant or given in mesh vertices (to use on a right-hand side).
    :definition: $\int_{\Omega} K_{ij} \nabla_j q$
    """
    name = 'dw_permeability_r'
    arg_types = ('material', 'virtual', 'index')
    geometry = [(Volume, 'virtual')]
    use_caches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_permeability_r )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, index = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
        else:
            raise StopIteration

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                        mat = mat, ap = ap,
                        assumed_shapes = [(1, n_qp, dim, dim),
                                        (n_el, n_qp, dim, dim)],
                        mode_in = None )
        mat_qp = nm.ascontiguousarray( mat_qp[...,index:index+1] )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, mat_qp, vg, ap.econn, chunk )
            yield out, chunk, status

class DiffusionVelocityTerm( Term ):
    r""":description: Diffusion velocity averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} -K_{ij} \nabla_j r
    / \int_{T_K} 1$
    """
    name = 'de_diffusion_velocity'
    arg_types = ('material','parameter')
    geometry = [(Volume, 'parameter')]
    use_caches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.de_diffusion_velocity )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, parameter = self.get_args( **kwargs )
        ap, vg = parameter.get_approximation( self.get_current_group(),
                                              'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim, 1)
        else:
            raise StopIteration

        vec = parameter()
        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                        mat = nm.asarray( mat ), ap = ap,
                        assumed_shapes = [(1, n_qp, dim, dim),
                                          (n_el, n_qp, dim, dim)],
                        mode_in = None )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0,
                                    mat_qp, vg, ap.econn, chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status
