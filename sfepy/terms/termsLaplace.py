from terms import *
from utils import fix_scalar_constant, fix_scalar_in_el

##
# 28.11.2005, c
class LaplaceTerm( Term ):
    r""":description: Laplace term with $c$ constant or constant per element.
    :definition: $c \int_{\Omega}\nabla s \cdot \nabla r$
    or $\sum_{K \in \Tcal_h}\int_{T_K} c_K\ \nabla s \cdot \nabla r$
    """
    name = 'dw_laplace'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    symbolic = {'expression': 'c * div( grad( u ) )',
                'map' : {'u' : 'state', 'c' : 'material'}}

    ##
    # c: 28.11.2005, r: 27.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    ##
    # c: 27.03.2008, r: 27.03.2008
    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name )
        n_el, n_qp, dim, n_ep = self.data_shape
        if diff_var is None:
            return (chunk_size, 1, n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, n_ep, n_ep), 1
        else:
            raise StopIteration

    ##
    # c: 27.03.2008, r: 28.03.2008
    def build_c_fun_args( self, mat, state, ap, vg ):
        vec, indx = state()
        mat_arg = fix_scalar_constant( mat, nm.float64 )
        if mat_arg is None:
            mat_arg = fix_scalar_in_el( mat, self.data_shape[0], nm.float64 )
            self.function = terms.dw_st_pspg_p
        else:
            self.function = terms.dw_laplace
        fargs = vec, indx.start, mat_arg, vg, ap.econn

        return fargs
        
    ##
    # c: 28.11.2005, r: 28.03.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        material, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, ap )
        fargs = self.build_c_fun_args( material, state, ap, vg )
        
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

##
# 01.08.2006, c
class DiffusionTerm( LaplaceTerm ):
    r""":description: General diffusion term with permeability $K_{ij}$
    constant or given in mesh vertices.
    :definition: $\int_{\Omega} K_{ij} \nabla_i q  \nabla_j p$
    """
    name = 'dw_diffusion'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'mat_in_qp' : [['material']]}
    symbolic = {'expression': 'div( K * grad( u ) )',
                'map' : {'u' : 'state', 'K' : 'material'}}

    ##
    # c: 03.08.2006?, r: 28.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_diffusion )
        
    ##
    # c: 28.03.2008, r: 28.03.2008
    def build_c_fun_args( self, mat, state, ap, vg ):
        vec, indx = state()

        n_el, n_qp, dim, n_ep = self.data_shape
        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = nm.asarray( mat ), ap = ap,
                       assumed_shapes = [(1, n_qp, dim, dim),
                                        (n_el, n_qp, dim, dim)],
                       mode_in = None )
        return 1.0, vec, indx.start, mat_qp, vg, ap.econn

##
# 12.03.2007, c
class DiffusionIntegratedTerm( Term ):
    r""":description: Integrated general diffusion term with permeability
    $K_{ij}$  constant or given in mesh vertices.
    :definition: $\int_{\Omega} K_{ij} \nabla_i \bar{p}  \nabla_j r$
    """
    name = 'd_diffusion'
    arg_types = ('material', 'parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    use_caches = {'grad_scalar' : [['parameter_1'], ['parameter_2']],
                 'mat_in_qp' : [['material']]}

    ##
    # c: 12.03.2007, r: 28.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_diffusion )
        
    ##
    # c: 12.03.2007, r: 28.03.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par1, par2 = self.get_args( **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )
        shape = (chunk_size, 1, 1, 1)

        cache = self.get_cache( 'grad_scalar', 0 )
        gp1 = cache( 'grad', self.get_current_group(), 0, state = par1 )
        cache = self.get_cache( 'grad_scalar', 1 )
        gp2 = cache( 'grad', self.get_current_group(), 0, state = par2 )

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = mat, ap = ap,
                       assumed_shapes = [(1, n_qp, dim, dim),
                                        (n_el, n_qp, dim, dim)],
                       mode_in = None )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, 1.0, gp1, gp2, mat_qp, vg, chunk )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

##
# 23.04.2007, c
class PermeabilityRTerm( Term ):
    r""":description: Special-purpose diffusion-like term with permeability
    $K_{ij}$ constant or given in mesh vertices (to use on a right-hand side).
    :definition: $\int_{\Omega} K_{ij} \nabla_j q$
    """
    name = 'dw_permeability_r'
    arg_types = ('material', 'virtual', 'index')
    geometry = [(Volume, 'virtual')]
    use_caches = {'mat_in_qp' : [['material']]}

    ##
    # c: 23.04.2007, r: 28.03.2008
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_permeability_r )
        
    ##
    # c: 23.04.2007, r: 28.03.2008
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

##
# 07.09.2006, c
class DiffusionVelocityTerm( Term ):
    r""":description: Diffusion velocity averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} K_{ij} \nabla_j r
    / \int_{T_K} 1$
    """
    name = 'de_diffusion_velocity'
    arg_types = ('material','parameter')
    geometry = [(Volume, 'parameter')]
    use_caches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.de_diffusion_velocity )
        
    ##
    # c: 07.09.2006, r: 14.07.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, parameter = self.get_args( **kwargs )
        ap, vg = parameter.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim, 1)
        else:
            raise StopIteration

        vec, indx = parameter()
        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = nm.asarray( mat ), ap = ap,
                       assumed_shapes = [(1, n_qp, dim, dim),
                                        (n_el, n_qp, dim, dim)],
                       mode_in = None )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, indx.start,
                                    mat_qp, vg, ap.econn, chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status
