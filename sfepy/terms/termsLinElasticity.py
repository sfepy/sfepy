from sfepy.terms.terms import *
from sfepy.terms.terms_base import VectorVector, VectorVectorTH

## expr = """
## e = 1/2 * (grad( vec( u ) ) + grad( vec( u ) ).T)
## D = map( D_sym )
## s = D * e
## div( s )
## """

## """
## e[i,j] = 1/2 * (der[j]( u[i] ) + der[i]( u[j] ))
## map =
## D[i,j,k,l]
## s[i,j] = D[i,j,k,l] * e[k,l]
## """

class LinearElasticTerm( VectorVector, Term ):
    r""":description: General linear elasticity term, with $D_{ijkl}$ given in
    the usual matrix form exploiting symmetry: in 3D it is $6\times6$ with the
    indices ordered as $[11, 22, 33, 12, 13, 23]$, in 2D it is $3\times3$ with
    the indices ordered as $[11, 22, 12]$. Can be evaluated. Can use derivatives.

    :definition: $\int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})$
    """
    name = 'dw_lin_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter_1'), (Volume, 'parameter_2')])
    modes = ('weak', 'eval')
##     symbolic = {'expression': expr,
##                 'map' : {'u' : 'state', 'D_sym' : 'material'}}

    def check_mat_shape( self, mat ):
        dim = self.data_shape[2]
        sym = (dim + 1) * dim / 2
        assert_(mat.shape == (self.data_shape[0], self.data_shape[1], sym, sym))

    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
        if diff_var is None:
            cache = self.get_cache( 'cauchy_strain', 0 )
            strain = cache( 'strain', self.get_current_group(), 0,
                            state = state, get_vector = self.get_vector )
        else:
            strain = aux

        self.check_mat_shape( mat )

        return (1.0, strain, mat, vg), shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par1, par2 = self.get_args( **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )

        self.check_mat_shape( mat )

        cache = self.get_cache( 'cauchy_strain', 0 )
        strain1 = cache( 'strain', self.get_current_group(), 0,
                         state = par1, get_vector = self.get_vector )
        cache = self.get_cache( 'cauchy_strain', 1 )
        strain2 = cache( 'strain', self.get_current_group(), 0,
                         state = par2, get_vector = self.get_vector )

        return (1.0, strain1, strain2, mat, vg), (chunk_size, 1, 1, 1), 0

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = terms.dw_lin_elastic
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['state', {'strain' : (1,1)}]]}
        else:
            self.function = terms.d_lin_elastic
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'cauchy_strain' : [['parameter_1'],
                                                  ['parameter_2']]}

class LinearElasticIsotropicTerm( VectorVector, Term ):
    r""":description: Isotropic linear elasticity term.
    :definition: $\int_{\Omega}  D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})$
    with $D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}$ 
    """
    name = 'dw_lin_elastic_iso'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_elastic_iso )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )
        lam, mu = map( nm.float64, [mat[ii] for ii in ['lambda', 'mu']] )
        return (vec, 0, lam, mu, vg, ap.econn), shape, mode

class LinearElasticTHTerm( VectorVectorTH, Term ):
    r""":definition: $\int_{\Omega} \left [\int_0^t
    \Hcal_{ijkl}(t-\tau)\,\tdiff{e_{kl}(\ul{u}(\tau))}{\tau}
    \difd{\tau} \right]\,e_{ij}(\ul{v})$"""
    name = 'dw_lin_elastic_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'cauchy_strain' : [['state', {'strain' : (-1,-1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_lin_elastic )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        n_qp = self.data_shape[1]

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat_qp = mats[0][nm.newaxis,:,:].repeat( n_qp, 0 )
            return (ts.dt, aux, mat_qp, vg), shape, mode

        else:
            cache = self.get_cache( 'cauchy_strain', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,:].repeat( n_qp, 0 )
                    strain = cache( 'strain', self.get_current_group(), ii,
                                    state = state, get_vector = self.get_vector )
                    yield ii, (ts.dt, strain, mat_qp, vg)
            return iter_kernel, shape, mode

class CauchyStrainTerm( Term ):
    r""":description: Cauchy strain tensor averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} \ull{e}(\ul{w}) /
    \int_{T_K} 1$
    """
    name = 'de_cauchy_strain'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.de_cauchy_strain )

    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name )
        n_el, n_qp, dim, n_ep = self.data_shape
        
        if diff_var is None:
            return chunk_size, 1, dim * (dim + 1) / 2, 1
        else:
            raise StopIteration

    def build_c_fun_args( self, state, ap, vg, **kwargs ):
        vec = state()
        return vec, 0, vg, ap.econn
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        parameter, = self.get_args( ['parameter'], **kwargs )
        ap, vg = parameter.get_approximation( self.get_current_group(), 'Volume' )

        shape = self.get_shape( diff_var, chunk_size, ap )
        fargs = self.build_c_fun_args( parameter, ap, vg, **kwargs )
        
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk,) )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status

class CauchyStressTerm( CauchyStrainTerm ):
    r""":description: Cauchy stress tensor averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h:
    \int_{T_K} D_{ijkl} e_kl(\ul{w}) / \int_{T_K} 1$
    """
    name = 'de_cauchy_stress'
    arg_types = ('material', 'parameter')
    geometry = [(Volume, 'parameter')]
    use_caches = {'cauchy_strain' : [['parameter']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.de_cauchy_stress )

    def build_c_fun_args( self, state, ap, vg, **kwargs ):
        mat, = self.get_args( ['material'], **kwargs )
        mat_qp = mat[nm.newaxis,:,:].repeat( self.data_shape[1], 0 )
        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.get_current_group(), 0,
                        state = state, get_vector = self.get_vector )
        return strain, mat_qp, vg

