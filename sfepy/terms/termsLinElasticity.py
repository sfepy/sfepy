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
    r"""
    :Description:
    General linear elasticity term, with :math:`D_{ijkl}` given in
    the usual matrix form exploiting symmetry: in 3D it is :math:`6\times6`
    with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it is
    :math:`3\times3` with the indices ordered as :math:`[11, 22, 12]`. Can be
    evaluated. Can use derivatives.

    :definition:
    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
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
    r"""
    :Description:
    Isotropic linear elasticity term.

    :Definition:
    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u}) \mbox{ with }
        D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
        \lambda \ \delta_{ij} \delta_{kl}

    :Arguments:
        material_1 : :math:`\lambda`,
        material_2 : :math:`\mu`,
        virtual :    :math:`\ul{v}`,
        state :      :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_iso'
    arg_types = ('material_1', 'material_2', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    function = staticmethod(terms.dw_lin_elastic_iso)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        lam, mu, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )

        return (vec, 0, lam, mu, vg, ap.econn), shape, mode

class LinearElasticTHTerm( VectorVectorTH, Term ):
    r"""
    :Description:
    Fading memory linear elastic (viscous) term. Can use derivatives.

    :Definition:
    .. math::
        \int_{\Omega} \left [\int_0^t
        \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{u}(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v})
    """
    name = 'dw_lin_elastic_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'cauchy_strain' : [['state', {'strain' : (-1,-1)}]]}

    function = staticmethod(terms.dw_lin_elastic)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        n_el, n_qp = self.data_shape[:2]

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat = mats[0]
            mat = nm.tile(mat, (n_el, n_qp, 1, 1))
            return (ts.dt, aux, mat, vg), shape, mode

        else:
            cache = self.get_cache( 'cauchy_strain', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    strain = cache( 'strain', self.get_current_group(), ii,
                                    state = state, get_vector = self.get_vector )
                    mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                    yield ii, (ts.dt, strain, mat, vg)
            return iter_kernel, shape, mode

class LinearElasticETHTerm(VectorVector, Term):
    r"""
    :Description:
    This term has the same definition as dw_lin_elastic_th, but assumes an
    exponential approximation of the convolution kernel resulting in much
    higher efficiency. Can use derivatives.

    :Definition:
    .. math::
        \int_{\Omega} \left [\int_0^t
        \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{u}(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v})

    :Arguments:
    material_0 : :math:`\Hcal_{ijkl}(0)`,
    material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
    """
    name = 'dw_lin_elastic_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'cauchy_strain' : [['state']],
                  'exp_history' : [['material_0', 'material_1', 'state']]}

    function = staticmethod(terms.dw_lin_elastic)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat0, mat1, virtual, state = self.get_args(**kwargs)
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if diff_var is None:
            cache = self.get_cache('cauchy_strain', 0)
            strain = cache('strain', self.get_current_group(), 0,
                           state=state, get_vector=self.get_vector)

            cache = self.get_cache('exp_history', 0)
            increment = cache('increment', self.get_current_group(), 0,
                              decay=mat1, values=strain)
            history = cache('history', self.get_current_group(), 0)

            fargs = (ts.dt, history + increment, mat0, vg)
            if ts.step == 0: # Just init the history in step 0.
                raise StopIteration

        else:
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            fargs = (ts.dt, aux, mat0, vg)

##        self.check_mat_shape( mat )

        return fargs, shape, mode

class LinearPrestressTerm(VectorVector, Term):
    r"""
    :Description:
    Linear prestress term, with the prestress :math:`\sigma_{ij}` given in
    the usual vector form exploiting symmetry: in 3D it has 6 components
    with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has
    3 components with the indices ordered as :math:`[11, 22, 12]`. Can be
    evaluated.

    :Definition:
    .. math::
        \int_{\Omega} \sigma_{ij} e_{ij}(\ul{v})
    """
    name = 'dw_lin_prestress'
    arg_types = (('material', 'virtual'),
                 ('material', 'parameter'))
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter')])
    modes = ('weak', 'eval')

    def check_mat_shape(self, mat):
        dim = self.data_shape[2]
        sym = (dim + 1) * dim / 2
        assert_(mat.shape == (self.data_shape[0], self.data_shape[1], sym, 1))

    def get_fargs_weak(self, diff_var=None, chunk_size=None, **kwargs):
        mat, virtual = self.get_args(**kwargs)
        ap, vg = virtual.get_approximation(self.get_current_group(), 'Volume')

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        if diff_var is not None:
            raise StopIteration

        self.check_mat_shape(mat)

        return (mat, vg), shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par = self.get_args(**kwargs)
        ap, vg = par.get_approximation(self.get_current_group(), 'Volume')

        self.set_data_shape(ap)

        self.check_mat_shape(mat)

        cache = self.get_cache('cauchy_strain', 0)
        strain = cache('strain', self.get_current_group(), 0,
                       state=par, get_vector=self.get_vector)

        return (strain, mat, vg), (chunk_size, 1, 1, 1), 0

    def d_lin_prestress(self, out, strain, mat, vg, chunk):
        aux = (mat[chunk] * strain[chunk]).sum(axis=2)
        aux.shape = aux.shape + (1,)

        status = vg.integrate_chunk(out, aux, chunk)
        return status

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_lin_prestress
            use_method_with_name(self, self.get_fargs_weak, 'get_fargs')

        else:
            self.function = self.d_lin_prestress
            use_method_with_name(self, self.get_fargs_eval, 'get_fargs')
            self.use_caches = {'cauchy_strain' : [['parameter']]}

class CauchyStrainTerm( Term ):
    r"""
    :Description:
    Cauchy strain tensor averaged in elements.
    
    :Definition:
    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \ull{e}(\ul{w}) /
        \int_{T_K} 1
    """
    name = 'de_cauchy_strain'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]

    function = staticmethod(terms.de_cauchy_strain)

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
    r"""
    :Description:
    Cauchy stress tensor averaged in elements.

    :Definition:
    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} D_{ijkl} e_{kl}(\ul{w}) / \int_{T_K} 1
    """
    name = 'de_cauchy_stress'
    arg_types = ('material', 'parameter')
    geometry = [(Volume, 'parameter')]
    use_caches = {'cauchy_strain' : [['parameter']]}

    function = staticmethod(terms.de_cauchy_stress)

    def build_c_fun_args( self, state, ap, vg, **kwargs ):
        mat, = self.get_args( ['material'], **kwargs )
        cache = self.get_cache( 'cauchy_strain', 0 )
        strain = cache( 'strain', self.get_current_group(), 0,
                        state = state, get_vector = self.get_vector )
        return strain, mat, vg

class CauchyStrainQTerm(Term):
    r"""
    :Description:
    Cauchy strain tensor in quadrature points, given in the usual vector form
    exploiting symmetry: in 3D it has 6 components with the indices ordered as
    :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has 3 components with the
    indices ordered as :math:`[11, 22, 12]`. The last three (non-diagonal)
    components are doubled so that it is energetically conjugate to the Cauchy
    stress tensor with the same storage.
    
    :Definition:
    .. math::
        \ull{e}(\ul{w})|_{qp}
    """
    name = 'dq_cauchy_strain'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'cauchy_strain' : [['parameter']]}

    def get_strain(self, diff_var=None, chunk_size=None, **kwargs):
        if diff_var is not None:
            raise StopIteration

        par, = self.get_args(['parameter'], **kwargs)

        cache = self.get_cache('cauchy_strain', 0)
        strain = cache('strain', self.get_current_group(), 0,
                       state=par, get_vector=self.get_vector)

        shape = (chunk_size,) + strain.shape[1:]

        return strain, shape

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        strain, shape = self.get_strain(diff_var, chunk_size, **kwargs)
        
        for out, chunk in self.char_fun(chunk_size, shape):
            yield strain[chunk], chunk, 0

class CauchyStressQTerm(CauchyStrainQTerm):
    r"""
    :Description:
    Cauchy stress tensor in quadrature points, given in the usual vector form
    exploiting symmetry: in 3D it has 6 components with the indices ordered as
    :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has 3 components with the
    indices ordered as :math:`[11, 22, 12]`.
    
    :Definition:
    .. math::
        D_{ijkl} e_{kl}(\ul{w})|_{qp}
    """
    name = 'dq_cauchy_stress'
    arg_types = ('material', 'parameter')
    geometry = [(Volume, 'parameter')]
    use_caches = {'cauchy_strain' : [['parameter']]}

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        strain, shape = self.get_strain(diff_var, chunk_size, **kwargs)

        mat, = self.get_args(['material'], **kwargs)

        for out, chunk in self.char_fun(chunk_size, shape):
            mc = mat[chunk]
            sc = strain[chunk]

            mc.shape = (mc.shape[0] * mc.shape[1],) + mc.shape[2:]
            sc.shape = (sc.shape[0] * sc.shape[1],) + sc.shape[2:]

            stress = nm.sum(mc * sc, axis=1)
            stress.shape = strain[chunk].shape

            yield stress, chunk, 0
