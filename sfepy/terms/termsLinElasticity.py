import numpy as nm

from sfepy.base.base import use_method_with_name, assert_
from sfepy.terms.terms import Term, terms
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

class LinearElasticTerm(Term):
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

    :Arguments 1:
        material : :math:`D_{ijkl}`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`

    :Arguments 2:
        material    : :math:`D_{ijkl}`,
        parameter_1 : :math:`\ul{w}`,
        parameter_2 : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')
##     symbolic = {'expression': expr,
##                 'map' : {'u' : 'state', 'D_sym' : 'material'}}

    def check_shapes(self, mat, virtual, state):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = (dim + 1) * dim / 2
        assert_(mat.shape == (n_el, n_qp, sym, sym))

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mode == 'weak':
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            if diff_var is None:
                strain = self.get(state, 'cauchy_strain')
                fmode = 0

            else:
                strain = aux
                fmode = 1

            return 1.0, strain, mat, vg, fmode

        elif mode == 'eval':
            strain1 = self.get(virtual, 'cauchy_strain')
            strain2 = self.get(state, 'cauchy_strain')

            return 1.0, strain1, strain2, mat, vg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_lin_elastic

        else:
            self.function = terms.d_lin_elastic

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
        virtual    : :math:`\ul{v}`,
        state      : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_iso'
    arg_types = ('material_1', 'material_2', 'virtual', 'state')

    function = staticmethod(terms.dw_lin_elastic_iso)

    def check_shapes(self, lam, mu, virtual, state):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        assert_(lam.shape == (n_el, n_qp, 1, 1))
        assert_(mu.shape == (n_el, n_qp, 1, 1))

    def get_fargs(self, lam, mu, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mode == 'weak':
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            if diff_var is None:
                strain = self.get(state, 'cauchy_strain')
                fmode = 0

            else:
                strain = aux
                fmode = 1

            return strain, lam, mu, vg, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class LinearElasticTHTerm( VectorVectorTH, Term ):
    r"""
    :Description:
    Fading memory linear elastic (viscous) term. Can use derivatives.

    :Definition:
    .. math::
        \int_{\Omega} \left [\int_0^t
        \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{u}(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v})

    :Arguments:
        ts       : :class:`TimeStepper` instance,
        material : :math:`\Hcal_{ijkl}(\tau)`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    use_caches = {'cauchy_strain' : [['state', {'strain' : (-1,-1)}]]}

    function = staticmethod(terms.dw_lin_elastic)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)

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
                    strain = cache('strain', self, ii,
                                   state=state, get_vector=self.get_vector)
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
        ts         : :class:`TimeStepper` instance,
        material_0 : :math:`\Hcal_{ijkl}(0)`,
        material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`),
        virtual    : :math:`\ul{v}`,
        state      : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'virtual', 'state')
    use_caches = {'cauchy_strain' : [['state']],
                  'exp_history' : [['material_0', 'material_1', 'state']]}

    function = staticmethod(terms.dw_lin_elastic)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat0, mat1, virtual, state = self.get_args(**kwargs)
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if diff_var is None:
            cache = self.get_cache('cauchy_strain', 0)
            strain = cache('strain', self, 0,
                           state=state, get_vector=self.get_vector)

            cache = self.get_cache('exp_history', 0)
            increment = cache('increment', self, 0,
                              decay=mat1, values=strain)
            history = cache('history', self, 0)

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

    :Arguments 1:
        material : :math:`\sigma_{ij}`,
        virtual  : :math:`\ul{v}`

    :Arguments 2:
        material : :math:`\sigma_{ij}`,
        parameter : :math:`\ul{u}`
    """
    name = 'dw_lin_prestress'
    arg_types = (('material', 'virtual'),
                 ('material', 'parameter'))
    modes = ('weak', 'eval')

    def check_mat_shape(self, mat):
        dim = self.data_shape[2]
        sym = (dim + 1) * dim / 2
        assert_(mat.shape == (self.data_shape[0], self.data_shape[1], sym, 1))

    def get_fargs_weak(self, diff_var=None, chunk_size=None, **kwargs):
        mat, virtual = self.get_args(**kwargs)
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        if diff_var is not None:
            raise StopIteration

        self.check_mat_shape(mat)

        return (mat, vg), shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par = self.get_args(**kwargs)
        ap, vg = self.get_approximation(par)

        self.set_data_shape(ap)

        self.check_mat_shape(mat)

        cache = self.get_cache('cauchy_strain', 0)
        strain = cache('strain', self, 0,
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

class LinearStrainFiberTerm(VectorVector, Term):
    r"""
    :Description:
    Linear (pre)strain fiber term. Given fiber orientation :math:`\nu_{i}`.

    :Definition:
    .. math::
        \int_{\Omega} D_{ijkl} e_{ij}(\ul{v}) \left(\nu_i \nu_j\right)

    :Arguments:
        material_1 : :math:`D_{ijkl}`,
        material_2 : :math:`\nu_i`,
        virtual  : :math:`\ul{v}`

    """
    name = 'dw_lin_strain_fib'
    arg_types = ('material_1', 'material_2', 'virtual')

    function = staticmethod(terms.dw_lin_strain_fib)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, virtual = self.get_args(**kwargs)
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1)
        else:
            raise StopIteration

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, mat1, mat2, vg, chunk )
            yield out, chunk, status

class CauchyStrainTerm(Term):
    r"""
    :Description:
    Cauchy strain tensor averaged in elements.

    :Definition:
    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \ull{e}(\ul{w}) /
        \int_{T_K} 1

    :Arguments:
        parameter : :math:`\ul{w}`
    """
    name = 'de_cauchy_strain'
    arg_types = ('parameter',)

    function = staticmethod(terms.de_cauchy_strain)

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        strain = self.get(parameter, 'cauchy_strain')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return strain, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, dim * (dim + 1) / 2, 1), parameter.dtype

class CauchyStressTerm(Term):
    r"""
    :Description:
    Cauchy stress tensor averaged in elements.

    :Definition:
    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} D_{ijkl} e_{kl}(\ul{w}) / \int_{T_K} 1

    :Arguments:
        material  : :math:`D_{ijkl}`,
        parameter : :math:`\ul{w}`
    """
    name = 'de_cauchy_stress'
    arg_types = ('material', 'parameter')

    function = staticmethod(terms.de_cauchy_stress)

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        strain = self.get(parameter, 'cauchy_strain')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return strain, mat, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, dim * (dim + 1) / 2, 1), parameter.dtype

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

    :Arguments:
        parameter : :math:`\ul{w}`
   """
    name = 'dq_cauchy_strain'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, strain):
        out[:] = strain

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return (self.get(parameter, 'cauchy_strain'),)

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, n_qp, dim * (dim + 1) / 2, 1), parameter.dtype

class CauchyStressQTerm(Term):
    r"""
    :Description:
    Cauchy stress tensor in quadrature points, given in the usual vector form
    exploiting symmetry: in 3D it has 6 components with the indices ordered as
    :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has 3 components with the
    indices ordered as :math:`[11, 22, 12]`.

    :Definition:
    .. math::
        D_{ijkl} e_{kl}(\ul{w})|_{qp}

    :Arguments:
        material  : :math:`D_{ijkl}`,
        parameter : :math:`\ul{w}`
    """
    name = 'dq_cauchy_stress'
    arg_types = ('material', 'parameter')

    @staticmethod
    def function(out, mat, strain):
        mc = mat.reshape((mat.shape[0] * mat.shape[1],) + mat.shape[2:])
        sc = strain.reshape((strain.shape[0] * strain.shape[1],)
                            + strain.shape[2:])

        stress = nm.sum(mc * sc, axis=1)
        out[:] = stress.reshape(strain.shape)

        return 0

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return mat, self.get(parameter, 'cauchy_strain')

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, n_qp, dim * (dim + 1) / 2, 1), parameter.dtype
