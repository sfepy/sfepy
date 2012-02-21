import numpy as nm

from sfepy.base.base import use_method_with_name, assert_
from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_th import THTerm, ETHTerm
from sfepy.terms.terms_base import VectorVector

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
            if diff_var is None:
                strain = self.get(state, 'cauchy_strain')
                fmode = 0

            else:
                strain = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return 1.0, strain, mat, vg, fmode

        elif mode == 'eval':
            strain1 = self.get(virtual, 'cauchy_strain')
            strain2 = self.get(state, 'cauchy_strain')

            return 1.0, strain1, strain2, mat, vg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_lin_elastic

        else:
            self.function = terms.d_lin_elastic

class LinearElasticIsotropicTerm(Term):
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
            if diff_var is None:
                strain = self.get(state, 'cauchy_strain')
                fmode = 0

            else:
                strain = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return strain, lam, mu, vg, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class LinearElasticTHTerm(THTerm):
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

    function = staticmethod(terms.dw_lin_elastic)

    def get_fargs(self, ts, mats, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        if mode == 'weak':
            if diff_var is None:
                def iter_kernel():
                    for ii, mat in enumerate(mats):
                        strain = self.get(state, 'cauchy_strain',
                                          step=-ii)
                        mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                        yield ii, (ts.dt, strain, mat, vg, 0)
                fargs = iter_kernel

            else:
                strain = nm.array([0], ndmin=4, dtype=nm.float64)
                mat = nm.tile(mats[0], (n_el, n_qp, 1, 1))
                fargs = ts.dt, strain, mat, vg, 1

            return fargs

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class LinearElasticETHTerm(ETHTerm):
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

    function = staticmethod(terms.dw_lin_elastic)

    def get_fargs(self, ts, mat0, mat1, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _, key = self.get_mapping(state, return_key=True)

        if diff_var is None:
            strain = self.get(state, 'cauchy_strain')

            key += tuple(self.arg_names[ii] for ii in [1, 2, 4])
            data = self.get_eth_data(key, state, mat1, strain)

            fargs = (ts.dt, data.history + data.values, mat0, vg, 0)

        else:
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            fargs = (ts.dt, aux, mat0, vg, 1)

        return fargs

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
    Evaluate Cauchy strain tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22,
    12]`. The last three (non-diagonal) components are doubled so that it is
    energetically conjugate to the Cauchy stress tensor with the same storage.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:
    .. math::
        \int_{\Omega} \ull{e}(\ul{w})

    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \ull{e}(\ul{w}) /
        \int_{T_K} 1

    .. math::
        \ull{e}(\ul{w})|_{qp}

    :Arguments:
        parameter : :math:`\ul{w}`
    """
    name = 'ev_cauchy_strain'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, strain, vg, fmode):
        if fmode == 2:
            out[:] = strain
            status = 0

        else:
            status = terms.de_cauchy_strain(out, strain, vg, fmode)

        return status

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        strain = self.get(parameter, 'cauchy_strain')

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return strain, vg, fmode

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim * (dim + 1) / 2, 1), parameter.dtype

class CauchyStressTerm(Term):
    r"""
    :Description:
    Evaluate Cauchy stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:
    .. math::
        \int_{\Omega} D_{ijkl} e_{kl}(\ul{w})

    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} D_{ijkl} e_{kl}(\ul{w}) / \int_{T_K} 1

    .. math::
        D_{ijkl} e_{kl}(\ul{w})|_{qp}

    :Arguments:
        material  : :math:`D_{ijkl}`,
        parameter : :math:`\ul{w}`
    """
    name = 'ev_cauchy_stress'
    arg_types = ('material', 'parameter')

    @staticmethod
    def function(out, coef, strain, mat, vg, fmode):
        if fmode == 2:
            out[:] = dot_sequences(mat, strain)
            status = 0

        else:
            status = terms.de_cauchy_stress(out, strain, mat, vg, fmode)

        if coef is not None:
            out *= coef

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        strain = self.get(parameter, 'cauchy_strain')

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return None, strain, mat, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim * (dim + 1) / 2, 1), parameter.dtype

class CauchyStressTHTerm(CauchyStressTerm, THTerm):
    r"""
    :Description:
    Evaluate fading memory Cauchy stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:
    .. math::
        \int_{\Omega} \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau))
        \difd{\tau}

    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau))
        \difd{\tau} / \int_{T_K} 1

    .. math::
        \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau)) \difd{\tau}|_{qp}

    :Arguments:
        ts        : :class:`TimeStepper` instance,
        material  : :math:`\Hcal_{ijkl}(\tau)`,
        parameter : :math:`\ul{w}`
    """
    name = 'ev_cauchy_stress_th'
    arg_types = ('ts', 'material', 'parameter')

    def get_fargs(self, ts, mats, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)
        def iter_kernel():
            for ii, mat in enumerate(mats):
                strain = self.get(state, 'cauchy_strain',
                                  step=-ii)
                mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                yield ii, (ts.dt, strain, mat, vg, fmode)

        return iter_kernel

    def get_eval_shape(self, ts, mats, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        out = CauchyStressTerm.get_eval_shape(self, mats, parameter, mode,
                                              term_mode, diff_var, **kwargs)
        return out

class CauchyStressETHTerm(CauchyStressTerm, ETHTerm):
    r"""
    :Description:
    Evaluate fading memory Cauchy stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Assumes an exponential approximation of the convolution kernel resulting in
    much higher efficiency.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:
    .. math::
        \int_{\Omega} \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau))
        \difd{\tau}

    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau))
        \difd{\tau} / \int_{T_K} 1

    .. math::
        \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau)) \difd{\tau}|_{qp}

    :Arguments:
        ts         : :class:`TimeStepper` instance,
        material_0 : :math:`\Hcal_{ijkl}(0)`,
        material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`),
        parameter  : :math:`\ul{w}`
    """
    name = 'ev_cauchy_stress_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'parameter')

    def get_fargs(self, ts, mat0, mat1, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _, key = self.get_mapping(state, return_key=True)

        strain = self.get(state, 'cauchy_strain')

        key += tuple(self.arg_names[1:])
        data = self.get_eth_data(key, state, mat1, strain)

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return ts.dt, data.history + data.values, mat0, vg, fmode

    def get_eval_shape(self, ts, mat0, mat1, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        out = CauchyStressTerm.get_eval_shape(self, mat0, parameter, mode,
                                              term_mode, diff_var, **kwargs)
        return out
