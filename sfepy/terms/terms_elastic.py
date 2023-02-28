import numpy as nm

from sfepy.linalg import dot_sequences
from sfepy.homogenization.utils import iter_sym
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_th import THTerm, ETHTerm

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
    General linear elasticity term, with :math:`D_{ijkl}` given in
    the usual matrix form exploiting symmetry: in 3D it is :math:`6\times6`
    with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it is
    :math:`3\times3` with the indices ordered as :math:`[11, 22, 12]`. Can be
    evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})

    :Arguments 1:
        - material : :math:`D_{ijkl}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`

    :Arguments 2:
        - material    : :math:`D_{ijkl}`
        - parameter_1 : :math:`\ul{w}`
        - parameter_2 : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'S, S', 'virtual' : ('D', 'state'),
                  'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'}
    modes = ('weak', 'eval')
##     symbolic = {'expression': expr,
##                 'map' : {'u' : 'state', 'D_sym' : 'material'}}

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

class LinearElasticIsotropicTerm(LinearElasticTerm):
    r"""
    Isotropic linear elasticity term.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})\\ \mbox{ with } \\
        D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
        \lambda \ \delta_{ij} \delta_{kl}

    :Arguments:
        - material_1: :math:`\lambda`
        - material_2: :math:`\mu`
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_iso'
    arg_types = (('material_1', 'material_2', 'virtual', 'state'),
                 ('material_1', 'material_2', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D',
                  'parameter_1' : 'D', 'parameter_2' : 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']

    def get_fargs(self, lam, mu, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.mechanics.matcoefs import stiffness_from_lame

        mat = stiffness_from_lame(self.region.dim, lam, mu)[:, :, 0, 0, :, :]
        return LinearElasticTerm.get_fargs(self, mat, virtual, state,
                                           mode=mode, term_mode=term_mode,
                                           diff_var=diff_var, **kwargs)

    def get_eval_shape(self, mat1, mat2, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        return LinearElasticTerm.get_eval_shape(self, None, None, state)

class SDLinearElasticTerm(Term):
    r"""
    Sensitivity analysis of the linear elastic term.

    :Definition:

    .. math::
        \int_{\Omega} \hat{D}_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})

    .. math::
        \hat{D}_{ijkl} = D_{ijkl}(\nabla \cdot \ul{\Vcal})
        - D_{ijkq}{\partial \Vcal_l \over \partial x_q}
        - D_{iqkl}{\partial \Vcal_j \over \partial x_q}

    :Arguments:
        - material    : :math:`D_{ijkl}`
        - parameter_w : :math:`\ul{w}`
        - parameter_u : :math:`\ul{u}`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_lin_elastic'
    arg_types = ('material', 'parameter_w', 'parameter_u',
                 'parameter_mv')
    arg_shapes = {'material' : 'S, S',
                  'parameter_w' : 'D', 'parameter_u' : 'D',
                  'parameter_mv' : 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']
    function = terms.d_sd_lin_elastic

    def get_fargs(self, mat, par_w, par_u, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        grad_w = self.get(par_w, 'grad')
        grad_u = self.get(par_u, 'grad')
        grad_mv = self.get(par_mv, 'grad')

        return 1.0, grad_w, grad_u, grad_mv, mat, vg

    def get_eval_shape(self, mat, par_w, par_u, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)

        return (n_el, 1, 1, 1), par_u.dtype

class LinearElasticTHTerm(THTerm):
    r"""
    Fading memory linear elastic (viscous) term. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} \left [\int_0^t
        \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{u}(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v})

    :Arguments:
        - ts       : :class:`TimeStepper` instance
        - material : :math:`\Hcal_{ijkl}(\tau)`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    arg_shapes = {'material' : '.: N, S, S',
                  'virtual' : ('D', 'state'), 'state' : 'D'}

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
    This term has the same definition as dw_lin_elastic_th, but assumes an
    exponential approximation of the convolution kernel resulting in much
    higher efficiency. Can use derivatives.

    :Definition:

    .. math::
        \int_{\Omega} \left [\int_0^t
        \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{u}(\tau)) \difd{\tau}
        \right]\,e_{ij}(\ul{v})

    :Arguments:
        - ts         : :class:`TimeStepper` instance
        - material_0 : :math:`\Hcal_{ijkl}(0)`
        - material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'virtual', 'state')
    arg_shapes = {'material_0' : 'S, S', 'material_1' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}

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

class LinearPrestressTerm(Term):
    r"""
    Linear prestress term, with the prestress :math:`\sigma_{ij}` given either
    in the usual vector form exploiting symmetry: in 3D it has 6 components
    with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it has
    3 components with the indices ordered as :math:`[11, 22, 12]`, or in the
    matrix (possibly non-symmetric) form. Can be evaluated.

    :Definition:

    .. math::
        \int_{\Omega} \sigma_{ij} e_{ij}(\ul{v})

    :Arguments 1:
        - material : :math:`\sigma_{ij}`
        - virtual  : :math:`\ul{v}`

    :Arguments 2:
        - material : :math:`\sigma_{ij}`
        - parameter : :math:`\ul{u}`
    """
    name = 'dw_lin_prestress'
    arg_types = (('material', 'virtual'),
                 ('material', 'parameter'))
    arg_shapes = [{'material' : 'S, 1', 'virtual' : ('D', None),
                  'parameter' : 'D'},
                  {'material' : 'D, D'}]
    modes = ('weak', 'eval')

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        sh = mat.shape
        is_nonsym = sh[2] == sh[3] == vg.dim and not(vg.dim == 1)

        if is_nonsym:
            mat = mat.reshape(sh[:2] + (vg.dim**2, 1))

        if mode == 'weak':
            return mat, vg

        else:
            if is_nonsym:
                strain = self.get(virtual, 'grad').transpose((0,1,3,2))
                nel, nqp, nr, nc = strain.shape
                strain = strain.reshape((nel, nqp, nr*nc, 1))
            else:
                strain = self.get(virtual, 'cauchy_strain')

            fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)
            return strain, mat, vg, fmode

    def get_eval_shape(self, mat, virtual,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, 1, 1), virtual.dtype

    def d_lin_prestress(self, out, strain, mat, vg, fmode):
        aux = dot_sequences(mat, strain, mode='ATB')
        if fmode == 2:
            out[:] = aux
            status = 0

        else:
            status = vg.integrate(out, aux, fmode)

        return status

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_lin_prestress

        else:
            self.function = self.d_lin_prestress

class LinearStrainFiberTerm(Term):
    r"""
    Linear (pre)strain fiber term with the unit direction vector
    :math:`\ul{d}`.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl} e_{ij}(\ul{v}) \left(d_k d_l\right)

    :Arguments:
        - material_1 : :math:`D_{ijkl}`
        - material_2 : :math:`\ul{d}`
        - virtual  : :math:`\ul{v}`
    """
    name = 'dw_lin_strain_fib'
    arg_types = ('material_1', 'material_2', 'virtual')
    arg_shapes = {'material_1' : 'S, S', 'material_2' : 'D, 1',
                  'virtual' : ('D', None)}

    function = staticmethod(terms.dw_lin_strain_fib)

    def get_fargs(self, mat1, mat2, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        omega = nm.empty(mat1.shape[:3] + (1,), dtype=nm.float64)
        for ii, (ir, ic) in enumerate(iter_sym(mat2.shape[2])):
            omega[..., ii, 0] = mat2[..., ir, 0] * mat2[..., ic, 0]

        return mat1, omega, vg

class CauchyStrainTerm(Term):
    r"""
    Evaluate Cauchy strain tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22,
    12]`. The last three (non-diagonal) components are doubled so that it is
    energetically conjugate to the Cauchy stress tensor with the same storage.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\cal{D}} \ull{e}(\ul{w})

    :Arguments:
        - parameter : :math:`\ul{w}`
    """
    name = 'ev_cauchy_strain'
    arg_types = ('parameter',)
    arg_shapes = {'parameter' : 'D'}
    integration = ('cell', 'facet_extra')

    @staticmethod
    def function(out, strain, vg, fmode):
        if fmode == 2:
            out[:] = strain
            status = 0

        else:
            status = terms.de_cauchy_strain(out, strain, vg.cmap, fmode)

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

        return (n_el, n_qp, dim * (dim + 1) // 2, 1), parameter.dtype


class CauchyStressTerm(Term):
    r"""
    Evaluate Cauchy stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\cal{D}} D_{ijkl} e_{kl}(\ul{w})

    :Arguments:
        - material  : :math:`D_{ijkl}`
        - parameter : :math:`\ul{w}`
    """
    name = 'ev_cauchy_stress'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : 'S, S', 'parameter' : 'D'}
    integration = ('cell', 'facet_extra')

    @staticmethod
    def function(out, coef, strain, mat, vg, fmode):
        if fmode == 2:
            out[:] = dot_sequences(mat, strain)
            status = 0

        else:
            status = terms.de_cauchy_stress(out, strain, mat, vg.cmap, fmode)

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

        return (n_el, n_qp, dim * (dim + 1) // 2, 1), parameter.dtype

class CauchyStressTHTerm(CauchyStressTerm, THTerm):
    r"""
    Evaluate fading memory Cauchy stress tensor.

    It is given in the usual vector form exploiting symmetry: in 3D it has 6
    components with the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in
    2D it has 3 components with the indices ordered as :math:`[11, 22, 12]`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\Omega} \int_0^t \Hcal_{ijkl}(t-\tau)\,e_{kl}(\ul{w}(\tau))
        \difd{\tau}

    :Arguments:
        - ts        : :class:`TimeStepper` instance
        - material  : :math:`\Hcal_{ijkl}(\tau)`
        - parameter : :math:`\ul{w}`
    """
    name = 'ev_cauchy_stress_th'
    arg_types = ('ts', 'material', 'parameter')
    arg_shapes = {'material' : '.: N, S, S', 'parameter' : 'D'}

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

    :Arguments:
        - ts         : :class:`TimeStepper` instance
        - material_0 : :math:`\Hcal_{ijkl}(0)`
        - material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
        - parameter  : :math:`\ul{w}`
    """
    name = 'ev_cauchy_stress_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'parameter')
    arg_shapes = {'material_0' : 'S, S', 'material_1' : '1, 1',
                  'parameter' : 'D'}

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

class NonsymElasticTerm(Term):
    r"""
    Elasticity term with non-symmetric gradient. The indices of matrix
    :math:`D_{ijkl}` are ordered as
    :math:`[11, 12, 13, 21, 22, 23, 31, 32, 33]` in 3D and as
    :math:`[11, 12, 21, 22]` in 2D.

    :Definition:

    .. math::
        \int_{\Omega} \ull{D} \nabla\ul{u} : \nabla\ul{v}

    :Arguments 1:
        - material : :math:`\ull{D}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`

    :Arguments 2:
        - material    : :math:`\ull{D}`
        - parameter_1 : :math:`\ul{w}`
        - parameter_2 : :math:`\ul{u}`
    """
    name = 'dw_nonsym_elastic'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : 'D2, D2', 'virtual' : ('D', 'state'),
                  'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'}
    modes = ('weak', 'eval')
    geometries = ['2_3', '2_4', '3_4', '3_8']

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mode == 'weak':
            if diff_var is None:
                grad = self.get(state, 'grad').transpose((0,1,3,2))
                nel, nqp, nr, nc = grad.shape
                grad = grad.reshape((nel,nqp,nr*nc,1))
                fmode = 0

            else:
                grad = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return grad, mat, vg, fmode

        elif mode == 'eval':
            grad1 = self.get(virtual, 'grad').transpose((0,1,3,2))
            grad2 = self.get(state, 'grad').transpose((0,1,3,2))
            nel, nqp, nr, nc = grad1.shape

            return 1.0,\
                   grad1.reshape((nel,nqp,nr*nc,1)),\
                   grad2.reshape((nel,nqp,nr*nc,1)),\
                   mat, vg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_nonsym_elastic

        else:
            self.function = terms.d_lin_elastic

def _build_wave_strain_op(vec, bf):
    dim = len(vec)

    if dim == 2:
        n0, n1 = vec
        nmat = nm.array([[n0, 0],
                         [0, n1],
                         [n1, n0]], dtype=nm.float64)

    else:
        n0, n1, n2 = vec
        nmat = nm.array([[n0, 0, 0],
                         [0, n1, 0],
                         [0, 0, n2],
                         [n1, n0, 0],
                         [n2, 0, n0],
                         [0, n2, n1]], dtype=nm.float64)

    out = nm.einsum('ik,cqkj->cqij', nmat, bf)
    return out

from sfepy.base.compat import block

def _build_cauchy_strain_op(bfg):
    dim = bfg.shape[2]
    if dim == 2:
        g1, g2 = bfg[..., 0:1, :], bfg[..., 1:2, :]
        zz = nm.zeros_like(g1)
        out = block([[g1, zz],
                     [zz, g2],
                     [g2, g1]])

    else:
        g1, g2, g3 = bfg[..., 0:1, :], bfg[..., 1:2, :], bfg[..., 2:3, :]
        zz = nm.zeros_like(g1)
        out = block([[g1, zz, zz],
                     [zz, g2, zz],
                     [zz, zz, g3],
                     [g2, g1, zz],
                     [g3, zz, g1],
                     [zz, g3, g2]])

    return out

class ElasticWaveTerm(Term):
    r"""
    Elastic dispersion term involving the wave strain :math:`g_{ij}`,
    :math:`g_{ij}(\ul{u}) = \frac{1}{2}(u_i \kappa_j + \kappa_i u_j)`, with the
    wave vector :math:`\ul{\kappa}`. :math:`D_{ijkl}` is given in the usual
    matrix form exploiting symmetry: in 3D it is :math:`6\times6` with the
    indices ordered as :math:`[11, 22, 33, 12, 13, 23]`, in 2D it is
    :math:`3\times3` with the indices ordered as :math:`[11, 22, 12]`.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ g_{ij}(\ul{v}) g_{kl}(\ul{u})

    :Arguments:
        - material_1 : :math:`D_{ijkl}`
        - material_2 : :math:`\ul{\kappa}`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_elastic_wave'
    arg_types = ('material_1', 'material_2', 'virtual', 'state')
    arg_shapes = {'material_1' : 'S, S', 'material_2' : '.: D',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']

    @staticmethod
    def function(out, out_qp, geo, fmode):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, mat, kappa, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import create_adof_conn, expand_basis

        geo, _ = self.get_mapping(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        ebf = expand_basis(geo.bf, dim)

        mat = Term.tile_mat(mat, n_el)
        gmat = _build_wave_strain_op(kappa, ebf)

        if diff_var is None:
            econn = state.field.get_econn('cell', self.region)
            adc = create_adof_conn(nm.arange(state.n_dof, dtype=nm.int32),
                                   econn, n_c, 0)
            vals = state()[adc]
            # Same as nm.einsum('qij,cj->cqi', gmat[0], vals)[..., None]
            aux = dot_sequences(gmat, vals[:, None, :, None])
            out_qp = dot_sequences(gmat, dot_sequences(mat, aux), 'ATB')
            fmode = 0

        else:
            out_qp = dot_sequences(gmat, dot_sequences(mat, gmat), 'ATB')
            fmode = 1

        return out_qp, geo, fmode

class ElasticWaveCauchyTerm(Term):
    r"""
    Elastic dispersion term involving the wave strain :math:`g_{ij}`,
    :math:`g_{ij}(\ul{u}) = \frac{1}{2}(u_i \kappa_j + \kappa_i u_j)`, with the
    wave vector :math:`\ul{\kappa}` and the elastic strain :math:`e_{ij}`.
    :math:`D_{ijkl}` is given in the usual matrix form exploiting symmetry: in
    3D it is :math:`6\times6` with the indices ordered as :math:`[11, 22, 33,
    12, 13, 23]`, in 2D it is :math:`3\times3` with the indices ordered as
    :math:`[11, 22, 12]`.

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ g_{ij}(\ul{v}) e_{kl}(\ul{u})\\
        \int_{\Omega} D_{ijkl}\ g_{ij}(\ul{u}) e_{kl}(\ul{v})

    :Arguments 1:
        - material_1 : :math:`D_{ijkl}`
        - material_2 : :math:`\ul{\kappa}`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`

    :Arguments 2:
        - material_1 : :math:`D_{ijkl}`
        - material_2 : :math:`\ul{\kappa}`
        - state      : :math:`\ul{u}`
        - virtual    : :math:`\ul{v}`
    """
    name = 'dw_elastic_wave_cauchy'
    arg_types = (('material_1', 'material_2', 'virtual', 'state'),
                 ('material_1', 'material_2', 'state', 'virtual'))
    arg_shapes = {'material_1' : 'S, S', 'material_2' : '.: D',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']
    modes = ('ge', 'eg')

    @staticmethod
    def function(out, out_qp, geo, fmode):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, mat, kappa, gvar, evar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import create_adof_conn, expand_basis

        geo, _ = self.get_mapping(evar)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(gvar)

        ebf = expand_basis(geo.bf, dim)

        mat = Term.tile_mat(mat, n_el)
        gmat = _build_wave_strain_op(kappa, ebf)
        emat = _build_cauchy_strain_op(geo.bfg)

        if diff_var is None:
            avar = evar if self.mode == 'ge' else gvar
            econn = avar.field.get_econn('cell', self.region)
            adc = create_adof_conn(nm.arange(avar.n_dof, dtype=nm.int32),
                                   econn, n_c, 0)
            vals = avar()[adc]

            if self.mode == 'ge':
                # Same as aux = self.get(avar, 'cauchy_strain'),
                aux = dot_sequences(emat, vals[:, None, :, None])
                out_qp = dot_sequences(gmat, dot_sequences(mat, aux), 'ATB')

            else:
                aux = dot_sequences(gmat, vals[:, None, :, None])
                out_qp = dot_sequences(emat, dot_sequences(mat, aux), 'ATB')

            fmode = 0

        else:
            if self.mode == 'ge':
                out_qp = dot_sequences(gmat, dot_sequences(mat, emat), 'ATB')

            else:
                out_qp = dot_sequences(emat, dot_sequences(mat, gmat), 'ATB')

            fmode = 1

        return out_qp, geo, fmode
