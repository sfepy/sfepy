from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms_hyperelastic_tl import HyperElasticTLBase
from sfepy.mechanics.tensors import dim2sym
from sfepy.homogenization.utils import iter_sym
from six.moves import range

def create_omega(fdir):
    r"""
    Create the fibre direction tensor :math:`\omega_{ij} = d_i d_j`.
    """
    n_el, n_qp, dim, _ = fdir.shape
    sym = dim2sym(dim)
    omega = nm.empty((n_el, n_qp, sym, 1), dtype=nm.float64)
    for ii, (ir, ic) in enumerate(iter_sym(dim)):
        omega[..., ii, 0] = fdir[..., ir, 0] * fdir[..., ic, 0]

    return omega

def compute_fibre_strain(green_strain, omega):
    """
    Compute the Green strain projected to the fibre direction.
    """
    eps = nm.zeros_like(omega[..., :1, :])
    for ii in range(omega.shape[2]):
        eps[..., 0, 0] += omega[..., ii, 0] * green_strain[..., ii, 0]

    return eps

def _setdefault_fibre_data(self, state):
    """
    Returns `fibre_data` :class:``Struct`` for storing/caching fibre-related
    data.
    """
    cache = self.set_default('fibre_cache', {})

    _, _, key = self.get_mapping(state, return_key=True)

    data_key = key + ('fibre_data',)
    if data_key in cache:
        fibre_data = cache[data_key]

    else:
        fibre_data = Struct()
        cache[data_key] = fibre_data

    return fibre_data

class FibresActiveTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic active fibres term. Effective stress
    :math:`S_{ij} = A f_{\rm max} \exp{\left\{-(\frac{\epsilon -
    \varepsilon_{\rm opt}}{s})^2\right\}} d_i d_j`,
    where :math:`\epsilon = E_{ij} d_i d_j` is the Green strain
    :math:`\ull{E}` projected to the fibre direction :math:`\ul{d}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material_1 : :math:`f_{\rm max}`
        - material_2 : :math:`\varepsilon_{\rm opt}`
        - material_3 : :math:`s`
        - material_4 : :math:`\ul{d}`
        - material_5 : :math:`A`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_tl_fib_a'
    arg_types = ('material_1', 'material_2', 'material_3',
                 'material_4', 'material_5', 'virtual', 'state')
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'material_3' : '1, 1', 'material_4' : 'D, 1',
                  'material_5' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    family_data_names = ['green_strain']

    def get_fargs(self, mat1, mat2, mat3, mat4, mat5, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        fibre_data = _setdefault_fibre_data(self, state)

        n_el, _, _, _, _ = self.get_data_shape(state)
        mat1 = HyperElasticTLBase.tile_mat(mat1, n_el)
        mat2 = HyperElasticTLBase.tile_mat(mat2, n_el)
        mat3 = HyperElasticTLBase.tile_mat(mat3, n_el)
        mat4 = HyperElasticTLBase.tile_mat(mat4, n_el)
        mat5 = HyperElasticTLBase.tile_mat(mat5, n_el)

        fargs = HyperElasticTLBase.get_fargs(self,
                                             (mat1, mat2, mat3, mat4, mat5),
                                             virtual, state,
                                             mode, term_mode, diff_var,
                                             fibre_data=fibre_data,
                                             **kwargs)
        return fargs

    @staticmethod
    def stress_function(out, pars, green_strain,
                        fibre_data=None):
        fmax, eps_opt, s, fdir, act = pars

        omega = fibre_data.get('omega', None)
        if omega is None:
            omega = fibre_data.omega = create_omega(fdir)

        eps = fibre_data.eps = compute_fibre_strain(green_strain, omega)

        tau = fibre_data.tau = act * fmax * nm.exp(-((eps - eps_opt) / s)**2.0)

        out[:] = omega * tau

    @staticmethod
    def tan_mod_function(out, pars, green_strain,
                         fibre_data=None):
        fmax, eps_opt, s, fdir, act = pars

        omega, eps, tau = fibre_data.omega, fibre_data.eps, fibre_data.tau

        for ir in range(omega.shape[2]):
            for ic in range(omega.shape[2]):
                out[..., ir, ic] = omega[..., ir, 0] * omega[..., ic, 0]

        out[:] *= -2.0 * ((eps - eps_opt) / (s**2.0)) * tau

    def get_eval_shape(self, mat1, mat2, mat3, mat4, mat5, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim2sym(dim)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, sym, 1), state.dtype

class FibresExponentialTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic fibres term with an exponential response. Effective stress
    :math:`S_{ij} = \max\left(0, \sigma \left[ \exp{\left\{k (\epsilon -
    \epsilon_0)\right\}} - 1 \right]\right) d_i d_j`, where :math:`\epsilon =
    E_{ij} d_i d_j` is the Green strain :math:`\ull{E}` projected to the fibre
    direction :math:`\ul{d}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material_1 : :math:`\sigma`
        - material_3 : :math:`k`
        - material_3 : :math:`\epsilon_{0}`
        - material_4 : :math:`\ul{d}`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_tl_fib_e'
    arg_types = ('material_1', 'material_2', 'material_3',
                 'material_4', 'virtual', 'state')
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'material_3' : '1, 1', 'material_4' : 'D, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    family_data_names = ['green_strain']

    def get_fargs(self, mat1, mat2, mat3, mat4, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        fibre_data = _setdefault_fibre_data(self, state)

        n_el, _, _, _, _ = self.get_data_shape(state)
        mat1 = HyperElasticTLBase.tile_mat(mat1, n_el)
        mat2 = HyperElasticTLBase.tile_mat(mat2, n_el)
        mat3 = HyperElasticTLBase.tile_mat(mat3, n_el)
        mat4 = HyperElasticTLBase.tile_mat(mat4, n_el)

        fargs = HyperElasticTLBase.get_fargs(self,
                                             (mat1, mat2, mat3, mat4),
                                             virtual, state,
                                             mode, term_mode, diff_var,
                                             fibre_data=fibre_data,
                                             **kwargs)
        return fargs

    @staticmethod
    def stress_function(out, pars, green_strain,
                        fibre_data=None):
        sigma, k, eps0, fdir = pars

        omega = fibre_data.get('omega', None)
        if omega is None:
            omega = fibre_data.omega = create_omega(fdir)

        eps = fibre_data.eps = compute_fibre_strain(green_strain, omega)

        tau = fibre_data.tau = nm.maximum(
            0,
            sigma * (nm.exp(k * (eps - eps0)) - 1.0),
        )

        out[:] = omega * tau

    @staticmethod
    def tan_mod_function(out, pars, green_strain,
                         fibre_data=None):
        sigma, k, eps0, fdir = pars

        omega, eps, tau = fibre_data.omega, fibre_data.eps, fibre_data.tau

        for ir in range(omega.shape[2]):
            for ic in range(omega.shape[2]):
                out[..., ir, ic] = omega[..., ir, 0] * omega[..., ic, 0]

        tan_mod = nm.where(
            tau > 0.0,
            sigma * k * nm.exp(k * (eps - eps0)),
            0.0,
        )
        out[:] *= tan_mod

    def get_eval_shape(self, mat1, mat2, mat3, mat4, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim2sym(dim)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, sym, 1), state.dtype
