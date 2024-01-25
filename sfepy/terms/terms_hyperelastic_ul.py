from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import terms
from sfepy.terms.terms_hyperelastic_base import\
    HyperElasticBase, HyperElasticFamilyData

_msg_missing_data = 'missing family data!'

class HyperElasticULFamilyData(HyperElasticFamilyData):
    """
    Family data for UL formulation.
    """
    family_function = staticmethod(terms.dq_finite_strain_ul)
    cache_name = 'ul_common'
    data_names = ('mtx_f', 'det_f', 'sym_b', 'tr_b', 'in2_b', 'green_strain')

class HyperElasticULBase(HyperElasticBase):
    """
    Base class for all hyperelastic terms in UL formulation family.

    The subclasses should have the following static method attributes:
    - `stress_function()` (the stress)
    - `tan_mod_function()` (the tangent modulus)
    """
    weak_function = staticmethod(terms.dw_he_rtm)
    hyperelastic_mode = 1
    get_family_data = HyperElasticULFamilyData()

class NeoHookeanULTerm(HyperElasticULBase):
    r"""
    Hyperelastic neo-Hookean term. Effective stress :math:`\tau_{ij} = \mu
    J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

    :Definition:

    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        - material : :math:`\mu`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_ul_he_neohook'
    family_data_names = ['det_f', 'tr_b', 'sym_b']

    stress_function = staticmethod(terms.dq_ul_he_stress_neohook)
    tan_mod_function = staticmethod(terms.dq_ul_he_tan_mod_neohook)

class MooneyRivlinULTerm(HyperElasticULBase):
    r"""
    Hyperelastic Mooney-Rivlin term.

    :Definition:

    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        - material : :math:`\kappa`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_ul_he_mooney_rivlin'
    family_data_names = ['det_f', 'tr_b', 'sym_b', 'in2_b']

    stress_function = staticmethod(terms.dq_ul_he_stress_mooney_rivlin)
    tan_mod_function = staticmethod(terms.dq_ul_he_tan_mod_mooney_rivlin)

class BulkPenaltyULTerm(HyperElasticULBase):
    r"""
    Hyperelastic bulk penalty term. Stress :math:`\tau_{ij} = K(J-1)\; J
    \delta_{ij}`.

    :Definition:

    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        - material : :math:`K`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_ul_bulk_penalty'
    family_data_names = ['det_f']

    stress_function = staticmethod(terms.dq_ul_he_stress_bulk)
    tan_mod_function = staticmethod(terms.dq_ul_he_tan_mod_bulk)

class BulkPressureULTerm(HyperElasticULBase):
    r"""
    Hyperelastic bulk pressure term. Stress :math:`S_{ij} = -p J \delta_{ij}`.

    :Definition:

    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        - virtual : :math:`\ul{v}`
        - state   : :math:`\ul{u}`
        - state_p : :math:`p`
    """

    name = 'dw_ul_bulk_pressure'
    arg_types = ('virtual', 'state', 'state_p')
    arg_geometry_types = {('state_p', None) : {'facet_extra' : 'facet'}}
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D', 'state_p' : 1}
    family_data_names = ['det_f', 'sym_b']

    family_function = staticmethod(terms.dq_finite_strain_ul)
    weak_function = staticmethod(terms.dw_he_rtm)
    weak_dp_function = staticmethod(terms.dw_ul_volume)

    stress_function = staticmethod(terms.dq_ul_stress_bulk_pressure)
    tan_mod_u_function = staticmethod(terms.dq_ul_tan_mod_bulk_pressure_u)

    def compute_data(self, family_data, mode, **kwargs):
        det_f, sym_b = family_data.det_f, family_data.sym_b
        p_qp = family_data.p_qp

        if mode == 0:
            out = nm.empty_like(sym_b)
            fun = self.stress_function

        elif mode == 1:
            shape = list(sym_b.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)
            fun = self.tan_mod_u_function

        else:
            raise ValueError('bad mode! (%d)' % mode)

        fun(out, p_qp, det_f)

        return out

    def get_fargs(self, virtual, state, state_p,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgv, _ = self.get_mapping(state)

        name = state.name
        fd = self.get_family_data(state, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])
        fd.p_qp = self.get(state_p, 'val')

        if mode == 'weak':
            if diff_var != state_p.name:
                if diff_var is None:
                    stress = self.compute_data(fd, 0, **kwargs)
                    self.stress_cache = stress
                    tan_mod = nm.array([0], ndmin=4, dtype=nm.float64)

                    fmode = 0

                else:
                    stress = self.stress_cache
                    if stress is None:
                        stress = self.compute_data(fd, 0, **kwargs)

                    tan_mod = self.compute_data(fd, 1, **kwargs)
                    fmode = 1

                fargs = (self.weak_function,
                         stress, tan_mod, fd.mtx_f, fd.det_f, vgv, fmode, 1)

            else:
                vgs, _ = self.get_mapping(state_p)

                fargs =  (self.weak_dp_function, fd.det_f, vgs, vgv, 1, -1)

            return fargs

        elif mode == 'el_avg':
            if term_mode == 'strain':
                out_qp = fd.green_strain

            elif term_mode == 'stress':
                out_qp = self.compute_data(fd, 0, **kwargs)

            else:
                raise ValueError('unsupported term mode in %s! (%s)'
                                 % (self.name, term_mode))

            return self.integrate, out_qp, vgv, 1

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, virtual, state, state_p,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = (dim + 1) * dim // 2

        return (n_el, 1, sym, 1), state.dtype

class VolumeULTerm(HyperElasticULBase):
    r"""
    Volume term (weak form) in the updated Lagrangian formulation.

    :Definition:

    .. math::
         \begin{array}{l}
         \int_{\Omega} q J(\ul{u}) \\
         \mbox{volume mode: vector for } K \from \Ical_h: \int_{T_K}
         J(\ul{u}) \\
         \mbox{rel\_volume mode: vector for } K \from \Ical_h:
         \int_{T_K} J(\ul{u}) / \int_{T_K} 1
         \end{array}

    :Arguments:
        - virtual : :math:`q`
        - state   : :math:`\ul{u}`
    """
    name = 'dw_ul_volume'
    arg_types = ('virtual', 'state')
    arg_geometry_types = {('virtual', None) : {'facet_extra' : 'facet'}}
    arg_shapes = {'virtual' : (1, None), 'state' : 'D'}
    family_data_names = ['mtx_f', 'det_f']

    function = staticmethod(terms.dw_ul_volume)
    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgs, _ = self.get_mapping(virtual)
        vgv, _ = self.get_mapping(state)

        name = state.name
        fd = self.get_family_data(state, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])

        if mode == 'weak':
            if diff_var is None:
                fmode = 0

            else:
                fmode = 1

        elif mode == 'eval':
            if term_mode == 'volume':
                fmode = 2

            elif term_mode == 'rel_volume':
                fmode = 3

            else:
                raise ValueError('unsupported term evaluation mode in %s! (%s)'
                                 % (self.name, term_mode))

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

        return fd.det_f, vgs, vgv, 0, fmode

    def get_eval_shape(self, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

class CompressibilityULTerm(HyperElasticULBase):
    r"""
    Compressibility term for the updated Lagrangian formulation

    :Definition:

    .. math::
        \int_{\Omega} 1\over \gamma p \, q

    :Arguments:
        - material : :math:`\gamma`
        - virtual  : :math:`q`
        - state    : :math:`p`
        - parameter_u  : :math:`\ul(u)`
    """
    name = 'dw_ul_compressible'
    arg_types = ('material', 'virtual', 'state', 'parameter_u')
    arg_geometry_types = {('virtual', None) : {'facet_extra' : 'facet'},
                          ('state', None) : {'facet_extra' : 'facet'}}
    arg_shapes = {'material' : '1, 1', 'virtual' : (1, 'state'), 'state' : 1,
                  'parameter_u' : 'D'}
    family_data_names = ['mtx_f', 'det_f']

    function = staticmethod(terms.dw_volume_dot_scalar)

    def get_fargs(self, bulk, virtual, state, parameter_u,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgp, _ = self.get_mapping(virtual)
        vgs, _ = self.get_mapping(state)
        vgu, _ = self.get_mapping(parameter_u)

        name = parameter_u.name
        fd = self.get_family_data(parameter_u, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])

        coef = nm.divide(bulk, fd.det_f)

        if mode == 'weak':
            if diff_var is None:
                val_qp = self.get(state, 'val')
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return coef, val_qp, vgp, vgs, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))
