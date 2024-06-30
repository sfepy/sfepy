from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_hyperelastic_base import\
    HyperElasticBase, HyperElasticFamilyData

class HyperElasticTLFamilyData(HyperElasticFamilyData):
    """
    Family data for TL formulation.
    """
    family_function = staticmethod(terms.dq_finite_strain_tl)
    cache_name = 'tl_common'
    data_names = ('mtx_f', 'det_f', 'sym_c', 'tr_c', 'in2_c', 'sym_inv_c',
                  'green_strain')

class HyperElasticTLBase(HyperElasticBase):
    """
    Base class for all hyperelastic terms in TL formulation family.

    The subclasses should have the following static method attributes:
    - `stress_function()` (the stress)
    - `tan_mod_function()` (the tangent modulus)

    The common (family) data are cached in the evaluate cache of state
    variable.
    """
    weak_function = staticmethod(terms.dw_he_rtm)
    hyperelastic_mode = 0
    get_family_data = HyperElasticTLFamilyData()

class NeoHookeanTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic neo-Hookean term. Effective stress
    :math:`S_{ij} = \mu J^{-\frac{2}{3}}(\delta_{ij} -
    \frac{1}{3}C_{kk}C_{ij}^{-1})`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`\mu`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_neohook'
    family_data_names = ['det_f', 'tr_c', 'sym_inv_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_neohook)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_neohook)

class GenYeohTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic generalized Yeoh term [1]. Effective stress
    :math:`S_{ij} = 2 p K (I_1 - 3)^{p-1} J^{-\frac{2}{3}}(\delta{ij} -
    \frac{1}{3}C_{kk}C{ij}^{-1})`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`p, K`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`

    [1] Travis W. Hohenberger, Richard J. Windslow, Nicola M. Pugno,
    James J. C. Busfield. Aconstitutive Model For Both Lowand High
    Strain Nonlinearities In Highly Filled Elastomers And
    Implementation With User-Defined Material Subroutines In
    Abaqus. Rubber Chemistry And Technology, Vol. 92, No. 4,
    Pp. 653-686 (2019)
    """
    name = 'dw_tl_he_genyeoh'
    family_data_names = ['det_f', 'tr_c', 'sym_inv_c']
    arg_shapes = {'material' : '1, 2',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['3_4', '3_8']

    @staticmethod
    def _get_single_stress(i_1, i_3, c_inv, coef, exp):
        _bracket = (i_3**(-1 / 3) * i_1 - 3.)
        _bracket_pow = _bracket**(exp - 1) \
            if (_bracket > 0. or exp >= 1.) else 1.
        out = 2 * exp * coef * _bracket_pow \
            * i_3**(-1 / 3) * (nm.eye(3) - i_1 * c_inv / 3)
        return out

    def stress_function(self, out, mat, *fargs, **kwargs):
        det_f, tr_c, inv_c = fargs
        mat = HyperElasticBase.tile_mat(mat, det_f.shape[0])
        coef, exp = mat[:, :, :, :1], mat[:, :, :, 1:]

        i_3 = det_f**2
        ident = nm.array([1., 1, 1, 0, 0, 0])

        n_cells, n_qps, _, _ = out.shape
        for cell in range(n_cells):
            for qp in range(n_qps):
                _inv_c = inv_c[cell, qp]
                _c_inv = nm.array([
                    [_inv_c[0], _inv_c[3], _inv_c[4]],
                    [_inv_c[3], _inv_c[1], _inv_c[5]],
                    [_inv_c[4], _inv_c[5], _inv_c[2]],
                ])[:, :, 0]

                _val = self._get_single_stress(
                    tr_c[cell, qp, 0, 0], det_f[cell, qp, 0, 0]**2, _c_inv,
                    coef[cell, qp, 0, 0], exp[cell, qp, 0, 0])
                out[cell, qp, :, 0] = [
                    _val[0, 0], _val[1, 1], _val[2, 2],
                    _val[0, 1], _val[0, 2], _val[1, 2]]
        return out

    @staticmethod
    def _get_single_tan_mod(i_1, i_3, c_inv, coef, exp):
        krond = nm.eye(3)
        _bracket = i_3**(-1/3) * i_1 - 3
        _bracket_p1 = _bracket**(exp - 1) \
            if (_bracket > 0. or exp >= 1.) else 1.
        _bracket_p2 = _bracket**(exp - 2) \
            if (_bracket > 0. or exp >=  2.) else 0.

        tan_mod = nm.zeros((3, 3, 3, 3))
        for ii, jj, kk, ll in zip(*[ind.flatten()
                                    for ind in nm.indices(tan_mod.shape)]):
            tan_mod[ii, jj, kk, ll] = 4 / 3 * exp * coef * (
                3 * (exp - 1) * _bracket_p2 * i_3**(-2/3)
                *krond[ii, jj] * krond[kk, ll]
                -(
                    (exp - 1) * _bracket_p2 * i_1 * i_3**(-2/3)
                    +_bracket_p1 * i_3**(-1/3)
                ) * (krond[ii, jj] * c_inv[kk, ll]
                     +c_inv[ii, jj] * krond[kk, ll])
                +(
                    (exp - 1) * _bracket_p2 * i_1**2 * i_3**(-2/3)
                    +_bracket_p1 * i_3**(-1/3) * i_1
                ) / 3 * c_inv[ii, jj] * c_inv[kk, ll]
                +.5 * _bracket_p1 * i_3**(-1/3) * i_1 * (
                    c_inv[ii, kk] * c_inv[jj, ll]
                    +c_inv[ii, ll] * c_inv[jj, kk])
            )

        return tan_mod

    def tan_mod_function(self, out, mat, *fargs, **kwargs):
        det_f, tr_c, inv_c = fargs
        mat = HyperElasticBase.tile_mat(mat, det_f.shape[0])
        coef, exp = mat[:, :, :, :1], mat[:, :, :, 1:]

        n_cells, n_qps, _, _ = out.shape
        for cell in range(n_cells):
            for qp in range(n_qps):
                _inv_c = inv_c[cell, qp, :, 0]
                _c_inv = nm.array([
                    [_inv_c[0], _inv_c[3], _inv_c[4]],
                    [_inv_c[3], _inv_c[1], _inv_c[5]],
                    [_inv_c[4], _inv_c[5], _inv_c[2]],
                ])

                _dh = self._get_single_tan_mod(
                    tr_c[cell, qp, 0, 0], det_f[cell, qp, 0, 0]**2, _c_inv,
                    coef[cell, qp, 0, 0], exp[cell, qp, 0, 0],
                )

                out[cell, qp] = nm.array([
                    [_dh[0, 0, 0, 0], _dh[0, 0, 1, 1], _dh[0, 0, 2, 2],
                     _dh[0, 0, 0, 1], _dh[0, 0, 0, 2], _dh[0, 0, 1, 2]],
                    [_dh[1, 1, 0, 0], _dh[1, 1, 1, 1], _dh[1, 1, 2, 2],
                     _dh[1, 1, 0, 1], _dh[1, 1, 0, 2], _dh[1, 1, 1, 2]],
                    [_dh[2, 2, 0, 0], _dh[2, 2, 1, 1], _dh[2, 2, 2, 2],
                     _dh[2, 2, 0, 1], _dh[2, 2, 0, 2], _dh[2, 2, 1, 2]],
                    [_dh[0, 1, 0, 0], _dh[0, 1, 1, 1], _dh[0, 1, 2, 2],
                     _dh[0, 1, 0, 1], _dh[0, 1, 0, 2], _dh[0, 1, 1, 2]],
                    [_dh[0, 2, 0, 0], _dh[0, 2, 1, 1], _dh[0, 2, 2, 2],
                     _dh[0, 2, 0, 1], _dh[0, 2, 0, 2], _dh[0, 2, 1, 2]],
                    [_dh[1, 2, 0, 0], _dh[1, 2, 1, 1], _dh[1, 2, 2, 2],
                     _dh[1, 2, 0, 1], _dh[1, 2, 0, 2], _dh[1, 2, 1, 2]],
                ])

class OgdenTLTerm(HyperElasticTLBase):
    r"""
    Single term of the hyperelastic Ogden model [1] with the strain energy
    density

    .. math::
        W = \frac{\mu}{\alpha} \, \left(
            \lambda_1^{\alpha} + \lambda_2^{\alpha} + \lambda_3^{\alpha}
            - 3 \right) \; ,

    where :math:`\lambda_k, k=1, 2, 3` are the principal stretches, whose
    squares are the principal values of the right Cauchy-Green deformation
    tensor :math:`\mathbf{C}`.

    Effective stress (2nd Piola-Kirchhoff) is [2]

    .. math::
        S_{ij} = 2 \, \frac{\partial W}{\partial C_{ij}} =
        \sum_{k=1}^3 S^{(k)} \, N^{(k)}_i \, N^{(k)}_j \; ,

    where the principal stresses are

    .. math::
        S^{(k)} = J^{-2/3} \, \left(
            \mu \, \bar\lambda^{\alpha - 2}
            -\sum_{j=1}^3 \frac{\mu}{3}
            \frac{\lambda_j^{\alpha}}{\lambda_k^2} \right) \; ,
        \quad k = 1, 2, 3 \; .

    and :math:`\mathbf{N}^{(k)}`, :math:`k=1, 2, 3` are the eigenvectors of
    :math:`\mathbf{C}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`p, K`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`

    [1] Ogden, R. W. Large deformation isotropic elasticity - on the
    correlation of theory and experiment for incompressible rubberlike solids.
    Proceedings of the Royal Society A, Vol. 326, No. 1567, Pp. 565-584 (1972),
    DOI `10.1098/rspa.1972.0026 <https://doi.org/10.1098/rspa.1972.0026>`_.

    [2] Steinmann, P., Hossain, M., Possart, G. Hyperelastic models for
    rubber-like materials: Consistent tangent operators and suitability for
    Treloar's data. Archive of Applied Mechanics, Vol. 82, No. 9, Pp. 1183-1217
    (2012), DOI `10.1007/s00419-012-0610-z
    <https://dx.doi.org/10.1007/s00419-012-0610-z>`_.
    """
    name = 'dw_tl_he_ogden'
    family_data_names = ['det_f', 'sym_c', 'tr_c', 'sym_inv_c']
    arg_shapes = {'material' : '1, 2',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['3_4', '3_8']

    @staticmethod
    def _get_single_stress(lbds, nks, det_f, coef, exp):
        a_p = sum(lbds**exp)
        s_k = [
            coef * det_f**(-2. / 3) *(
                lbdi**(exp - 2) - a_p / lbdi**2 / 3)
            for lbdi, ni in zip(lbds, nks.T)]
        out = nm.sum(
            [ski * nm.outer(nki, nki) for ski, nki in zip(s_k, nks.T)], axis=0)
        return out

    def stress_function(self, out, mat, *fargs, **kwargs):
        det_f, sym_c, _, _ = fargs
        mat = HyperElasticBase.tile_mat(mat, det_f.shape[0])
        coef, exp = mat[:, :, :, :1], mat[:, :, :, 1:]

        # compute principal stretches and directions
        c_mats = sym_c[:, :, [[0, 3, 4], [3, 1, 5], [4, 5, 2]], 0]
        lbds, nks = nm.linalg.eigh(c_mats)
        lbds = lbds**.5

        # evaluate stress
        n_cells, n_qps, _, _ = out.shape
        for cell in range(n_cells):
            for qp in range(n_qps):
                _val = self._get_single_stress(
                    lbds[cell, qp], nks[cell, qp], det_f[cell, qp, 0, 0],
                    coef[cell, qp, 0, 0], exp[cell, qp, 0, 0])
                out[cell, qp, :, 0] = [
                    _val[0, 0], _val[1, 1], _val[2, 2],
                    _val[0, 1], _val[0, 2], _val[1, 2]]
        return out

    @staticmethod
    def _get_single_tan_mod(total_lbds, nks, det_f, coef, exp):
        lbds = det_f**(-1/3) * total_lbds
        _bracket = [
            (det_f**(-1 / 3) * lbdi)**(exp - 2)
            -sum([1. / 3 * lbdj**exp / lbdi**2 for lbdj in lbds])
            for lbdi in lbds]
        s_k = nm.array([
            coef * det_f**(-2 / 3) * _bracket_ii
            for lbdi, _bracket_ii in zip(lbds, _bracket)])

        dj_dlbd = [det_f / lbdj for lbdj in lbds]

        dlbd_dlbd = det_f**(-1/3) * nm.array([[
            nm.eye(3)[ii, jj] - lbds[ii] / lbds[jj] / 3
            for jj in range(3)] for ii, dj_dlbdi in enumerate(dj_dlbd)])

        dsk_dlbd_1 = coef * nm.array([[
            -2 / 3 * det_f**(-5 / 3) * dj_dlbd[ii] * lbds[kk]**(exp - 2)
            +det_f**(-2/3) * (exp - 2) * lbds[kk]**(exp - 3) * dlbd_dlbd[kk, ii]
            for ii in range(3)] for kk in range(3)])

        dsk_dlbd_2 = coef / 3 * nm.array([[
            -2 / 3 * det_f**(-5 / 3) * dj_dlbd[ii] * nm.sum([
                lbds[jj]**exp / lbds[kk]**2
                for jj in range(3)])
            +det_f**(-2/3) * nm.sum([
                exp * lbds[jj]**(exp - 1) / lbds[kk]**2 * dlbd_dlbd[jj, ii]
                -2 * lbds[jj]**exp / lbds[kk]**3 * dlbd_dlbd[kk, ii]
                for jj in range(3)])
            for ii in range(3)] for kk in range(3)])
        dsk_dlbd = dsk_dlbd_1 - dsk_dlbd_2

        tan_mod = nm.zeros((3, 3, 3, 3))
        for mm, nn, pp, qq in zip(*[ind.flatten()
                                    for ind in nm.indices(tan_mod.shape)]):
            for ii, jj in zip(*[ind.flatten() for ind in nm.indices((3, 3))]):
                tan_mod[mm, nn, pp, qq] += 1. / lbds[jj] * dsk_dlbd[ii, jj] * (
                    nks[ii, mm] * nks[ii, nn] * nks[jj, pp] * nks[jj, qq])
                if ii != jj:
                    if lbds[ii] != lbds[jj]:
                        tan_mod[mm, nn, pp, qq] += (s_k[jj] - s_k[ii]) / \
                            (lbds[jj]**2 - lbds[ii]**2) * (
                                nks[ii, mm] * nks[jj, nn] * nks[ii, pp]
                                *nks[jj, qq]
                                +nks[ii, mm] * nks[jj, nn] * nks[jj, pp]
                                *nks[ii, qq])
                    else:
                        _val = 0.5 * (
                            dsk_dlbd[ii, ii] - dsk_dlbd[jj, ii]) / lbds[ii]
                        tan_mod[mm, nn, pp, qq] += _val * (
                                nks[ii, mm] * nks[jj, nn] * nks[ii, pp]
                                *nks[jj, qq]
                                +nks[ii, mm] * nks[jj, nn] * nks[jj, pp]
                                *nks[ii, qq])
        return tan_mod

    def tan_mod_function(self, out, mat, *fargs, **kwargs):
        det_f, sym_c, tr_c, inv_c = fargs
        mat = HyperElasticBase.tile_mat(mat, det_f.shape[0])
        coef, exp = mat[:, :, :, :1], mat[:, :, :, 1:]

        # compute principal stretches and directions
        c_mats = sym_c[:, :, [[0, 3, 4], [3, 1, 5], [4, 5, 2]], 0]
        lbds, nks = nm.linalg.eigh(c_mats)
        lbds = lbds**.5

        n_cells, n_qps, _, _ = out.shape
        for cell in range(n_cells):
            for qp in range(n_qps):
                _dh = self._get_single_tan_mod(
                    lbds[cell, qp], nks[cell, qp].T, det_f[cell, qp, 0, 0],
                    coef[cell, qp, 0, 0], exp[cell, qp, 0, 0],
                )

                out[cell, qp] = nm.array([
                    [_dh[0, 0, 0, 0], _dh[0, 0, 1, 1], _dh[0, 0, 2, 2],
                     _dh[0, 0, 0, 1], _dh[0, 0, 0, 2], _dh[0, 0, 1, 2]],
                    [_dh[1, 1, 0, 0], _dh[1, 1, 1, 1], _dh[1, 1, 2, 2],
                     _dh[1, 1, 0, 1], _dh[1, 1, 0, 2], _dh[1, 1, 1, 2]],
                    [_dh[2, 2, 0, 0], _dh[2, 2, 1, 1], _dh[2, 2, 2, 2],
                     _dh[2, 2, 0, 1], _dh[2, 2, 0, 2], _dh[2, 2, 1, 2]],
                    [_dh[0, 1, 0, 0], _dh[0, 1, 1, 1], _dh[0, 1, 2, 2],
                     _dh[0, 1, 0, 1], _dh[0, 1, 0, 2], _dh[0, 1, 1, 2]],
                    [_dh[0, 2, 0, 0], _dh[0, 2, 1, 1], _dh[0, 2, 2, 2],
                     _dh[0, 2, 0, 1], _dh[0, 2, 0, 2], _dh[0, 2, 1, 2]],
                    [_dh[1, 2, 0, 0], _dh[1, 2, 1, 1], _dh[1, 2, 2, 2],
                     _dh[1, 2, 0, 1], _dh[1, 2, 0, 2], _dh[1, 2, 1, 2]],
                ])

class MooneyRivlinTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic Mooney-Rivlin term. Effective stress
    :math:`S_{ij} = \kappa J^{-\frac{4}{3}} (C_{kk} \delta_{ij} - C_{ij}
    - \frac{2}{3 } I_2 C_{ij}^{-1})`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`\kappa`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_mooney_rivlin'
    family_data_names = ['det_f', 'tr_c', 'sym_inv_c', 'sym_c', 'in2_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_mooney_rivlin)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_mooney_rivlin)

class BulkPenaltyTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic bulk penalty term. Stress
    :math:`S_{ij} = K(J-1)\; J C_{ij}^{-1}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`K`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """

    name = 'dw_tl_bulk_penalty'
    family_data_names = ['det_f', 'sym_inv_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_bulk)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_bulk)

class BulkActiveTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic bulk active term. Stress :math:`S_{ij} = A J C_{ij}^{-1}`,
    where :math:`A` is the activation in :math:`[0, F_{\rm max}]`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`A`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """

    name = 'dw_tl_bulk_active'
    family_data_names = ['det_f', 'sym_inv_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_bulk_active)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_bulk_active)

class BulkPressureTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic bulk pressure term. Stress
    :math:`S_{ij} = -p J C_{ij}^{-1}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(p) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - virtual : :math:`\ul{v}`
        - state   : :math:`\ul{u}`
        - state_p : :math:`p`
    """
    name = 'dw_tl_bulk_pressure'
    arg_types = ('virtual', 'state', 'state_p')
    arg_geometry_types = {('state_p', None) : {'facet_extra' : 'facet'}}
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D', 'state_p' : 1}
    family_data_names = ['det_f', 'sym_inv_c']

    weak_function = staticmethod(terms.dw_he_rtm)
    weak_dp_function = staticmethod(terms.dw_tl_volume)

    stress_function = staticmethod(terms.dq_tl_stress_bulk_pressure)
    tan_mod_u_function = staticmethod(terms.dq_tl_tan_mod_bulk_pressure_u)

    def compute_data(self, family_data, mode, **kwargs):
        det_f, sym_inv_c = family_data.det_f, family_data.sym_inv_c
        p_qp = family_data.p_qp

        if mode == 0:
            out = nm.empty_like(sym_inv_c)
            fun = self.stress_function

        elif mode == 1:
            shape = list(sym_inv_c.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)
            fun = self.tan_mod_u_function

        else:
            raise ValueError('bad mode! (%d)' % mode)

        fun(out, p_qp, det_f, sym_inv_c)

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
                         stress, tan_mod, fd.mtx_f, fd.det_f, vgv, fmode, 0)

            else:
                vgs, _ = self.get_mapping(state_p)

                fargs =  (self.weak_dp_function,
                          fd.mtx_f, fd.sym_inv_c, fd.det_f, vgs, vgv, 1, -1)

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

class VolumeTLTerm(HyperElasticTLBase):
    r"""
    Volume term (weak form) in the total Lagrangian formulation.

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
    name = 'dw_tl_volume'
    arg_types = ('virtual', 'state')
    arg_geometry_types = {('virtual', None) : {'facet_extra' : 'facet'}}
    arg_shapes = {'virtual' : (1, None), 'state' : 'D'}
    family_data_names = ['mtx_f', 'det_f', 'sym_inv_c']

    function = staticmethod(terms.dw_tl_volume)

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

        elif (mode == 'eval') or (mode == 'el_avg'):
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

        return fd.mtx_f, fd.sym_inv_c, fd.det_f, vgs, vgv, 0, fmode

    def get_eval_shape(self, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

class DiffusionTLTerm(HyperElasticTLBase):
    r"""
    Diffusion term in the total Lagrangian formulation with
    linearized deformation-dependent permeability
    :math:`\ull{K}(\ul{u}) = J \ull{F}^{-1} \ull{k} f(J) \ull{F}^{-T}`,
    where :math:`\ul{u}` relates to the previous time step :math:`(n-1)`
    and
    :math:`f(J) = \max\left(0, \left(1 + \frac{(J - 1)}{N_f}\right)\right)^2`
    expresses the dependence on volume compression/expansion.

    :Definition:

    .. math::
        \int_{\Omega} \ull{K}(\ul{u}^{(n-1)}) : \pdiff{q}{\ul{X}}
        \pdiff{p}{\ul{X}}

    :Arguments:
        - material_1 : :math:`\ull{k}`
        - material_2 : :math:`N_f`
        - virtual    : :math:`q`
        - state      : :math:`p`
        - parameter  : :math:`\ul{u}^{(n-1)}`
    """
    name = 'dw_tl_diffusion'
    arg_types = ('material_1', 'material_2', 'virtual', 'state', 'parameter')
    arg_shapes = {'material_1' : 'D, D', 'material_2' : '1, 1',
                  'virtual' : (1, 'state'), 'state' : 1, 'parameter' : 'D'}
    family_data_names = ['mtx_f', 'det_f']

    function = staticmethod(terms.dw_tl_diffusion)

    def get_fargs(self, perm, ref_porosity, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgv, _ = self.get_mapping(parameter)

        name = parameter.name
        fd = self.get_family_data(parameter, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])
        grad = self.get(state, 'grad')

        if mode == 'weak':
            if diff_var is None:
                fmode = 0

            else:
                fmode = 1

        elif mode == 'el_avg':
            if term_mode == 'diffusion_velocity':
                fmode = 2

            else:
                raise ValueError('unsupported term evaluation mode in %s! (%s)'
                                 % (self.name, term_mode))

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

        return grad, perm, ref_porosity, fd.mtx_f, fd.det_f, vgv, fmode

    def get_eval_shape(self, perm, ref_porosity, virtual, state, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, dim, 1), state.dtype

class HyperElasticSurfaceTLFamilyData(HyperElasticFamilyData):
    """
    Family data for TL formulation applicable for surface terms.
    """
    family_function = staticmethod(terms.dq_tl_finite_strain_surface)
    cache_name = 'tl_surface_common'
    data_names = ('mtx_f', 'det_f', 'inv_f')

    def __call__(self, state, region, integral, geometry_type,
                 step=0, derivative=None):
        sg, _ = state.field.get_mapping(region,
                                        integral, geometry_type[0],
                                        get_saved=True)
        sd = state.field.extra_data[f'sd_{region.name}']

        vec = state(step=step, derivative=derivative)

        st_shape = state.get_data_shape(integral, geometry_type[0],
                                        region.name)
        data = self.init_data_struct(st_shape, name='surface_family_data')

        fargs = tuple([getattr(data, k) for k in self.data_names])
        fargs = fargs + (vec, sg, sd.fis, state.field.econn)
        fargs = Term.translate_fargs_mapping(self.family_function,
                                             list(fargs))

        self.family_function(*fargs)

        return data

class HyperElasticSurfaceTLBase(HyperElasticTLBase):
    """
    Base class for all hyperelastic surface terms in TL formulation family.
    """
    get_family_data = HyperElasticSurfaceTLFamilyData()

class SurfaceFluxTLTerm(HyperElasticSurfaceTLBase):
    r"""
    Surface flux term in the total Lagrangian formulation, consistent with
    :class:`DiffusionTLTerm`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{\nu} \cdot \ull{K}(\ul{u}^{(n-1)}) \pdiff{p}{\ul{X}}

    :Arguments:
        - material_1 : :math:`\ull{k}`
        - material_2 : :math:`N_f`
        - parameter_1 : :math:`p`
        - parameter_2 : :math:`\ul{u}^{(n-1)}`
    """
    name = 'ev_tl_surface_flux'
    arg_types = ('material_1', 'material_2', 'parameter_1', 'parameter_2')
    arg_shapes = {'material_1' : 'D, D', 'material_2' : '1, 1',
                  'parameter_1' : 1, 'parameter_2' : 'D'}
    family_data_names = ['det_f', 'inv_f']
    integration = 'facet_extra'

    function = staticmethod(terms.d_tl_surface_flux)

    def get_fargs(self, perm, ref_porosity, pressure, displacement,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(displacement)
        name = displacement.name
        fd = self.get_family_data(displacement, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])
        grad = self.get(pressure, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 0)

        return grad, perm, ref_porosity, fd.inv_f, fd.det_f, sg, fmode

    def get_eval_shape(self, perm, ref_porosity, pressure, displacement,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_en, n_c = self.get_data_shape(displacement)

        return (n_fa, 1, 1, 1), pressure.dtype

class SurfaceTractionTLTerm(HyperElasticSurfaceTLBase):
    r"""
    Surface traction term in the total Lagrangian formulation, expressed
    using :math:`\ul{\nu}`, the outward unit normal vector w.r.t. the
    undeformed surface, :math:`\ull{F}(\ul{u})`, the deformation gradient,
    :math:`J = \det(\ull{F})`, and :math:`\ull{\sigma}` a given traction,
    often equal to a given pressure, i.e.
    :math:`\ull{\sigma} = \pi \ull{I}`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{\nu} \cdot \ull{F}^{-1} \cdot \ull{\sigma} \cdot
        \ul{v} J

    :Arguments:
        - material : :math:`\ull{\sigma}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_surface_traction'
    arg_types = ('opt_material', 'virtual', 'state')
    arg_shapes = [{'opt_material' : 'D, D', 'virtual' : ('D', 'state'),
                   'state' : 'D'},
                  {'opt_material' : None}]
    family_data_names = ['det_f', 'inv_f']
    integration = 'facet_extra'

    function = staticmethod(terms.dw_tl_surface_traction)

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)
        sd = virtual.field.extra_data[f'sd_{self.region.name}']
        bf = virtual.field.eval_basis(sd.bkey, 0, self.integral)

        name = state.name
        fd = self.get_family_data(state, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])

        if mat is None:
            eye = nm.eye(sg.dim, dtype=nm.float64)
            mat = nm.tile(eye, ((1, sg.n_qp, 1, 1)))

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return mat, fd.det_f, fd.inv_f, bf, sg, sd.fis, fmode

class VolumeSurfaceTLTerm(HyperElasticSurfaceTLBase):
    r"""
    Volume of a :math:`D`-dimensional domain, using a surface integral in the
    total Lagrangian formulation, expressed using :math:`\ul{\nu}`, the outward
    unit normal vector w.r.t. the undeformed surface, :math:`\ull{F}(\ul{u})`,
    the deformation gradient, and :math:`J = \det(\ull{F})`. Uses the
    approximation of :math:`\ul{u}` for the deformed surface coordinates
    :math:`\ul{x}`.

    :Definition:

    .. math::
        1 / D \int_{\Gamma} \ul{\nu} \cdot \ull{F}^{-1} \cdot \ul{x} J

    :Arguments:
        - parameter : :math:`\ul{u}`
    """
    name = 'ev_tl_volume_surface'
    arg_types = ('parameter',)
    arg_shapes = {'parameter' : 'D'}
    family_data_names = ['det_f', 'inv_f']
    integration = 'facet_extra'

    function = staticmethod(terms.d_tl_volume_surface)

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)
        sd = parameter.field.extra_data[f'sd_{self.region.name}']
        bf = parameter.field.eval_basis(sd.bkey, 0, self.integral)

        name = parameter.name
        fd = self.get_family_data(parameter, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])

        asc = nm.ascontiguousarray

        coors0 = parameter.field.get_coor()
        coors = asc(coors0 + parameter().reshape(coors0.shape))

        return coors, fd.det_f, fd.inv_f, bf, sg, asc(sd.econn)

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, 1, 1), parameter.dtype
