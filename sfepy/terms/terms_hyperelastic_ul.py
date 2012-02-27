import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import terms
from sfepy.terms.terms_hyperelastic_base import HyperElasticBase

_msg_missing_data = 'missing family data!'

class HyperElasticULBase(HyperElasticBase):
    """
    Base class for all hyperelastic terms in UL formulation family.

    The subclasses should have the following static method attributes:
    - `stress_function()` (the stress)
    - `tan_mod_function()` (the tangent modulus)

    The common (family) data are cached in the evaluate cache of state
    variable.
    """
    family_function = staticmethod(terms.dq_finite_strain_ul)
    weak_function = staticmethod(terms.dw_he_rtm)
    fd_cache_name = 'ul_common'
    hyperelastic_mode = 1

    def compute_family_data(self, state):
        ap, vg = self.get_approximation(state, get_saved=True)

        vec = self.get_vector(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim * (dim + 1) / 2

        shapes = {
            'mtx_f' : (n_el, n_qp, dim, dim),
            'det_f' : (n_el, n_qp, 1, 1),
            'sym_b' : (n_el, n_qp, sym, 1),
            'tr_b' : (n_el, n_qp, 1, 1),
            'in2_b' : (n_el, n_qp, 1, 1),
            'green_strain' : (n_el, n_qp, sym, 1),
        }
        data = Struct(name='ul_family_data')
        for key, shape in shapes.iteritems():
            setattr(data, key, nm.zeros(shape, dtype=nm.float64))

        self.family_function(data.mtx_f,
                             data.det_f,
                             data.sym_b,
                             data.tr_b,
                             data.in2_b,
                             data.green_strain,
                             vec, vg, ap.econn)

        return data

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

# class NeoHookeanULHTerm(NeoHookeanULTerm):
#     r"""
#     Hyperelastic neo-Hookean term.  Geometrical configuration given by
#     parameter :math:`\ul{w}`.  Effective stress :math:`\tau_{ij} = \mu
#     J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

#     :Definition:

#     .. math::
#         \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

#     :Arguments 1:
#         - material : :math:`\mu`
#         - virtual  : :math:`\ul{v}`
#         - state    : :math:`\ul{u}`
#         - state_u  : :math:`\ul{w}`
#     """
#     name = 'dw_ul_he_neohook_h'
#     arg_types = ('material', 'virtual', 'state', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']]}

#     def __init__(self, *args, **kwargs):
#         HyperElasticULBase.__init__(self, *args, **kwargs)
#         self.call_mode = 1

# class NeoHookeanULEHTerm(NeoHookeanULTerm):
#     r"""
#     Hyperelastic neo-Hookean term.
#     Geometrical configuration given by parameter :math:`\ul{w}`.
#     Effective stress :math:`\tau_{ij} = \mu J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

#     :Definition:

#     .. math::
#         \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

#     :Arguments:
#         - material    : :math:`\mu`
#         - parameter_1 : :math:`\ul{v}`
#         - parameter_2 : :math:`\ul{u}`
#         - state_u     : :math:`\ul{w}`
#     """
#     name = 'd_ul_he_neohook_h'
#     arg_types = ('material', 'parameter_1', 'parameter_2', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']]}

#     def __init__(self, *args, **kwargs):
#         HyperElasticULBase.__init__(self, *args, **kwargs)
#         self.call_mode = 2

# class NeoHookeanULEvalTerm(Term):

#     name = 'de_ul_he_neohook'
#     arg_types = ('material', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']]}
#     function = {'stress': terms.dq_ul_he_stress_neohook,
#                 'element_contribution' : terms.de_he_rtm}

#     def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
#         mat, state_u = self.get_args( ['material', 'state_u'], **kwargs )
#         ap, vg = self.get_approximation(state_u)

#         dim = ap.dim[0]
#         sym = (dim + 1) * dim / 2
#         shape = (chunk_size, 1, sym, 1)

#         cache = self.get_cache('finite_strain_ul', 0)
#         detF, trB, B = cache(['detF', 'trB', 'B'], self, 0, state=state_u)

#         stress = nm.empty_like(B)
#         fun = self.function['stress']
#         fun(stress, mat, detF, trB, B)

#         fun = self.function['element_contribution']
#         for out, chunk in self.char_fun(chunk_size, shape):
#             status = fun(out, stress, detF, vg, chunk, 1)
#             out1 = nm.sum(out,0).reshape((sym,))

#             yield out1, chunk, status

# class BulkPenaltyULHTerm(BulkPenaltyULTerm):
#     r"""
#     Hyperelastic bulk penalty term.
#     Geometrical configuration given by parameter :math:`\ul{w}`.
#     Stress :math:`\tau_{ij} = K(J-1)\; J \delta_{ij}`.

#     :Definition:

#     .. math::
#         \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

#     :Arguments:
#         - material : :math:`K`
#         - virtual  : :math:`\ul{v}`
#         - state    : :math:`\ul{u}`
#         - state_u  : :math:`\ul{w}`
#     """
#     name = 'dw_ul_bulk_penalty_h'
#     arg_types = ('material', 'virtual', 'state', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']]}

#     def __init__(self, *args, **kwargs):
#         HyperElasticULBase.__init__(self, *args, **kwargs)
#         self.call_mode = 1

# class BulkPenaltyULEHTerm(BulkPenaltyULTerm):
#     r"""
#     Hyperelastic bulk penalty term.
#     Geometrical configuration given by parameter :math:`\ul{w}`.
#     Stress :math:`\tau_{ij} = K(J-1)\; J \delta_{ij}`.

#     :Definition:

#     .. math::
#         \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

#     :Arguments:
#         - material    : :math:`K`
#         - parameter_1 : :math:`\ul{v}`
#         - parameter_2 : :math:`\ul{u}`
#         - state_u  : :math:`\ul{w}`
#     """
#     name = 'd_ul_bulk_penalty_h'
#     arg_types = ('material', 'parameter_1', 'parameter_2', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']]}

#     def __init__(self, *args, **kwargs):
#         HyperElasticULBase.__init__(self, *args, **kwargs)
#         self.call_mode = 2

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

        fd = self.get_family_data(state, 'ul_common', self.family_data_names)
        fd.p_qp = self.get(state_p, 'val')

        if mode == 'weak':
            ig = self.char_fun.ig

            if diff_var != state_p.name:
                if diff_var is None:
                    stress = self.compute_data(fd, 0, **kwargs)
                    self.stress_cache[ig] = stress
                    tan_mod = nm.array([0], ndmin=4, dtype=nm.float64)

                    fmode = 0

                else:
                    stress = self.stress_cache[ig]
                    if stress is None:
                        stress = self.compute_data(fd, 0, **kwargs)

                    tan_mod = self.compute_data(fd, 1, **kwargs)
                    fmode = 1

                fargs = (self.weak_function,
                         stress, tan_mod, fd.mtx_f, fd.det_f, vgv, fmode, 1)

            else:
                vgs, _ = self.get_mapping(state_p)

                fargs =  (self.weak_dp_function,
                          -vgs.bf, fd.det_f, vgv, 1, 1)

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
        sym = dim * (dim + 1) / 2

        return (n_el, 1, sym, 1), state.dtype

# class BulkPressureULHTerm(BulkPressureULTerm):
#     r"""
#     Hyperelastic bulk pressure term. Stress
#     :math:`S_{ij} = -p J \delta_{ij}`.

#     :Definition:

#     .. math::
#         \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

#     :Arguments:
#         - virtual : :math:`\ul{v}`
#         - state   : :math:`\ul{u}`
#         - state_p : :math:`p`
#         - state_u : :math:`w`
#     """
#     name = 'dw_ul_bulk_pressure_h'
#     arg_types = ('virtual', 'state', 'state_p', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']],
#                   'state_in_volume_qp' : [['state_p']]}

#     def __call__(self, diff_var=None, chunk_size=None, **kwargs):
#         term_mode, = self.get_kwargs(['term_mode'], **kwargs)
#         virtual, state, state_p, state_u = self.get_args(**kwargs)
#         apv, vgv = self.get_approximation(virtual)
#         aps, vgs = self.get_approximation(state_p)

#         self.set_data_shape(apv, aps)
#         shape, mode = self.get_shape_grad(diff_var, chunk_size)

#         cache = self.get_cache('finite_strain_ul', 0)
#         family_data = cache(['detF', 'B'], self, 0, state=state_u)

#         ig = self.char_fun.ig

#         if term_mode is None:

#             if mode < 2:
#                 stress = self.crt_data.stress[ig]
#                 if stress is None:
#                     stress = self.compute_crt_data(family_data, 0, **kwargs)
#                     self.crt_data.stress[ig] = stress
#                 tan_mod = self.crt_data.tan_mod[ig]
#                 if tan_mod is None:
#                     tan_mod = self.compute_crt_data(family_data, 1, **kwargs)
#                     self.crt_data.tan_mod[ig] = tan_mod

#                 fun = self.function['element_contribution']

#                 mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

#                 if mode == 0:
#                     vec = self.get_vector(state)
#                     for out, chunk in self.char_fun(chunk_size, shape):
#                         out2 = nm.zeros(out.shape[:-1] + (out.shape[-2],),
#                                         dtype=nm.float64)
#                         status1 = fun(out2, stress, tan_mod,
#                                       mtxF, detF, vgv, chunk, 1, 1)
#                         status2 = terms.he_residuum_from_mtx(out, out2, vec, apv.econn, chunk)
#                         yield out, chunk, status1 or status2

#                 else:
#                     for out, chunk in self.char_fun(chunk_size, shape):
#                         status = fun(out, stress, tan_mod,
#                                      mtxF, detF, vgv, chunk, 1, 1)
#                         yield out, chunk, status

#             else:
#                 from sfepy.base.base import debug
#                 debug()
#                 # fun = self.function['element_contribution_dp']

#                 # mtxF, B, detF = cache(['F', 'B', 'detF'],
#                 #                       self, 0, state=state_u)

#                 # bf = aps.get_base('v', 0, self.integral)
#                 # for out, chunk in self.char_fun(chunk_size, shape):
#                 #     status = fun(out, bf, mtxF, B, detF, vgv, 1, chunk, 1)
#                 #     yield -out, chunk, status

#         elif term_mode == 'd_eval':
#             raise NotImplementedError

# class BulkPressureULEHTerm(BulkPressureULTerm):
#     r"""
#     Hyperelastic bulk pressure term. Stress
#     :math:`S_{ij} = -p J \delta_{ij}`.

#     :Definition:

#     .. math::
#         \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

#     :Arguments:
#         - virtual : :math:`\ul{v}`
#         - state   : :math:`\ul{u}`
#         - state_p : :math:`p`
#         - state_u : :math:`w`
#     """
#     name = 'd_ul_bulk_pressure_h'
#     arg_types = ('virtual', 'state', 'state_p', 'state_u')
#     use_caches = {'finite_strain_ul' : [['state_u']],
#                   'state_in_volume_qp' : [['state_p']]}

#     def __call__(self, diff_var=None, chunk_size=None, **kwargs):
#         term_mode, = self.get_kwargs(['term_mode'], **kwargs)
#         par1, par2, state_p, state_u = self.get_args(**kwargs)
#         apv, vgv = self.get_approximation(par1)
#         aps, vgs = self.get_approximation(state_p)

#         self.set_data_shape(apv, aps)
#         n_el, n_qp, dim, n_ep = self.data_shape_v
#         shape0 = (1, dim * n_ep, dim * n_ep)
#         shape = (chunk_size, 1, 1, 1)

#         cache = self.get_cache('finite_strain_ul', 0)
#         family_data = cache(['detF', 'B'], self, 0, state=state_u)

#         ig = self.char_fun.ig
#         p1 = self.get_vector(par1)
#         p2 = self.get_vector(par2)

#         stress = self.crt_data.stress[ig]
#         if stress is None:
#             stress = self.compute_crt_data(family_data, 0, **kwargs)
#             self.crt_data.stress[ig] = stress
#         tan_mod = self.crt_data.tan_mod[ig]
#         if tan_mod is None:
#             tan_mod = self.compute_crt_data(family_data, 1, **kwargs)
#             self.crt_data.tan_mod[ig] = tan_mod

#         fun = self.function['element_contribution']
#         mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

#         for out, chunk in self.char_fun( chunk_size, shape ):
#             out2 = nm.zeros((out.shape[0],) + shape0, dtype=nm.float64)
#             status1 = fun(out2, stress, tan_mod,
#                           mtxF, detF, vgv, chunk, 1, 1)
#             status2 = terms.he_eval_from_mtx(out, out2, p1, p2, apv.econn, chunk)
#             out0 = nm.sum(out)

#             yield out0, chunk, status1 or status2

# class BulkPressureULEvalTerm(Term):

#     name = 'de_ul_bulk_pressure'
#     arg_types = ('state_u', 'state_p')
#     use_caches = {'finite_strain_ul' : [['state_u']],
#                   'state_in_volume_qp' : [['state_p']]}

#     function = {'stress': terms.dq_ul_stress_bulk_pressure,
#                 'element_contribution' : terms.de_he_rtm}

#     def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
#         state_u, state_p = self.get_args( ['state_u', 'state_p'], **kwargs )
#         ap, vg = self.get_approximation(state_u)

#         dim = ap.dim[0]
#         sym = (dim + 1) * dim / 2
#         shape = (chunk_size, 1, sym, 1)

#         cache = self.get_cache( 'finite_strain_ul', 0 )
#         detF, B = cache(['detF', 'B'], self, 0, state=state_u)
#         cache = self.get_cache('state_in_volume_qp', 0)
#         p_qp = cache('state', self, 0, state=state_p, get_vector=self.get_vector)

#         stress = nm.empty_like(B)
#         fun = self.function['stress']
#         fun(stress, p_qp, detF)

#         fun = self.function['element_contribution']
#         for out, chunk in self.char_fun(chunk_size, shape):
#             status = fun(out, stress, detF, vg, chunk, 1)
#             out1 = nm.sum(out,0).reshape((sym,))

#             yield out1, chunk, status

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
    family_data_names = ['mtx_f', 'det_f']

    function = staticmethod(terms.dw_ul_volume)
    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgv, _ = self.get_mapping(virtual)
        vgs, _ = self.get_mapping(state)

        fd = self.get_family_data(state, 'ul_common', self.family_data_names)

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

        return vgv.bf, fd.det_f, vgs, 0, fmode

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
    family_data_names = ['mtx_f', 'det_f']

    function = staticmethod(terms.dw_mass_scalar)

    def get_fargs(self, bulk, virtual, state, parameter_u,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgp, _ = self.get_mapping(virtual)
        vgu, _ = self.get_mapping(parameter_u)

        fd = self.get_family_data(parameter_u, 'ul_common', self.family_data_names)

        coef = nm.divide(bulk, fd.det_f)

        if mode == 'weak':
            if diff_var is None:
                val_qp = self.get(state, 'val')
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return coef, val_qp, vgp.bf, vgp, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))
