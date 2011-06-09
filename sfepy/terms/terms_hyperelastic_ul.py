import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import terms, Term
from sfepy.terms.terms_hyperelastic_base \
     import CouplingVectorScalarHE, HyperElasticBase
from sfepy.terms.terms_base import VectorVector, ScalarScalar, InstantaneousBase

class HyperElasticULBase( HyperElasticBase ):
    """Base class for all hyperelastic terms in UL formulation. This is not a
    proper Term!
    """
    use_caches = {'finite_strain_ul' : [['state']]}
    mode = 'ul'

class NeoHookeanULTerm( VectorVector, HyperElasticULBase ):
    r"""
    :Description:
    Hyperelastic neo-Hookean term. Effective stress
    :math:`\tau_{ij} = \mu J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        material : :math:`\mu`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_ul_he_neohook'
    arg_types = ('material', 'virtual', 'state')

    family_data_names = ['detF', 'trB', 'B']
    term_function = {'stress' : terms.dq_ul_he_stress_neohook,
                     'tangent_modulus' : terms.dq_ul_he_tan_mod_neohook}

    def compute_crt_data( self, family_data, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        detF, trB, B = family_data

        if mode == 0:
            out = nm.empty_like( B )
            fun = self.term_function['stress']
        else:
            shape = list( B.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.term_function['tangent_modulus']

        fun( out, mat, detF, trB, B )

        return out

class NeoHookeanULHTerm(NeoHookeanULTerm):
    r"""
    :Description:
    Hyperelastic neo-Hookean term.
    Geometrical configuration given by parameter :math:`\ul{w}`.
    Effective stress :math:`\tau_{ij} = \mu J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments 1:
        material : :math:`\mu`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`,
        state_u  : :math:`\ul{w}`
    """
    name = 'dw_ul_he_neohook_h'
    arg_types = ('material', 'virtual', 'state', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']]}

    def __init__(self, *args, **kwargs):
        HyperElasticULBase.__init__(self, *args, **kwargs)
        self.call_mode = 1

class NeoHookeanULEHTerm(NeoHookeanULTerm):
    r"""
    :Description:
    Hyperelastic neo-Hookean term.
    Geometrical configuration given by parameter :math:`\ul{w}`.
    Effective stress :math:`\tau_{ij} = \mu J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        material    : :math:`\mu`,
        parameter_1 : :math:`\ul{v}`,
        parameter_2 : :math:`\ul{u}`,
        state_u     : :math:`\ul{w}`
    """
    name = 'd_ul_he_neohook_h'
    arg_types = ('material', 'parameter_1', 'parameter_2', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']]}

    def __init__(self, *args, **kwargs):
        HyperElasticULBase.__init__(self, *args, **kwargs)
        self.call_mode = 2

class NeoHookeanULEvalTerm(Term):

    name = 'de_ul_he_neohook'
    arg_types = ('material', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']]}
    function = {'stress': terms.dq_ul_he_stress_neohook,
                'element_contribution' : terms.de_he_rtm}

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, state_u = self.get_args( ['material', 'state_u'], **kwargs )
        ap, vg = self.get_approximation(state_u)

        dim = ap.dim[0]
        sym = (dim + 1) * dim / 2
        shape = (chunk_size, 1, sym, 1)

        cache = self.get_cache('finite_strain_ul', 0)
        detF, trB, B = cache(['detF', 'trB', 'B'], self, 0, state=state_u)

        stress = nm.empty_like(B)
        fun = self.function['stress']
        fun(stress, mat, detF, trB, B)

        fun = self.function['element_contribution']
        for out, chunk in self.char_fun(chunk_size, shape):
            status = fun(out, stress, detF, vg, chunk, 1)
            out1 = nm.sum(out,0).reshape((sym,))

            yield out1, chunk, status

class BulkPenaltyULTerm( VectorVector, HyperElasticULBase ):
    r"""
    :Description:
    Hyperelastic bulk penalty term. Stress
    :math:`\tau_{ij} = K(J-1)\; J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        material : :math:`K`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_ul_bulk_penalty'
    arg_types = ('material', 'virtual', 'state')

    family_data_names = ['detF', 'B']
    term_function = {'stress' : terms.dq_ul_he_stress_bulk,
                     'tangent_modulus' : terms.dq_ul_he_tan_mod_bulk}

    def compute_crt_data( self, family_data, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        detF, B = family_data

        if mode == 0:
            out = nm.empty_like( B )
            fun = self.term_function['stress']
        else:
            shape = list( B.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.term_function['tangent_modulus']

        fun( out, mat, detF )

        return out

class BulkPenaltyULHTerm(BulkPenaltyULTerm):
    r"""
    :Description:
    Hyperelastic bulk penalty term.
    Geometrical configuration given by parameter :math:`\ul{w}`.
    Stress :math:`\tau_{ij} = K(J-1)\; J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        material : :math:`K`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`,
        state_u  : :math:`\ul{w}`
    """
    name = 'dw_ul_bulk_penalty_h'
    arg_types = ('material', 'virtual', 'state', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']]}

    def __init__(self, *args, **kwargs):
        HyperElasticULBase.__init__(self, *args, **kwargs)
        self.call_mode = 1

class BulkPenaltyULEHTerm(BulkPenaltyULTerm):
    r"""
    :Description:
    Hyperelastic bulk penalty term.
    Geometrical configuration given by parameter :math:`\ul{w}`.
    Stress :math:`\tau_{ij} = K(J-1)\; J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        material    : :math:`K`,
        parameter_1 : :math:`\ul{v}`,
        parameter_2 : :math:`\ul{u}`,
        state_u  : :math:`\ul{w}`
    """
    name = 'd_ul_bulk_penalty_h'
    arg_types = ('material', 'parameter_1', 'parameter_2', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']]}

    def __init__(self, *args, **kwargs):
        HyperElasticULBase.__init__(self, *args, **kwargs)
        self.call_mode = 2

class MooneyRivlinULTerm( VectorVector, HyperElasticULBase ):
    r"""
    :Description:
    Hyperelastic Mooney-Rivlin term.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        material : :math:`\kappa`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_ul_he_mooney_rivlin'
    arg_types = ('material', 'virtual', 'state')

    family_data_names = ['detF', 'trB', 'B', 'in2B']
    term_function = {'stress' : terms.dq_ul_he_stress_mooney_rivlin,
                     'tangent_modulus' : terms.dq_ul_he_tan_mod_mooney_rivlin}

    def compute_crt_data( self, family_data, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        detF, trB, B, in2B = family_data

        if mode == 0:
            out = nm.empty_like( B )
            fun = self.term_function['stress']
        else:
            shape = list( B.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.term_function['tangent_modulus']

        fun( out, mat, detF, trB, B, in2B )

        return out

class BulkPressureULTerm(CouplingVectorScalarHE, HyperElasticULBase):
    r"""
    :Description:
    Hyperelastic bulk pressure term. Stress
    :math:`S_{ij} = -p J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{u}`,
        state_p : :math:`p`
    """

    name = 'dw_ul_bulk_pressure'
    arg_types = ('virtual', 'state', 'state_p')
    use_caches = {'finite_strain_ul' : [['state']],
                  'state_in_volume_qp' : [['state_p']]}

    term_function = {'stress' : terms.dq_ul_stress_bulk_pressure,
                     'tangent_modulus_u' : terms.dq_ul_tan_mod_bulk_pressure_u}

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.function = {
            'element_contribution' : terms.dw_he_rtm,
            'element_contribution_dp' : terms.dw_ul_volume,
        }
        igs = self.region.igs
        self.crt_data = Struct(stress={}.fromkeys(igs, None),
                               tan_mod={}.fromkeys(igs, None))

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        term_mode, = self.get_kwargs(['term_mode'], **kwargs)
        virtual, state, state_p = self.get_args(**kwargs)
        apv, vgv = self.get_approximation(virtual)
        aps, vgs = self.get_approximation(state_p)

        self.set_data_shape(apv, aps)
        shape, mode = self.get_shape_grad(diff_var, chunk_size)

        cache = self.get_cache('finite_strain_ul', 0)
        family_data = cache(['detF', 'B'], self, 0, state=state)

        if term_mode is None:

            if mode < 2:
                ig = self.char_fun.ig

                crt_data = self.compute_crt_data(family_data, mode, **kwargs)
                if mode == 0:
                    self.crt_data.stress[ig] = stress = crt_data
                    self.crt_data.tan_mod[ig] = nm.array([0], ndmin=4)
                else:
                    self.crt_data.tan_mod[ig] = crt_data

                    stress = self.crt_data.stress[ig]
                    if stress is None:
                        stress = self.compute_crt_data(family_data, 0, **kwargs)

                fun = self.function['element_contribution']

                mtxF, detF = cache(['F', 'detF'], self, 0, state=state)

                for out, chunk in self.char_fun(chunk_size, shape):
                    status = fun(out, stress, self.crt_data.tan_mod[ig], mtxF, detF,
                                 vgv, chunk, mode, 1)
                    yield out, chunk, status

            else:
                fun = self.function['element_contribution_dp']

                mtxF, B, detF = cache(['F', 'B', 'detF'],
                                      self, 0, state=state)

                bf = aps.get_base('v', 0, self.integral)
                for out, chunk in self.char_fun(chunk_size, shape):
                    status = fun(out, bf, detF, vgv, 1, chunk, 1)
                    yield -out, chunk, status

        elif term_mode == 'd_eval':
            raise NotImplementedError

        elif term_mode in ['strain', 'stress']:

            if term_mode == 'strain':
                out_qp = cache('E', self, 0, state=state)

            elif term_mode == 'stress':
                out_qp = self.compute_crt_data(family_data, 0, **kwargs)

            shape = (chunk_size, 1) + out_qp.shape[2:]
            for out, chunk in self.char_fun(chunk_size, shape):
                status = vgv.integrate_chunk(out, out_qp[chunk], chunk)
                out1 = out / vgv.variable(2)[chunk]

            yield out1, chunk, status

    def compute_crt_data(self, family_data, mode, **kwargs):
        detF, B = family_data

        p, = self.get_args(['state_p'], **kwargs)

        cache = self.get_cache('state_in_volume_qp', 0)
        p_qp = cache('state', self, 0, state=p, get_vector=self.get_vector)

        if mode == 0:
            out = nm.empty_like(B)
            fun = self.term_function['stress']
        elif mode == 1:
            shape = list(B.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)
            fun = self.term_function['tangent_modulus_u']
        else:
            raise ValueError('bad mode! (%d)' % mode)

        fun(out, p_qp, detF)

        return out

class BulkPressureULHTerm(BulkPressureULTerm):
    r"""
    :Description:
    Hyperelastic bulk pressure term. Stress
    :math:`S_{ij} = -p J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{u}`,
        state_p : :math:`p`,
        state_u : :math:`w`
    """
    name = 'dw_ul_bulk_pressure_h'
    arg_types = ('virtual', 'state', 'state_p', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']],
                  'state_in_volume_qp' : [['state_p']]}

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        term_mode, = self.get_kwargs(['term_mode'], **kwargs)
        virtual, state, state_p, state_u = self.get_args(**kwargs)
        apv, vgv = self.get_approximation(virtual)
        aps, vgs = self.get_approximation(state_p)

        self.set_data_shape(apv, aps)
        shape, mode = self.get_shape_grad(diff_var, chunk_size)

        cache = self.get_cache('finite_strain_ul', 0)
        family_data = cache(['detF', 'B'], self, 0, state=state_u)

        ig = self.char_fun.ig

        if term_mode is None:

            if mode < 2:
                stress = self.crt_data.stress[ig]
                if stress is None:
                    stress = self.compute_crt_data(family_data, 0, **kwargs)
                    self.crt_data.stress[ig] = stress
                tan_mod = self.crt_data.tan_mod[ig]
                if tan_mod is None:
                    tan_mod = self.compute_crt_data(family_data, 1, **kwargs)
                    self.crt_data.tan_mod[ig] = tan_mod

                fun = self.function['element_contribution']

                mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

                if mode == 0:
                    vec = self.get_vector(state)
                    for out, chunk in self.char_fun(chunk_size, shape):
                        out2 = nm.zeros(out.shape[:-1] + (out.shape[-2],),
                                        dtype=nm.float64)
                        status1 = fun(out2, stress, tan_mod,
                                      mtxF, detF, vgv, chunk, 1, 1)
                        status2 = terms.he_residuum_from_mtx(out, out2, vec, apv.econn, chunk)
                        yield out, chunk, status1 or status2

                else:
                    for out, chunk in self.char_fun(chunk_size, shape):
                        status = fun(out, stress, tan_mod,
                                     mtxF, detF, vgv, chunk, 1, 1)
                        yield out, chunk, status

            else:
                from sfepy.base.base import debug
                debug()
                # fun = self.function['element_contribution_dp']

                # mtxF, B, detF = cache(['F', 'B', 'detF'],
                #                       self, 0, state=state_u)

                # bf = aps.get_base('v', 0, self.integral)
                # for out, chunk in self.char_fun(chunk_size, shape):
                #     status = fun(out, bf, mtxF, B, detF, vgv, 1, chunk, 1)
                #     yield -out, chunk, status

        elif term_mode == 'd_eval':
            raise NotImplementedError

class BulkPressureULEHTerm(BulkPressureULTerm):
    r"""
    :Description:
    Hyperelastic bulk pressure term. Stress
    :math:`S_{ij} = -p J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J

    :Arguments:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{u}`,
        state_p : :math:`p`,
        state_u : :math:`w`
    """
    name = 'd_ul_bulk_pressure_h'
    arg_types = ('virtual', 'state', 'state_p', 'state_u')
    use_caches = {'finite_strain_ul' : [['state_u']],
                  'state_in_volume_qp' : [['state_p']]}

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        term_mode, = self.get_kwargs(['term_mode'], **kwargs)
        par1, par2, state_p, state_u = self.get_args(**kwargs)
        apv, vgv = self.get_approximation(par1)
        aps, vgs = self.get_approximation(state_p)

        self.set_data_shape(apv, aps)
        n_el, n_qp, dim, n_ep = self.data_shape_v
        shape0 = (1, dim * n_ep, dim * n_ep)
        shape = (chunk_size, 1, 1, 1)

        cache = self.get_cache('finite_strain_ul', 0)
        family_data = cache(['detF', 'B'], self, 0, state=state_u)

        ig = self.char_fun.ig
        p1 = self.get_vector(par1)
        p2 = self.get_vector(par2)

        stress = self.crt_data.stress[ig]
        if stress is None:
            stress = self.compute_crt_data(family_data, 0, **kwargs)
            self.crt_data.stress[ig] = stress
        tan_mod = self.crt_data.tan_mod[ig]
        if tan_mod is None:
            tan_mod = self.compute_crt_data(family_data, 1, **kwargs)
            self.crt_data.tan_mod[ig] = tan_mod

        fun = self.function['element_contribution']
        mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

        for out, chunk in self.char_fun( chunk_size, shape ):
            out2 = nm.zeros((out.shape[0],) + shape0, dtype=nm.float64)
            status1 = fun(out2, stress, tan_mod,
                          mtxF, detF, vgv, chunk, 1, 1)
            status2 = terms.he_eval_from_mtx(out, out2, p1, p2, apv.econn, chunk)
            out0 = nm.sum(out)

            yield out0, chunk, status1 or status2

class BulkPressureULEvalTerm(Term):

    name = 'de_ul_bulk_pressure'
    arg_types = ('state_u', 'state_p')
    use_caches = {'finite_strain_ul' : [['state_u']],
                  'state_in_volume_qp' : [['state_p']]}

    function = {'stress': terms.dq_ul_stress_bulk_pressure,
                'element_contribution' : terms.de_he_rtm}

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        state_u, state_p = self.get_args( ['state_u', 'state_p'], **kwargs )
        ap, vg = self.get_approximation(state_u)

        dim = ap.dim[0]
        sym = (dim + 1) * dim / 2
        shape = (chunk_size, 1, sym, 1)

        cache = self.get_cache( 'finite_strain_ul', 0 )
        detF, B = cache(['detF', 'B'], self, 0, state=state_u)
        cache = self.get_cache('state_in_volume_qp', 0)
        p_qp = cache('state', self, 0, state=state_p, get_vector=self.get_vector)

        stress = nm.empty_like(B)
        fun = self.function['stress']
        fun(stress, p_qp, detF)

        fun = self.function['element_contribution']
        for out, chunk in self.char_fun(chunk_size, shape):
            status = fun(out, stress, detF, vg, chunk, 1)
            out1 = nm.sum(out,0).reshape((sym,))

            yield out1, chunk, status

class VolumeULTerm(CouplingVectorScalarHE, InstantaneousBase, Term):
    r"""
    :Description:
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
        virtual : :math:`q`,
        state   : :math:`\ul{u}`
    """
    name = 'dw_ul_volume'
    arg_types = ('virtual', 'state')
    use_caches = {'finite_strain_ul' : [['state',
                                         {'F' : (2, 2),
                                          'detF' : (2, 2)}]]}

    function = staticmethod(terms.dw_ul_volume)

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        virtual, state = self.get_args( **kwargs )
        term_mode = kwargs.get('term_mode')

        apv, vgv = self.get_approximation(state)
        aps, vgs = self.get_approximation(virtual)

        self.set_data_shape(apv, aps)

        cache = self.get_cache('finite_strain_ul', 0)
        ih = self.arg_steps[state.name] # issue 104!
        mtxF, detF = cache(['F', 'detF'], self, ih, state=state)

        if term_mode == 'volume':
            n_el, _, _, _ = self.data_shape_s
            shape, mode = (n_el, 1, 1, 1), 2

        elif term_mode == 'rel_volume':
            n_el, _, _, _ = self.data_shape_s
            shape, mode = (n_el, 1, 1, 1), 3

        else:
            shape, mode = self.get_shape_div(diff_var, chunk_size)

        bf = aps.get_base('v', 0, self.integral)

        return (bf, detF, vgv, 0), shape, mode

class CompressibilityULTerm(ScalarScalar, Term):
    r"""
    :Description:
    Compressibility term in the updated Lagrangian formulation

    :Definition:
    .. math::
        \int_{\Omega} 1\over \gamma p \, q

    :Arguments:
        material: :math:`\gamma`,
        virtual    : :math:`q`,
        state      : :math:`p`,
    """
    name = 'dw_ul_compressible'
    arg_types = ('material', 'virtual', 'state', 'state_u')
    use_caches = {'finite_strain_ul': [['state_u']]}

    function = staticmethod(terms.dw_mass_scalar)

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        bulk, virtual, state, state_u = self.get_args(**kwargs)

        ap, vg = self.get_approximation(virtual)

        self.set_data_shape(ap)
        shape, mode = self.get_shape( diff_var, chunk_size )

        cache = self.get_cache('finite_strain_ul', 0)
        mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

        coef = nm.divide(bulk, detF)
        bf = ap.get_base('v', 0, self.integral)

        fargs = coef, state(), bf, vg, ap.econn

        return fargs, shape, mode
