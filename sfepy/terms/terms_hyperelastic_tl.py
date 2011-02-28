import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_hyperelastic_base \
     import CouplingVectorScalarTL, HyperElasticBase
from sfepy.terms.terms_base import VectorVector, ScalarScalar, InstantaneousBase
            
class HyperElasticTLBase( HyperElasticBase ):
    """Base class for all hyperelastic terms in TL formulation. This is not a
    proper Term!
    """
    use_caches = {'finite_strain_tl' : [['state']]}
    mode = 'tl'

class NeoHookeanTLTerm( VectorVector, HyperElasticTLBase ):
    r"""
    :Description:
    Hyperelastic neo-Hookean term. Effective stress
    :math:`S_{ij} = \mu J^{-\frac{2}{3}}(\delta_{ij} -
    \frac{1}{3}C_{kk}C_{ij}^{-1})`.

    :Definition:
    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        material : :math:`\mu`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_neohook'
    arg_types = ('material', 'virtual', 'state')

    family_data_names = ['detF', 'trC', 'invC']
    term_function = {'stress' : terms.dq_tl_he_stress_neohook,
                     'tangent_modulus' : terms.dq_tl_he_tan_mod_neohook}
    
    def compute_crt_data( self, family_data, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        detF, trC, invC = family_data

        if mode == 0:
            out = nm.empty_like( invC )
            fun = self.term_function['stress']
        else:
            shape = list( invC.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.term_function['tangent_modulus']

        fun( out, mat, detF, trC, invC )

        return out

class MooneyRivlinTLTerm( VectorVector, HyperElasticTLBase ):
    r"""
    :Description:
    Hyperelastic Mooney-Rivlin term. Effective stress
    :math:`S_{ij} = \kappa J^{-\frac{4}{3}} (C_{kk} \delta_{ij} - C_{ij}
    - \frac{2}{3 } I_2 C_{ij}^{-1})`.

    :Definition:
    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        material : :math:`\kappa`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_mooney_rivlin'
    arg_types = ('material', 'virtual', 'state')

    family_data_names = ['detF', 'trC', 'invC', 'C', 'in2C']
    term_function = {'stress' : terms.dq_tl_he_stress_mooney_rivlin,
                     'tangent_modulus' : terms.dq_tl_he_tan_mod_mooney_rivlin}

    def compute_crt_data( self, family_data, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        detF, trC, invC, vecC, in2C = family_data

        if mode == 0:
            out = nm.empty_like( invC )
            fun = self.term_function['stress']
        else:
            shape = list( invC.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.term_function['tangent_modulus']

        fun( out, mat, detF, trC, invC, vecC, in2C )

        return out

class BulkPenaltyTLTerm( VectorVector, HyperElasticTLBase ):
    r"""
    :Description:
    Hyperelastic bulk penalty term. Stress
    :math:`S_{ij} = K(J-1)\; J C_{ij}^{-1}`.

    :Definition:
    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        material : :math:`K`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """

    name = 'dw_tl_bulk_penalty'
    arg_types = ('material', 'virtual', 'state')

    family_data_names = ['detF', 'invC']
    term_function = {'stress' : terms.dq_tl_he_stress_bulk,
                     'tangent_modulus' : terms.dq_tl_he_tan_mod_bulk }

    def compute_crt_data( self, family_data, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        detF, invC = family_data
        
        if mode == 0:
            out = nm.empty_like( invC )
            fun = self.term_function['stress']
        else:
            shape = list( invC.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.term_function['tangent_modulus']

        fun( out, mat, detF, invC )

        return out

class BulkPressureTLTerm(CouplingVectorScalarTL, HyperElasticTLBase):
    r"""
    :Description:
    Hyperelastic bulk pressure term. Stress
    :math:`S_{ij} = -p J C_{ij}^{-1}`.

    :Definition:
    .. math::
        \int_{\Omega} S_{ij}(p) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{u}`,
        state_p : :math:`p`
    """

    name = 'dw_tl_bulk_pressure'
    arg_types = ('virtual', 'state', 'state_p')
    use_caches = {'finite_strain_tl' : [['state']],
                  'state_in_volume_qp' : [['state_p']]}

    term_function = {'stress' : terms.dq_tl_stress_bulk_pressure,
                     'tangent_modulus_u' : terms.dq_tl_tan_mod_bulk_pressure_u}

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.function = {
            'element_contribution' : terms.dw_he_rtm,
            'element_contribution_dp' : terms.dw_tl_volume,
        }
        igs = self.region.igs
        dummy = nm.array([0], ndmin=4)
        self.crt_data = Struct(stress={}.fromkeys(igs, None),
                               tan_mod=dummy)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        term_mode, = self.get_kwargs(['term_mode'], **kwargs)
        virtual, state, state_p = self.get_args(**kwargs)
        apv, vgv = self.get_approximation(virtual)
        aps, vgs = self.get_approximation(state_p)

        self.set_data_shape(apv, aps)
        shape, mode = self.get_shape_grad(diff_var, chunk_size)

        cache = self.get_cache('finite_strain_tl', 0)
        family_data = cache(['detF', 'invC'], self, 0, state=state)

        if term_mode is None:

            if mode < 2:
                ig = self.char_fun.ig

                crt_data = self.compute_crt_data(family_data, mode, **kwargs)
                if mode == 0:
                    self.crt_data.stress[ig] = stress = crt_data

                else:
                    self.crt_data.tan_mod = crt_data

                    stress = self.crt_data.stress[ig]
                    if stress is None:
                        stress = self.compute_crt_data(family_data, 0, **kwargs)

                fun = self.function['element_contribution']

                mtxF, detF = cache(['F', 'detF'], self, 0, state=state)

                for out, chunk in self.char_fun(chunk_size, shape):
                    status = fun(out, stress, self.crt_data.tan_mod, mtxF, detF,
                                 vgv, chunk, mode, 0)
                    yield out, chunk, status
            else:
                fun = self.function['element_contribution_dp']

                mtxF, invC, detF = cache(['F', 'invC', 'detF'],
                                         self, 0, state=state)

                bf = aps.get_base('v', 0, self.integral)
                for out, chunk in self.char_fun(chunk_size, shape):
                    status = fun(out, bf, mtxF, invC, detF, vgv, 1, chunk, 1)
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
        detF, invC = family_data

        p, = self.get_args(['state_p'], **kwargs)

        cache = self.get_cache('state_in_volume_qp', 0)
        p_qp = cache('state', self, 0, state=p, get_vector=self.get_vector)

        if mode == 0:
            out = nm.empty_like(invC)
            fun = self.term_function['stress']
        elif mode == 1:
            shape = list(invC.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)
            fun = self.term_function['tangent_modulus_u']
        else:
            raise ValueError('bad mode! (%d)' % mode)

        fun(out, p_qp, detF, invC)

        return out

class VolumeTLTerm(CouplingVectorScalarTL, InstantaneousBase, Term):
    r"""
    :Description:
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
        virtual : :math:`q`,
        state   : :math:`\ul{u}`
    """
    name = 'dw_tl_volume'
    arg_types = ('virtual', 'state')
    use_caches = {'finite_strain_tl' : [['state',
                                         {'F' : (2, 2),
                                          'invC' : (2, 2),
                                          'detF' : (2, 2)}]]}

    function = staticmethod(terms.dw_tl_volume)

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        virtual, state = self.get_args( **kwargs )
        term_mode = kwargs.get('term_mode')

        apv, vgv = self.get_approximation(state)
        aps, vgs = self.get_approximation(virtual)

        self.set_data_shape(apv, aps)

        cache = self.get_cache('finite_strain_tl', 0)
        ih = self.arg_steps[state.name] # issue 104!
        mtxF, invC, detF = cache(['F', 'invC', 'detF'], self, ih, state=state)

        if term_mode == 'volume':
            n_el, _, _, _ = self.data_shape_s
            shape, mode = (n_el, 1, 1, 1), 2

        elif term_mode == 'rel_volume':
            n_el, _, _, _ = self.data_shape_s
            shape, mode = (n_el, 1, 1, 1), 3

        else:
            shape, mode = self.get_shape_div(diff_var, chunk_size)

        bf = aps.get_base('v', 0, self.integral)

        return (bf, mtxF, invC, detF, vgv, 0), shape, mode

class DiffusionTLTerm(ScalarScalar, Term):
    r"""
    :Description:
    Diffusion term in the total Lagrangian formulation with
    linearized deformation-dependent permeability
    :math:`\ull{K}(\ul{u}) = J \ull{F}^{-1} \ull{k} f(J) \ull{F}^{-T}`,
    where :math:`\ul{u}` relates to the previous time step :math:`(n-1)`
    and
    :math:`f(J) = \max\left(0, \left(1 + \frac{(J - 1)}{N_f}\right)\right)^2`
    expresses the dependence on volume compression/expansion.

    :Definition:
    .. math::
        \int_{\Omega} \ull{K}(\ul{u}^{(n-1)}) : \pdiff{q}{X} \pdiff{p}{X}

    :Arguments:
        material_1 : :math:`\ull{k}`,
        material_2 : :math:`N_f`,
        virtual    : :math:`q`,
        state      : :math:`p`,
        parameter  : :math:`\ul{u}^{(n-1)}`
    """
    name = 'dw_tl_diffusion'
    arg_types = ('material_1', 'material_2', 'virtual', 'state', 'parameter')
    use_caches = {'grad_scalar' : [['state']],
                  'finite_strain_tl' : [['parameter',
                                         {'F' : (2, 2),
                                          'invC' : (2, 2),
                                          'detF' : (2, 2)}]]}

    function = staticmethod(terms.dw_tl_diffusion)

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        perm, ref_porosity, virtual, state, par = self.get_args(**kwargs)
        term_mode = kwargs.get('term_mode')

        apv, vgv = self.get_approximation(par)
        aps, vgs = self.get_approximation(virtual)

        self.set_data_shape(aps)

        cache = self.get_cache('finite_strain_tl', 0)
        # issue 104!
        if self.step == 0:
            ih = 0
        else:
            ih = 1
        mtxF, detF = cache(['F', 'detF'], self, ih, state=par)

        if term_mode == 'diffusion_velocity':
            n_el, n_qp, dim, n_ep = self.data_shape
            shape, mode = (n_el, 1, dim, 1), 2

        else:
            shape, mode = self.get_shape(diff_var, chunk_size)
            if self.step == 0: # Just init the history in step 0.
                raise StopIteration
        
        cache = self.get_cache('grad_scalar', 0)
        gp = cache('grad', self, 0, state=state, get_vector=self.get_vector)

        return (gp, perm, ref_porosity, mtxF, detF, vgv), shape, mode

class SurfaceTractionTLTerm(VectorVector, Term):
    r"""
    :Description:
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
        material : :math:`\ull{\sigma}`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_tl_surface_traction'
    arg_types = ('material', 'virtual', 'state')
    integration = 'surface_extra'
    use_caches = {'finite_strain_surface_tl' : [['state']]}

    function = staticmethod(terms.dw_tl_surface_traction)

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        trac_qp, virtual, state = self.get_args(**kwargs)
        ap, sg = self.get_approximation(virtual)
        sd = ap.surface_data[self.region.name]

        n_fa, n_qp = ap.get_s_data_shape(self.integral,
                                         self.region.name)[:2]
        n_el, dim, n_ep = ap.get_v_data_shape()
        self.data_shape = (n_fa, n_qp, dim, n_ep)
        shape, mode = self.get_shape(diff_var, chunk_size)

        cache = self.get_cache('finite_strain_surface_tl', 0)
        detF, invF = cache(['detF', 'invF'],
                           self, 0, state=state, data_shape=self.data_shape)

        bf = ap.get_base(sd.bkey, 0, self.integral)

        assert_(trac_qp.shape[2] == trac_qp.shape[3] == dim)

        return (trac_qp, detF, invF, bf, sg, sd.fis), shape, mode

    def needs_local_chunk(self):
        return True, False
