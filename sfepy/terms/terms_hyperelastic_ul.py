from sfepy.terms.terms import *
from sfepy.terms.terms_base import VectorVector
from sfepy.terms.terms_hyperelastic_base import HyperElasticBase

class HyperElasticULBase( HyperElasticBase ):
    """Base class for all hyperelastic terms in UL formulation. This is not a
    proper Term!
    """
    use_caches = {'finite_strain_ul' : [['state']]}
    
    def __init__(self, name, sign, **kwargs):
        HyperElasticBase.__init__(self, name, sign, mode='ul', **kwargs)

class NeoHookeanULTerm( VectorVector, HyperElasticULBase ):
    r"""
    :Description:
    Hyperelastic neo-Hookean term. Effective stress
    :math:`\tau_{ij} = \mu J^{-\frac{2}{3}}(b_{ij} - \frac{1}{3}b_{kk}\delta_{ij})`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J
    """
    name = 'dw_ul_he_neohook'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

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

class BulkPenaltyULTerm( VectorVector, HyperElasticULBase ):
    r"""
    :Description:
    Hyperelastic bulk penalty term. Stress
    :math:`\tau_{ij} = K(J-1)\; J \delta_{ij}`.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J
    """
    name = 'dw_ul_bulk_penalty'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

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

class MooneyRivlinULTerm( VectorVector, HyperElasticULBase ):
    r"""
    :Description:
    Hyperelastic Mooney-Rivlin term.

    :Definition:
    .. math::
        \int_{\Omega} \mathcal{L}\tau_{ij}(\ul{u}) e_{ij}(\delta\ul{v})/J
    """
    name = 'dw_ul_he_mooney_rivlin'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

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

