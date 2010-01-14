from sfepy.terms.terms import *
from sfepy.terms.terms_hyperelastic_base \
     import VectorVector, CouplingVectorScalarTL, HyperElasticBase
            
class HyperElasticTLBase( HyperElasticBase ):
    """Base class for all hyperelastic terms in TL formulation. This is not a
    proper Term!
    """
    use_caches = {'finite_strain_tl' : [['state']]}

    def __init__( self, region, name = None, sign = 1 ):
        HyperElasticBase.__init__( self, region, name, sign, mode_ul = 0 )

class NeoHookeanTLTerm( VectorVector, HyperElasticTLBase ):
    r""":description: Hyperelastic neo-Hookean term. Effective stress $S_{ij} =
    \mu J^{-\frac{2}{3}}(\delta_{ij} - \frac{1}{3}C_{kk}C_{ij}^{-1})$.
    :definition:
    $\int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})$
    """
    name = 'dw_tl_he_neohook'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

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
    r""":description: Hyperelastic Mooney-Rivlin term. Effective stress $S_{ij}
    = \kappa J^{-\frac{4}{3}} (C_{kk} \delta_{ij} - C_{ij} - \frac{2}{3
    } I_2 C_{ij}^{-1})$.
    :definition:
    $\int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})$
    """
    name = 'dw_tl_he_mooney_rivlin'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

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
    r""":description: Hyperelastic bulk penalty term. Stress $S_{ij}
    = K(J-1)\; J C_{ij}^{-1}$.
    :definition:
    $\int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})$
    """

    name = 'dw_tl_bulk_penalty'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

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
    r""":description: Hyperelastic bulk pressure term. Stress $S_{ij}
    = -p J C_{ij}^{-1}$.
    :definition:
    $\int_{\Omega} S_{ij}(p) \delta E_{ij}(\ul{u};\ul{v})$
    """

    name = 'dw_tl_bulk_pressure'
    arg_types = ('virtual', 'state', 'state_p')
    geometry = [(Volume, 'virtual')]
    use_caches = {'finite_strain_tl' : [['state']],
                  'state_in_volume_qp' : [['state_p']]}

    family_data_names = ['detF', 'invC']
    ## term_function = {'stress' : terms.dq_tl_stress_pressure,
    ##                  'tangent_modulus' : terms.dq_tl_tan_mod_pressure }

    def __init__(self, region, name=None, sign=1):
        Term.__init__(self, region, name, sign, function=None)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        call_mode, = self.get_kwargs(['call_mode'], **kwargs)
        virtual, state_u, state_p = self.get_args(**kwargs)
        apv, vgv = virtual.get_approximation(self.get_current_group(), 'Volume')
        aps, vgs = state_p.get_approximation(self.get_current_group(), 'Volume')

        self.set_data_shape(apv, aps)
        shape, mode = self.get_shape(diff_var, chunk_size)

        cache = self.get_cache('finite_strain_tl', 0)
        family_data = cache(self.family_data_names,
                            self.get_current_group(), 0, state=state_u)

        if call_mode is None:

            crt_data = self.compute_crt_data(family_data, mode, **kwargs)
            fun = self.function['element_contribution']

            invC, detF = cache(['invC', 'detF'],
                               self.get_current_group(), 0, state=state_u)
            for out, chunk in self.char_fun(chunk_size, shape):
                status = fun(out, crt_data, invC, detF, vgv, chunk, mode)
                yield out, chunk, status

        elif call_mode == 'd_eval':
            raise NotImplementedError

        elif call_mode in ['de_strain', 'de_stress']:

            if call_mode == 'de_strain':
                out_qp = cache('E', self.get_current_group(), 0, state=state_u)
            elif call_mode == 'de_stress':
                out_qp = self.compute_crt_data(family_data, 0, **kwargs)
                
            shape = (chunk_size, 1) + out_qp.shape[2:]
            for out, chunk in self.char_fun(chunk_size, shape):
                status = vg.integrate_chunk(out, out_qp[chunk], chunk)
                out1 = out / vg.variable(2)[chunk]

            yield out1, chunk, status

    def compute_crt_data(self, family_data, mode, **kwargs):
        detF, invC = family_data

        p, = self.get_kwargs(['state_p'], **kwargs)

        cache = self.get_cache('state_in_volume_qp', 0)
        p_qp = cache('state', self.get_current_group(), 0,
                     state=p, get_vector=self.get_vector)
        
        if mode == 0:
            out = nm.empty_like(invC)
            fun = self.term_function['stress']
        elif mode == 1:
            shape = list(invC.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)
            fun = self.term_function['tangent_modulus_u']
        else:
            fun = self.term_function['tangent_modulus_p']

        fun(out, p_qp, detF, invC)

        return out
