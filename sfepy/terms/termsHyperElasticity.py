from terms import *
from utils import fix_mat_qp_shape

class HyperElasticBase( Term ):
    """Base class for all hyperelastic terms in TL formulation. This is not a
    proper Term!

    HyperElasticBase.__call__() computes element contributions given either
    stress (-> rezidual) or tangent modulus (-> tangent sitffnes matrix),
    i.e. constitutive relation type (CRT) related data. The CRT data are
    computed in subclasses implementing particular CRT (e.g. neo-Hookean
    material), in self.compute_crt_data().
    """
    def __init__( self, region, name = None, sign = 1 ):
        Term.__init__( self, region, name, sign )

        self.function = {
            'element_contribution' : terms.dw_tl_he_rtm,
        }
        self.crt_data = Struct( stress = None,
                                tan_mod = nm.array( [0], ndmin = 4 ) )

    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name )
        n_el, n_qp, dim, n_ep = self.data_shape
        
        if diff_var is None:
            return (chunk_size, 1, dim * n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, dim * n_ep, dim * n_ep), 1
        else:
            raise StopIteration

    def mat_to_qp( self, mat, ap ):
        n_el, n_qp, dim, n_ep = self.data_shape

        mat = nm.asarray( mat )
        if mat.ndim == 0:
            mat = mat[nm.newaxis,nm.newaxis]

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                        mat = mat, ap = ap,
                        assumed_shapes =  [(1, n_qp, 1, 1),
                                           (n_el, n_qp, 1, 1)],
                        mode_in = None )

        mat_qp = fix_mat_qp_shape( mat_qp, n_el )

        return mat_qp

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        call_mode, = self.get_kwargs( ['call_mode'], **kwargs )
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, ap )

        cache = self.get_cache( 'finite_strain_tl', 0 )
        family_data = cache( self.family_data_names,
                             self.get_current_group(), 0, state = state )
##         print family_data

        if call_mode is None:

            out = self.compute_crt_data( family_data, ap, vg, mode, **kwargs )
            if mode == 0:
                self.crt_data.stress = out
            else:
                self.crt_data.tan_mod = out

            fun = self.function['element_contribution']

            mtxF = cache( 'F', self.get_current_group(), 0, state = state )
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = fun( out, self.crt_data.stress, self.crt_data.tan_mod,
                              mtxF, vg, chunk, mode )
                yield out, chunk, status

        elif call_mode == 'd_eval':
            raise NotImplementedError

        elif call_mode in ['de_strain', 'de_stress']:

            if call_mode == 'de_strain':
                out_qp = cache( 'E', self.get_current_group(), 0, state = state )
            elif call_mode == 'de_stress':
                out_qp = self.compute_crt_data( family_data, ap, vg, 0,
                                                **kwargs )
                
            shape = (chunk_size, 1) + out_qp.shape[2:]
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = vg.integrate_chunk( out, out_qp[chunk], chunk )
                out1 = out / vg.variable( 2 )[chunk]

            yield out1, chunk, status

class NeoHookeanTerm( HyperElasticBase ):
    r""":description: Hyperelastic neo-Hookean term. Effective stress $S_{ij} =
    \mu J^{-\frac{2}{3}}(\delta_{ij} - \frac{1}{3}C_{kk}C_{ij}^{-1})$.
    :definition:
    $\int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})$
    """
    name = 'dw_tl_he_neohook'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'finite_strain_tl' : [['state']],
                  'mat_in_qp' : [['material']]}

    family_data_names = ['detF', 'trC', 'invC']

    def __init__( self, region, name = name, sign = 1 ):
        HyperElasticBase.__init__( self, region, name, sign )

        self.function.update( {
            'stress' : terms.dq_tl_he_stress_neohook,
            'tangent_modulus' : terms.dq_tl_he_tan_mod_neohook,
        } )

    def compute_crt_data( self, family_data, ap, vg, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        mat_qp = self.mat_to_qp( mat, ap )

        detF, trC, invC = family_data

        if mode == 0:
            out = nm.empty_like( invC )
            fun = self.function['stress']
        else:
            shape = list( invC.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.function['tangent_modulus']

        fun( out, mat_qp, detF, trC, invC )

        return out

class MooneyRivlinTerm( HyperElasticBase ):
    r""":description: Hyperelastic Mooney-Rivlin term. Effective stress $S_{ij}
    = \kappa J^{-\frac{4}{3}} (C_{kk} \delta_{ij} - C_{ij} - \frac{2}{3
    } I_2 C_{ij}^{-1})$.
    :definition:
    $\int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})$
    """
    name = 'dw_tl_he_mooney_rivlin'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'finite_strain_tl' : [['state']],
                  'mat_in_qp' : [['material']]}

    family_data_names = ['detF', 'trC', 'invC', 'C', 'in2C']

    def __init__( self, region, name = name, sign = 1 ):
        HyperElasticBase.__init__( self, region, name, sign )

        self.function.update( {
            'stress' : terms.dq_tl_he_stress_mooney_rivlin,
            'tangent_modulus' : terms.dq_tl_he_tan_mod_mooney_rivlin,
        } )

    def compute_crt_data( self, family_data, ap, vg, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        mat_qp = self.mat_to_qp( mat, ap )

        detF, trC, invC, vecC, in2C = family_data

        if mode == 0:
            out = nm.empty_like( invC )
            fun = self.function['stress']
        else:
            shape = list( invC.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.function['tangent_modulus']

        fun( out, mat_qp, detF, trC, invC, vecC, in2C )

        return out

class BulkPenaltyTerm( HyperElasticBase ):
    r""":description: Hyperelastic bulk penalty term. Stress $S_{ij}
    = K(J-1)\; J C_{ij}^{-1}$.
    :definition:
    $\int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})$
    """
    name = 'dw_tl_bulk_penalty'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'finite_strain_tl' : [['state']],
                  'mat_in_qp' : [['material']]}

    family_data_names = ['detF', 'invC']

    def __init__( self, region, name = name, sign = 1 ):
        HyperElasticBase.__init__( self, region, name, sign )

        self.function.update( {
            'stress' : terms.dq_tl_he_stress_bulk,
            'tangent_modulus' : terms.dq_tl_he_tan_mod_bulk,
        } )

    def compute_crt_data( self, family_data, ap, vg, mode, **kwargs ):
        mat = self.get_args( ['material'], **kwargs )[0]

        mat_qp = self.mat_to_qp( mat, ap )

        detF, invC = family_data
        
        if mode == 0:
            out = nm.empty_like( invC )
            fun = self.function['stress']
        else:
            shape = list( invC.shape )
            shape[-1] = shape[-2]
            out = nm.empty( shape, dtype = nm.float64 )
            fun = self.function['tangent_modulus']

        fun( out, mat_qp, detF, invC )

        return out
