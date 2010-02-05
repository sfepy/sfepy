from sfepy.terms.terms import *
from sfepy.terms.terms_base import ScalarScalar

class AcousticAlphaSA1Term( ScalarScalar, Term ):
    r"""
    :Description:
    Sensitivity analysis acoustic term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} \partial_\alpha w_k\, \partial_k \ul{u}\, \partial_\alpha \ul{v},
        \alpha = 1,\dots,N-1
    """
    name = 'd_sa_acoustic_alpha'
    arg_types = ('parameter_1', 'parameter_2', 'parameter_3')
    geometry = [(Volume, 'parameter_1')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sa_acoustic_alpha )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2, par3 = self.get_args( ['parameter_1', 'parameter_2', 'parameter_3'], **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )

        fargs = ( par1(), par2(), par3(), vg, ap.econn, 0 )
        return fargs, (chunk_size, 1, 1, 1), 0
    
class AcousticAlphaSA2Term( ScalarScalar, Term ):
    r"""
    :Description:
    Sensitivity analysis acoustic term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} \dvg w \partial_\alpha \ul{u}\, \partial_\alpha \ul{v},
        \alpha = 1,\dots,N-1
    """
    name = 'd_sa_acoustic_alpha2'
    arg_types = ('parameter_1', 'parameter_2', 'parameter_3')
    geometry = [(Volume, 'parameter_1')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sa_acoustic_alpha )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2, par3 = self.get_args( ['parameter_1', 'parameter_2', 'parameter_3'], **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )

        fargs = ( par1(), par2(), par3(), vg, ap.econn, 1 )

        return fargs, (chunk_size, 1, 1, 1), 0

class AcousticZSATerm( ScalarScalar, Term ):
    r"""
    :Description:
    Sensitivity analysis acoustic term (transversal direction).

    :Definition:
    .. math::
        \int_{\Omega} \partial_z w_k\, \partial_k \ul{v}\,\partial_z \ul{u}
    """
    name = 'd_sa_acoustic_z'
    arg_types = ('parameter_1', 'parameter_2', 'parameter_3')
    geometry = [(Volume, 'parameter_1')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sa_acoustic_z )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2, par3 = self.get_args( ['parameter_1', 'parameter_2', 'parameter_3'], **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )

        fargs = ( par1(), par2(), par3(), vg, ap.econn, 0 )

        return fargs, (chunk_size, 1, 1, 1), 0

class AcousticZSA2Term( ScalarScalar, Term ):
    r"""
    :Description:
    Sensitivity analysis acoustic term (transversal direction).

    :Definition:
    .. math::
        \int_{\Omega} \dvg w \partial_z \ul{v}\,\partial_z \ul{u}
    """
    name = 'd_sa_acoustic_z2'
    arg_types = ('parameter_1', 'parameter_2', 'parameter_3')
    geometry = [(Volume, 'parameter_1')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sa_acoustic_z )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2, par3 = self.get_args( ['parameter_1', 'parameter_2', 'parameter_3'], **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )

        fargs = ( par1(), par2(), par3(), vg, ap.econn, 1 )

        return fargs, (chunk_size, 1, 1, 1), 0

class AcousticTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Acoustic term.

    :Definition:
    .. math::
        \int_{\Omega} (p_1 \partial_\alpha \ul{v}\,\partial_\alpha \ul{u} + p_2 \partial_z \ul{v}\,\partial_z \ul{u} ),
        \alpha = 1,\dots,N-1
    """
    name = 'dw_acoustic'
    arg_types = (('material', 'material', 'virtual', 'state'),
                 ('material', 'material', 'parameter_1', 'parameter_2'))
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter_1'), (Volume, 'parameter_2')])
    modes = ('weak', 'eval')
    functions = {'weak': terms.dw_acoustic,
                 'eval': terms.d_acoustic }
        
    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        return (state(), mat1, mat2, vg, ap.econn), shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, par1, par2 = self.get_args( **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )
        self.set_data_shape( ap )

        return (par1(), par2(), mat1, mat2, vg, ap.econn), (chunk_size, 1, 1, 1), 0

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.functions['weak']
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = self.functions['eval']
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )

class AcousticSurfaceTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Acoustic surface term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} \ul{n}\cdot \partial_{\alpha, 1/h z} \ul{y},
        \alpha = 1,\dots,N-1
    """
    name = 'd_acoustic_surface'
    arg_types = ('material', 'material', 'parameter')
    geometry = [(SurfaceExtra, 'parameter')]
        
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_acoustic_surface )
        self.dof_conn_type = 'surface'

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        mat1, mat2, par = self.get_args( **kwargs )
        ap, sg = par.get_approximation( self.get_current_group(), 'Surface' )

        self.set_data_shape( ap )

        fargs = (par(), mat1, mat2, sg, ap.econn)

        return fargs, (chunk_size, 1, 1, 1), 0

class AcousticIntegrateTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Integration of acoustic term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} m  \partial_\alpha \ul{u},
        \alpha = 1,\dots,N-1
    """
    name = 'dw_acoustic_integrate'
    arg_types = ('material', 'virtual')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_acoustic_integrate )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        return (mat, vg, ap.econn), shape, mode

class AcousticEvalAlphaTerm( Term ):
    r"""
    :Description:
    Evaluation of acoustic term (in-plane directions).

    :Definition:
    .. math::
        \int_{\Omega} \partial_{\alpha} \ul{y},
        \alpha = 1,\dots,N-1
    """
    name = 'd_acoustic_alpha'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_acoustic_alpha )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (1, 1, dim-1, 1)
        else:
            raise StopIteration

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, par(), vg, ap.econn, chunk )
            out1 = nm.sum( out, 0 )
            yield out1, chunk, status
