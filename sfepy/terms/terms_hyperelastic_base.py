from sfepy.terms.terms import *

class CouplingVectorScalarTL(Struct):
    
    def set_data_shape(self, apv, aps):
        """Set element data shape and checks dimensions of approximations."""
        if self.dof_conn_type == 'surface':
            self.data_shape_v = apv.get_s_data_shape(self.integral_name,
                                                     self.region.name)
            self.data_shape_s = aps.get_s_data_shape(self.integral_name,
                                                     self.region.name)
        else:
            self.data_shape_v = apv.get_v_data_shape(self.integral_name)
            self.data_shape_s = aps.get_v_data_shape(self.integral_name)

        assert_(aps.dim == (1,))
        assert_(apv.dim == (self.data_shape_v[2],))

    def get_shape_div(self, diff_var, chunk_size):
        n_el, n_qp, dim, n_eps = self.data_shape_s

        if diff_var is None:
            return (chunk_size, 1, n_eps, 1), 0
        elif diff_var == self.get_arg_name( 'state_p' ):
            return (chunk_size, 1, n_eps, n_eps), 2
        elif diff_var == self.get_arg_name( 'state' ):
            n_epv = self.data_shape_v[3]
            return (chunk_size, 1, n_eps, dim * n_epv), 1
        else:
            raise StopIteration

    def get_shape_grad(self, diff_var, chunk_size):
        n_el, n_qp, dim, n_epv = self.data_shape_v

        if diff_var is None:
            return (chunk_size, 1, dim * n_epv, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, dim * n_epv, dim * n_epv), 1
        elif diff_var == self.get_arg_name( 'state_p' ):
            n_eps = self.data_shape_s[3]
            return (chunk_size, 1, dim * n_epv, n_eps), 2
        else:
            raise StopIteration

class HyperElasticBase( Term ):
    """Base class for all hyperelastic terms in TL/UL formulation.

    **Note** This is not a proper Term!

    `HyperElasticBase.__call__()` computes element contributions given either
    stress (-> rezidual) or tangent modulus (-> tangent sitffnes matrix),
    i.e. constitutive relation type (CRT) related data. The CRT data are
    computed in subclasses implementing particular CRT (e.g. neo-Hookean
    material), in self.compute_crt_data().
    Mode: 0 - total formulation, 1 - updated formulation
    """
    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.mode_ul = {'tl' : 0, 'ul' : 1}[self.mode]
        self.function = {
            'finite_strain' : { 0: 'finite_strain_tl',
                                1: 'finite_strain_ul'},
            'element_contribution' : terms.dw_he_rtm,
        }

        self.crt_data = Struct( stress = None,
                                tan_mod = nm.array( [0], ndmin = 4 ) )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        call_mode, = self.get_kwargs( ['call_mode'], **kwargs )
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        cache = self.get_cache( self.function['finite_strain'][self.mode_ul], 0 )
        family_data = cache( self.family_data_names,
                             self.get_current_group(), 0, state = state )
##         print family_data

        if call_mode is None:

            out = self.compute_crt_data( family_data, mode, **kwargs )
            if mode == 0:
                self.crt_data.stress = out
            else:
                self.crt_data.tan_mod = out

            fun = self.function['element_contribution']

            mtxF, detF = cache( ['F', 'detF'],
                                self.get_current_group(), 0, state = state )
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = fun( out, self.crt_data.stress, self.crt_data.tan_mod,
                              mtxF, detF, vg, chunk, mode, self.mode_ul )
                yield out, chunk, status

        elif call_mode == 'd_eval':
            raise NotImplementedError

        elif call_mode in ['de_strain', 'de_stress']:

            if call_mode == 'de_strain':
                out_qp = cache( 'E', self.get_current_group(), 0, state = state )
            elif call_mode == 'de_stress':
                out_qp = self.compute_crt_data( family_data, 0, **kwargs )
                
            shape = (chunk_size, 1) + out_qp.shape[2:]
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = vg.integrate_chunk( out, out_qp[chunk], chunk )
                out1 = out / vg.variable( 2 )[chunk]

            yield out1, chunk, status

class DeformationGradientTerm(Term):
    r"""
    :Description:
    Deformation gradient :math:`F` in quadrature points for
    `call_mode='dq_def_grad'` (default) or the jacobian :math:`J` if
    `call_mode='dq_jacobian'`.

    :Definition:
    .. math::
        \ull{F} = \pdiff{\ul{x}}{\ul{X}}|_{qp}
        = \ull{I} + \pdiff{\ul{u}}{\ul{X}}|_{qp} \;, \\
        \ul{x} = \ul{X} + \ul{u} \;, J = \det{(\ull{F})}

    :Arguments:
        state : :math:`\ul{u}`
    """
    name = 'dq_def_grad'
    arg_types = ('state',)
    geometry = [(Volume, 'state')]

    function = staticmethod(terms.dq_def_grad)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        state, = self.get_args(**kwargs)
        call_mode = kwargs.get('call_mode', 'dq_def_grad')

        ap, vg = state.get_approximation(self.get_current_group(), 'Volume')
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral_name)

        if diff_var is None:
            if call_mode == 'dq_def_grad':
                shape = (chunk_size, n_qp, dim, dim)
                mode = 0

            elif call_mode == 'dq_jacobian':
                shape = (chunk_size, n_qp, 1, 1)
                mode = 1

        else:
            raise StopIteration

        vec = state()
        for out, chunk in self.char_fun(chunk_size, shape):
            status = self.function(out, vec, vg, ap.econn, chunk, mode)
            yield out, chunk, status
