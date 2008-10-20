from terms import *
from utils import choose_scalar_or_in_el

class MassTerm( Term ):
    r""":description: Inertial forces term (constant density).
    :definition: $\int_{\Omega} \rho \ul{v} \cdot \frac{\ul{u} -
    \ul{u}_0}{\dt}$
    :arguments: material.rho : $\rho$, ts.dt : $\dt$, parameter : $\ul{u}_0$"""
    name = 'dw_mass'
    arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name )
        n_el, n_qp, dim, n_ep = self.data_shape

        if diff_var is None:
            return (chunk_size, 1, dim * n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, dim * n_ep, dim * n_ep), 1
        else:
            raise StopIteration
        
    def build_c_fun_args( self, mat, state, ap, vg ):
        # terms.dw_mass_rho_in_el is missing
        mat, self.function = choose_scalar_or_in_el( mat, nm.float64,
                                                     terms.dw_mass,
                                                     NotImplemented )
        ts, state0 = self.get_args( ['ts', 'parameter'], **kwargs )

        dvec = state() - state0()
        rhodt = mat / ts.dt
        bf = ap.get_base( 'v', 0, self.integral_name )
        return rhodt, dvec, 0, bf, vg, ap.econn

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( ['material', 'virtual', 'state'],
                                             **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, ap )
        fargs = self.build_c_fun_args( mat, state, ap, vg )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status

class MassVectorTerm( MassTerm ):
    r""":description: Vector field mass matrix/rezidual.
    :definition: $\int_{\Omega} \rho\ \ul{v} \cdot \ul{u}$
    """
    name = 'dw_mass_vector'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    def build_c_fun_args( self, mat, state, ap, vg ):
        mat, self.function = choose_scalar_or_in_el( mat, nm.float64,
                                                     terms.dw_mass,
                                                     NotImplemented )
        vec = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        return mat, vec, 0, bf, vg, ap.econn

class MassScalarTerm( Term ):
    r""":description: Scalar field mass matrix/rezidual.
    :definition: $\int_{\Omega} q p$
    """
    name = 'dw_mass_scalar'
    arg_types = ('virtual', 'state')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_mass_scalar )
        
    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name )
        n_el, n_qp, dim, n_ep = self.data_shape

        if diff_var is None:
            return (chunk_size, 1, n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, n_ep, n_ep), 1
        else:
            raise StopIteration
        
    def build_c_fun_args( self, state, ap, vg, **kwargs ):
        vec = self.get_vector( state )
        bf = ap.get_base( 'v', 0, self.integral_name )

        if state.is_real():
            return vec, 0, bf, vg, ap.econn
        else:
            ac = nm.ascontiguousarray
            return [(ac( vec.real ), 0, bf, vg, ap.econn),
                    (ac( vec.imag ), 0, bf, vg, ap.econn)]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, ap )
        fargs = self.build_c_fun_args( state, ap, vg, **kwargs )

        if state.is_real():
            for out, chunk in self.char_fun( chunk_size, shape ):
                status = self.function( out, *fargs + (chunk, mode) )
                yield out, chunk, status
        else:
            # For mode == 1, the matrix is the same both for real and imaginary
            # part -> optimization possible.
            for out_real, chunk in self.char_fun( chunk_size, shape ):
                out_imag = nm.zeros_like( out_real )
                status1 = self.function( out_real, *fargs[0] + (chunk, mode) )
                status2 = self.function( out_imag, *fargs[1] + (chunk, mode) )
                yield out_real + 1j * out_imag, chunk, status1 or status2

class MassScalarVariableTerm( MassScalarTerm ):
    r""":description: Scalar field mass matrix/rezidual with coefficient $c$
    defined in nodes.
    :definition: $\int_{\Omega} c q p$
    """
    name = 'dw_mass_scalar_variable'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    use_caches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_mass_scalar_variable )
        
    def build_c_fun_args( self, state, ap, vg, **kwargs ):
        n_el, n_qp = self.data_shape[:2]
        
        mat, = self.get_args( ['material'], **kwargs )
        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                        mat = mat, ap = ap,
                        assumed_shapes = [(n_el, n_qp, 1, 1)],
                        mode_in = 'vertex' )

        vec = state()
        bf = ap.get_base( 'v', 0, self.integral_name )
        return mat_qp, vec, 0, bf, vg, ap.econn

class MassScalarFineCoarseTerm( Term ):
    r""":description: Scalar field mass matrix/rezidual for coarse to fine grid
    interpolation. Field $p_H$ belong to the coarse grid, test field $q_h$ to
    the fine grid.
    :definition: $\int_{\Omega} q_h p_H$
    """
    name = 'dw_mass_scalar_fine_coarse'
    arg_types = ('virtual', 'state', 'iemaps', 'pbase' )
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign,
                       terms.dw_mass_scalar_fine_coarse )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state, iemaps, pbase = self.get_args( **kwargs )
        apr, vgr = virtual.get_current_approximation()
        apc, vgc = virtual.get_current_approximation()
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, n_epr, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape()[3]
            shape = (chunk_size, 1, n_epr, n_epc)
            mode = 1
        else:
            raise StopIteration

        vec = state()

        cbfs = pbase[self.char_fun.ig]
        iemap = iemaps[self.char_fun.ig]
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, apr.bf['v'], cbfs,
                                    vgr, apc.econn, iemap, chunk, mode )
            
            yield out, chunk, status
