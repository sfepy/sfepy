from terms import *

class PiezoCouplingTerm( Term ):
    r""":description: Piezoelectric coupling term.
    :definition: $\int_{\Omega} g_{kij}\ e_{ij}(\ul{u}) \nabla_k q$,
    $\int_{\Omega} g_{kij}\ e_{ij}(\ul{v}) \nabla_k p$
    """
    name = 'dw_piezo_coupling'
    arg_types = ('material', 'virtual|state', 'state|virtual')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

    def set_arg_types( self ):
        if self.ats[1] == 'virtual':
            self.mode = 'grad'
#            self.function = terms.dw_piezo_coupling
            insert_as_static_method( self.__class__, 'get_shape',
                                     self.get_shape_grad )
        elif self.ats[2] == 'virtual':
            self.mode = 'div'
#            self.function = terms.dw_piezo_coupling
            insert_as_static_method( self.__class__, 'get_shape',
                                     self.get_shape_div )
        else:
            self.mode = 'eval'
            raise NotImplementedError

    def get_shape( self, diff_var, chunk_size, apr, apc = None ):
        """This is either get_shape_grad() or get_shape_div()."""
        pass
    
    def get_shape_grad( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name ) 
        n_el, n_qp, dim, n_epr = self.data_shape

        if diff_var is None:
            return (chunk_size, 1, dim * n_epr, 1 ), 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            return (chunk_size, 1, dim * n_epr, n_epc ), 1
        else:
            raise StopIteration

    def get_shape_div( self, diff_var, chunk_size, apr, apc = None ):
        self.data_shape = apr.get_v_data_shape( self.integral_name ) 
        n_el, n_qp, dim, n_epr = self.data_shape

        if diff_var is None:
            return (chunk_size, 1, n_epr, 1 ), 2
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            return (chunk_size, 1, n_epr, dim * n_epc ), 3
        else:
            raise StopIteration

    def build_c_fun_args( self, state, apc, vgr, **kwargs ):
        vec = state()
        bf = apc.get_base( 'v', 0, self.integral_name )
        return 1.0, vec, 0, bf, vgr, apc.econn

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )

        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = state.get_approximation( self.get_current_group(),
                                            'Volume' )

        shape, mode = self.get_shape( diff_var, chunk_size, apr, apc )
        debug()
        fargs = self.build_c_fun_args( state, apc, vgr, **kwargs )
        
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status
