from terms import *

class CouplingVectorScalar( Struct ):

    def set_data_shape( self, apr, apc ):
        """TODO: shape checking"""
        self.data_shape_r = apr.get_v_data_shape( self.integral_name ) 
        self.data_shape_c = apc.get_v_data_shape( self.integral_name ) 

##         if self.mode == 'grad':
##             dim_c = (1, 1)
##             dim_r = (self.data_shape_r[2], 1)
##         else: # 'div', 'eval'
##             dim_c = (self.data_shape_r[2], 1)
##             dim_r = (1, 1)
##         assert_( apc.field.dim == dim_c )
##         assert_( apr.field.dim == dim_r )

    def get_shape_div( self, diff_var, chunk_size, mode0 = 0 ):
        n_el, n_qp, dim, n_epr = self.data_shape_r

        if diff_var is None:
            return (chunk_size, 1, n_epr, 1 ), 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = self.data_shape_c[3]
            return (chunk_size, 1, n_epr, dim * n_epc ), 1
        else:
            raise StopIteration

    def get_shape_grad( self, diff_var, chunk_size ):
        n_el, n_qp, dim, n_epr = self.data_shape_r

        if diff_var is None:
            return (chunk_size, 1, dim * n_epr, 1 ), 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = self.data_shape_c[3]
            return (chunk_size, 1, dim * n_epr, n_epc ), 1
        else:
            raise StopIteration

    def _call( self, diff_var = None, chunk_size = None, **kwargs ):

        fargs, shape, mode = self.get_fargs( diff_var, chunk_size, **kwargs )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, *fargs + (chunk, mode) )
            yield out, chunk, status
        
