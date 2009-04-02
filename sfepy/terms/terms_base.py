from sfepy.terms.terms import *

class InstantaneousBase( Struct ):
    
    def _call( self, diff_var = None, chunk_size = None, **kwargs ):
        call_mode, = self.get_kwargs( ['call_mode'], **kwargs )

        # Imaginary mode marks problem in complex numbers.
        fargs, shape, mode = self.get_fargs( diff_var, chunk_size, **kwargs )

        if call_mode is None:
            if nm.isreal( mode ):
                for out, chunk in self.char_fun( chunk_size, shape ):
                    if self.dof_conn_type == 'surface':
                        chunk = self.char_fun.get_local_chunk()
                    status = self.function( out, *fargs + (chunk, mode) )
                    yield out, chunk, status
            else:
                # Assuming linear forms. Then for mode == 1, the matrix is the
                # same both for real and imaginary part.
                rmode = int( mode.real )
                for out_real, chunk in self.char_fun( chunk_size, shape ):
                    if self.dof_conn_type == 'surface':
                        chunk = self.char_fun.get_local_chunk()
                    status1 = self.function( out_real,
                                             *fargs[0] + (chunk, rmode ) )
                    if rmode == 0:
                        out_imag = nm.zeros_like( out_real )
                        status2 = self.function( out_imag,
                                                 *fargs[1] + (chunk, rmode) )
                        yield out_real + 1j * out_imag, chunk, status1 or status2

                    else:
                        yield out_real, chunk, status1

        elif call_mode == 'd_eval':
            if nm.isreal( mode ):
                for out, chunk in self.char_fun( chunk_size, shape ):
                    status = self.function( out, *fargs + (chunk,) )
                    out1 = nm.sum( out )
                    yield out1, chunk, status
            else:
                raise NotImplementedError

        else:
            msg = 'unknown call_mode for %s' % self.name
            raise ValueError( msg )

class TimeHistoryBase( Struct ):
    
    def _call( self, diff_var = None, chunk_size = None, **kwargs ):
        call_mode, = self.get_kwargs( ['call_mode'], **kwargs )

        fargs, shape, mode = self.get_fargs( diff_var, chunk_size, **kwargs )
        
        if call_mode is None:
            if mode == 1:
                for out, chunk in self.char_fun( chunk_size, shape ):
                    status = self.function( out, *fargs + (chunk, mode) )
                    yield out, chunk, status
            else:
                iter_kernel = fargs
                for out, chunk in self.char_fun( chunk_size, shape,
                                                 zero = True ):
                    out1 = nm.empty_like( out )
                    for ii, fargs in iter_kernel():
                        status = self.function( out1, *fargs + (chunk, mode) )
                        out += out1
                    yield out, chunk, status

        else:
            msg = 'unknown call_mode for %s' % self.name
            raise ValueError( msg )

class VectorVector( InstantaneousBase ):

    def set_data_shape( self, apr, apc = None ):
    
        if self.dof_conn_type == 'surface':
            self.data_shape = apr.get_s_data_shape( self.integral_name,
                                                    self.region.name )
        else:
            self.data_shape = apr.get_v_data_shape( self.integral_name )

        assert_( apr.dim == (self.data_shape[2], 1) )

    def get_shape( self, diff_var, chunk_size ):
        n_el, n_qp, dim, n_ep = self.data_shape
        
        if diff_var is None:
            return (chunk_size, 1, dim * n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, dim * n_ep, dim * n_ep), 1
        else:
            raise StopIteration

class VectorVectorTH( TimeHistoryBase, VectorVector ):
    pass

class ScalarScalar( InstantaneousBase ):

    def set_data_shape( self, apr, apc = None ):

        if self.dof_conn_type == 'surface':
            self.data_shape = apr.get_s_data_shape( self.integral_name,
                                                    self.region.name )
        else:
            self.data_shape = apr.get_v_data_shape( self.integral_name )

        assert_( apr.dim == (1, 1) )

    def get_shape( self, diff_var, chunk_size ):
        n_el, n_qp, dim, n_ep = self.data_shape
        
        if diff_var is None:
            return (chunk_size, 1, n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, n_ep, n_ep), 1
        else:
            raise StopIteration

class ScalarScalarTH( TimeHistoryBase, ScalarScalar ):
    pass

class VectorOrScalar( InstantaneousBase ):

    def set_data_shape( self, apr, apc = None ):

        if self.dof_conn_type == 'surface':
            self.data_shape = apr.get_s_data_shape( self.integral_name,
                                                    self.region.name )
        else:
            self.data_shape = apr.get_v_data_shape( self.integral_name )

        self.vdim = apr.dim[0]

    def get_shape( self, diff_var, chunk_size ):
        n_el, n_qp, dim, n_ep = self.data_shape
        vdim = self.vdim
        
        if diff_var is None:
            return (chunk_size, 1, vdim * n_ep, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, vdim * n_ep, vdim * n_ep), 1
        else:
            raise StopIteration
    
class VectorOrScalarTH( TimeHistoryBase, VectorOrScalar ):
    pass

class CouplingVectorScalar( InstantaneousBase ):

    def set_data_shape( self, apr, apc ):
        """Set element data shape and checks dimensions of approximations."""
        if self.dof_conn_type == 'surface':
            self.data_shape_r = apr.get_s_data_shape( self.integral_name,
                                                      self.region.name )
            self.data_shape_c = apc.get_s_data_shape( self.integral_name,
                                                      self.region.name )
        else:
            self.data_shape_r = apr.get_v_data_shape( self.integral_name )
            self.data_shape_c = apc.get_v_data_shape( self.integral_name )

        if self.mode == 'grad':
            dim_c = (1, 1)
            dim_r = (self.data_shape_r[2], 1)
        else: # 'div', 'eval'
            dim_c = (self.data_shape_r[2], 1)
            dim_r = (1, 1)
        assert_( apc.dim == dim_c )
        assert_( apr.dim == dim_r )

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


class CouplingVectorScalarTH( TimeHistoryBase, CouplingVectorScalar ):
    pass
