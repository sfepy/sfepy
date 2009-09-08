from sfepy.terms.terms import *

class LinearVolumeForceTerm( Term ):
    r""":description: Vector or scalar linear volume forces (weak form) --- a
    right-hand side source term.
    :definition: $\int_{\Omega} \ul{f} \cdot \ul{v}$ or $\int_{\Omega} f q$
    """
    name = 'dw_volume_lvf'
    arg_types = ('material', 'virtual')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_volume_lvf )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )
        vdim = virtual.field.shape[0]
        
        if diff_var is None:
            shape = (chunk_size, 1, vdim * n_ep, 1)
            mode = 0
        else:
            raise StopIteration

        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, bf, mat, vg, chunk )
            yield out, chunk, status
