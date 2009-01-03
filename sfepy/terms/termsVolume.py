from sfepy.terms.terms import *
from sfepy.terms.utils import fix_scalar_constant

class LinearVolumeForceTerm( Term ):
    r""":description: Vector or scalar linear volume forces (weak form) --- a
    right-hand side source term.
    :definition: $\int_{\Omega} \ul{f} \cdot \ul{v}$ or $\int_{\Omega} f q$
    """
    name = 'dw_volume_lvf'
    arg_types = ('material', 'virtual')
    geometry = [(Volume, 'virtual')]
    use_caches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_volume_lvf )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        force, virtual = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )
        vdim = virtual.field.dim[0]
        
        if diff_var is None:
            shape = (chunk_size, 1, vdim * n_ep, 1)
            mode = 0
        else:
            raise StopIteration

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat = fix_scalar_constant( force, nm.float64 )
        if mat is None:
            mat = nm.asarray( force, dtype = nm.float64 )
            if mat.ndim == 1:
                mat = nm.ascontiguousarray( mat[...,nm.newaxis] )
        else:
            try:
                mat = mat[...,nm.newaxis,nm.newaxis]
            except TypeError: # Old numpy...
                mat = nm.array( mat, ndmin = 2 )

        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = mat, ap = ap,
                       assumed_shapes = [(n_el, n_qp, vdim, 1)],
                       mode_in = None )
#        print mat_qp
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, bf, mat_qp, vg, chunk )
            yield out, chunk, status
