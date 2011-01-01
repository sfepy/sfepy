from sfepy.terms.terms import Term, terms

class LinearVolumeForceTerm( Term ):
    r"""
    :Description:
    Vector or scalar linear volume forces (weak form) --- a right-hand side
    source term.

    :Definition:
    .. math::
        \int_{\Omega} \ul{f} \cdot \ul{v} \mbox{ or } \int_{\Omega} f q

    :Arguments:
        material : :math:`\ul{f}` or :math:`f`,
        virtual  : :math:`\ul{v}` or :math:`q`
    """
    name = 'dw_volume_lvf'
    arg_types = ('material', 'virtual')

    function = staticmethod(terms.dw_volume_lvf)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)
        vdim = virtual.field.shape[0]
        
        if diff_var is None:
            shape = (chunk_size, 1, vdim * n_ep, 1)
            mode = 0
        else:
            raise StopIteration

        bf = ap.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, bf, mat, vg, chunk )
            yield out, chunk, status
