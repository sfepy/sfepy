from sfepy.terms.terms import *

class ElectricSourceTerm( Term ):
    r"""
    :Description:
    Electric source term.

    :Definition:
    .. math::
        \int_{\Omega} c s (\nabla \phi)^2

    :Arguments:
    material : :math:`c` (electric conductivity),
    virtual : :math:`s` (test function),
    parameter : :math:`\phi` (given electric potential)
    """
    name = 'dw_electric_source'
    arg_types = ('material', 'virtual', 'parameter')

    function = staticmethod(terms.dw_electric_source)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, parameter = self.get_args( **kwargs )
        apr, vgr = self.get_approximation(virtual)
        apc, vgc = self.get_approximation(parameter)
        n_el, n_qp, dim, n_ep = apr.get_v_data_shape(self.integral)

        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
            mode = 0
        else:
            raise StopIteration

        bfr = apr.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, parameter(), mat, bfr, vgc,
                                    apc.econn, chunk, mode )
            yield out, chunk, status
