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
    geometry = [(Volume, 'virtual'), (Volume, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_electric_source )

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, parameter = self.get_args( **kwargs )
        apr, vgr = virtual.get_approximation( self.get_current_group(),
                                              'Volume' )
        apc, vgc = parameter.get_approximation( self.get_current_group(),
                                                'Volume' )
        n_el, n_qp, dim, n_ep = apr.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, n_ep, 1)
            mode = 0
        else:
            raise StopIteration

        bfr = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, parameter(), mat, bfr, vgc,
                                    apc.econn, chunk, mode )
            yield out, chunk, status
