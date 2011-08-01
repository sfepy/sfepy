from sfepy.terms.terms import Term, terms

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

    def get_fargs(self, mat, virtual, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        grad = self.get(parameter, 'grad')

        return grad, mat, vg.bf, vg
