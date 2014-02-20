from sfepy.terms.terms import Term, terms

class LinearVolumeForceTerm(Term):
    r"""
    Vector or scalar linear volume forces (weak form) --- a right-hand side
    source term.

    :Definition:

    .. math::
        \int_{\Omega} \ul{f} \cdot \ul{v} \mbox{ or } \int_{\Omega} f q

    :Arguments:
        - material : :math:`\ul{f}` or :math:`f`
        - virtual  : :math:`\ul{v}` or :math:`q`
    """
    name = 'dw_volume_lvf'
    arg_types = ('material', 'virtual')
    arg_shapes = [{'material' : 'D, 1', 'virtual' : ('D', None)},
                  {'material' : '1, 1', 'virtual' : (1, None)}]

    function = staticmethod(terms.dw_volume_lvf)

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        return mat, vg
