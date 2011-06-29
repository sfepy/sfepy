import numpy as nm

from sfepy.terms.terms import Term

class LinearPointSpringTerm(Term):
    r"""
    :Description:
    Linear springs constraining movement of FE nodes in a region; to use as a
    relaxed Dirichlet boundary conditions.

    :Definition:
    .. math::
        \ul{f}^i = -k \ul{u}^i \quad \forall \mbox{ FE node } i \mbox{ in
        a region }

    :Arguments:
        material : :math:`k`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_point_lspring'
    arg_types = ('material', 'virtual', 'state')
    integration = 'point'

    def get_integral_info(self):
        """
        Get information on the term integral.

        Returns
        -------
        kind : 'v' or 's'
            The integral kind.
        """
        return 'v'

    @staticmethod
    def function(out, stiffness, vec, diff_var):
        if diff_var is None:
            out[:, 0, :, 0] = - stiffness * vec

        else:
            dim = vec.shape[1]
            eye = nm.eye(dim, dim, dtype=nm.float64)
            eye.shape = (1, 1) + eye.shape
            out[...] = - stiffness * eye

        return 0

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        n_dof, _, _, _, n_c = self.get_data_shape(state)

        vec = state.get_state_in_region(self.region)

        stiffness = mat['stiffness']

        return stiffness, vec, diff_var
