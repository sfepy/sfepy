import numpy as nm

from sfepy.base.base import assert_
from sfepy.terms.terms import Term

class PointTermBase(Term):
    """
    Common methods of point terms.
    """

    def get_integral_info(self):
        """
        Get information on the term integral.

        Returns
        -------
        kind : 'v' or 's'
            The integral kind.
        """
        return 'v'

class LinearPointSpringTerm(PointTermBase):
    r"""
    Linear springs constraining movement of FE nodes in a region; to use as a
    relaxed Dirichlet boundary conditions.

    :Definition:

    .. math::
        \ul{f}^i = -k \ul{u}^i \quad \forall \mbox{ FE node } i \mbox{ in
        a region }

    :Arguments:
        - material : :math:`k`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_point_lspring'
    arg_types = ('material', 'virtual', 'state')
    integration = 'point'

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
        vec = state.get_state_in_region(self.region)

        stiffness = mat['stiffness']

        return stiffness, vec, diff_var

class ConcentratedPointLoadTerm(PointTermBase):
    r"""
    Concentrated point load term.

    The load value must be given in form of a special material
    parameter (name prefixed with '.'), e.g. (in 2D)::

        'load' : ({'.val' : [0.0, 1.0]},)

    This term should be used with special care, as it bypasses the usual
    evaluation in quadrature points. It should only be used with nodal
    FE basis. The number of rows of the load must be equal to the number
    of nodes in the region and the number of columns equal to the field
    dimension.

    :Definition:

    .. math::
        \ul{f}^i = \ul{\bar f}^i \quad \forall \mbox{ FE node } i \mbox{ in
        a region }

    :Arguments:
        - material : :math:`\ul{\bar f}^i`
        - virtual  : :math:`\ul{v}`,
    """
    name = 'dw_point_load'
    arg_types = ('material', 'virtual')
    integration = 'point'

    @staticmethod
    def function(out, mat):
        out[:, 0, :, 0] = mat

        return 0

    def check_shapes(self, mat, virtual):
        mat = nm.asarray(mat)
        assert_(mat.shape[-1] == virtual.dim)

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return nm.asarray(mat),
