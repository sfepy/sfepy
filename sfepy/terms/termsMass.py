import numpy as nm

from sfepy.base.base import assert_
from sfepy.terms.terms import Term, terms

class MassVectorTerm(Term):
    r"""
    Vector field mass matrix/rezidual.

    :Definition:

    .. math::
        \int_{\Omega} \rho\ \ul{v} \cdot \ul{u}

    :Arguments:
        - material : :math:`\rho`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_mass_vector'
    arg_types = ('material', 'virtual', 'state')

    function = staticmethod(terms.dw_mass)

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if diff_var is None:
            val_qp = self.get(state, 'val')
            fmode = 0

        else:
            val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return mat, val_qp, vg.bf, vg, fmode

class MassScalarTerm(Term):
    r"""
    Scalar field mass matrix/rezidual weighted by a scalar function :math:`c`.

    :Definition:

    .. math::
        \int_{\Omega} q p, \int_{\Omega} c q p

    :Arguments 1:
        - material : :math:`c` (optional)
        - virtual : :math:`q`
        - state   : :math:`p`

    :Arguments 2:
        - material : :math:`c` (optional)
        - parameter_1 : :math:`r`
        - parameter_2 : :math:`p`
    """
    name = 'dw_mass_scalar'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')

    def check_shapes(self, material, virtual, state):
        assert_(virtual.n_components == 1)
        assert_(state.n_components == 1)

        if material is not None:
            n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
            assert_(material.shape[1:] == (n_qp, 1, 1))
            assert_((material.shape[0] == 1) or (material.shape[0] == n_el))

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        if material is None:
            material = nm.ones((1, n_qp, 1, 1), dtype=nm.float64)

        if mode == 'weak':
            if diff_var is None:
                val_qp = self.get(state, 'val')
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return material, val_qp, geo.bf, geo, fmode

        elif mode == 'eval':
            val_qp1 = self.get(virtual, 'val')
            val_qp2 = self.get(state, 'val')

            return material, val_qp1, val_qp2, geo.bf, geo

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, material, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.dw_mass_scalar

        else:
            self.function = terms.d_mass_scalar

class MassScalarSurfaceTerm(MassScalarTerm):
    r"""
    Scalar field mass matrix/rezidual on a surface weighted by a scalar
    function.

    :Definition:

    .. math::
        \int_{\Gamma} q p, \int_{\Gamma} c q p

    :Arguments:
        - material : :math:`c` (optional)
        - virtual : :math:`q`
        - state   : :math:`p`
    """
    name = 'dw_surface_mass_scalar'
    arg_types = ('opt_material', 'virtual', 'state')
    integration = 'surface'

    function = staticmethod(terms.dw_surf_mass_scalar)

    def set_arg_types(self):
        pass

class BCNewtonTerm(MassScalarSurfaceTerm):
    r"""
    Newton boundary condition term.

    :Definition:

    .. math::
        \int_{\Gamma} \alpha q (p - p_{\rm outer})

    :Arguments:
        - material_1 : :math:`\alpha`
        - material_2 : :math:`p_{\rm outer}`
        - virtual    : :math:`q`
        - state      : :math:`p`
    """
    name = 'dw_bc_newton'
    arg_types = ('material_1', 'material_2', 'virtual', 'state')

    def check_shapes(self, alpha, p_outer, virtual, state):
        pass

    def get_fargs(self, alpha, p_outer, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        fargs = MassScalarSurfaceTerm.get_fargs(self, virtual, state,
                                                mode, term_mode, diff_var,
                                                material=alpha, **kwargs)
        fargs = fargs[:1] + (fargs[1] - p_outer,) + fargs[2:]

        return fargs
