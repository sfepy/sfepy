import numpy as nm

from sfepy.base.base import assert_
from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_th import THTerm, ETHTerm

class DotProductVolumeTerm(Term):
    r"""
    Volume :math:`L^2(\Omega)` weighted dot product for both scalar and vector
    fields. Can be evaluated. Can use derivatives.

    :Definition:

    .. math::
        \int_\Omega q p \mbox{ , } \int_\Omega \ul{v} \cdot \ul{u}
        \mbox{ , }
        \int_\Omega p r \mbox{ , } \int_\Omega \ul{u} \cdot \ul{w} \\
        \int_\Omega c q p \mbox{ , } \int_\Omega c \ul{v} \cdot \ul{u}
        \mbox{ , }
        \int_\Omega c p r \mbox{ , } \int_\Omega c \ul{u} \cdot \ul{w}

    :Arguments 1:
        - material : :math:`c` (optional)
        - virtual  : :math:`q` or :math:`\ul{v}`
        - state    : :math:`p` or :math:`\ul{u}`

    :Arguments 2:
        - material    : :math:`c` (optional)
        - parameter_1 : :math:`p` or :math:`\ul{u}`
        - parameter_2 : :math:`r` or :math:`\ul{w}`
    """
    name = 'dw_volume_dot'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')

    @staticmethod
    def dw_dot(out, mat, val_qp, vgeo, sgeo, fun, fmode):
        status = fun(out, mat, val_qp, vgeo.bf, sgeo.bf, vgeo, fmode)
        return status

    @staticmethod
    def d_dot(out, mat, val1_qp, val2_qp, geo):
        if val1_qp.shape[2] > 1:
            vec = dot_sequences(val1_qp, val2_qp, mode='ATB')

        else:
            vec = val1_qp * val2_qp

        if mat is not None:
            status = geo.integrate(out, mat * vec)
        else:
            status = geo.integrate(out, vec)

        return status

    def check_shapes(self, mat, virtual, state):
        assert_(virtual.n_components == state.n_components)

        if mat is not None:
            n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
            assert_(mat.shape[1:] == (n_qp, 1, 1))
            assert_((mat.shape[0] == 1) or (mat.shape[0] == n_el))

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgeo, _ = self.get_mapping(virtual)

        if mode == 'weak':
            if mat is None:
                n_cell, n_qp, dim, n_n, n_c = self.get_data_shape(state)
                mat = nm.ones((n_cell, n_qp, 1, 1), dtype=nm.float64)

            sgeo, _ = self.get_mapping(state)

            if diff_var is None:
                val_qp = self.get(state, 'val')
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            if state.n_components > 1:
                if self.integration == 'volume':
                    fun = terms.dw_volume_dot_vector

                else:
                    fun = terms.dw_surface_dot_vector

            else:
                if self.integration == 'volume':
                    fun = terms.dw_volume_dot_scalar

                else:
                    fun = terms.dw_surface_dot_scalar

            return mat, val_qp, vgeo, sgeo, fun, fmode

        elif mode == 'eval':
            val1_qp = self.get(virtual, 'val')
            val2_qp = self.get(state, 'val')

            return mat, val1_qp, val2_qp, vgeo

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_cell, n_qp, dim, n_n, n_c = self.get_data_shape(state)

        return (n_cell, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = self.dw_dot

        else:
            self.function = self.d_dot

class DotProductSurfaceTerm(DotProductVolumeTerm):
    r"""
    Surface :math:`L^2(\Gamma)` dot product for both scalar and vector
    fields.

    :Definition:

    .. math::
        \int_\Gamma q p \mbox{ , } \int_\Gamma \ul{v} \cdot \ul{u}
        \mbox{ , }
        \int_\Gamma p r \mbox{ , } \int_\Gamma \ul{u} \cdot \ul{w} \\
        \int_\Gamma c q p \mbox{ , } \int_\Gamma c \ul{v} \cdot \ul{u}
        \mbox{ , }
        \int_\Gamma c p r \mbox{ , } \int_\Gamma c \ul{u} \cdot \ul{w}

    :Arguments 1:
        - material : :math:`c` (optional)
        - virtual  : :math:`q` or :math:`\ul{v}`
        - state    : :math:`p` or :math:`\ul{u}`

    :Arguments 2:
        - material    : :math:`c` (optional)
        - parameter_1 : :math:`p` or :math:`\ul{u}`
        - parameter_2 : :math:`r` or :math:`\ul{w}`
    """
    name = 'dw_surface_dot'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')
    integration = 'surface'

class BCNewtonTerm(DotProductSurfaceTerm):
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
        fargs = DotProductSurfaceTerm.get_fargs(self, alpha, virtual, state,
                                                mode, term_mode, diff_var,
                                                **kwargs)
        fargs = fargs[:1] + (fargs[1] - p_outer,) + fargs[2:]

        return fargs

class DotSProductVolumeOperatorWTHTerm(THTerm):
    r"""
    Fading memory volume :math:`L^2(\Omega)` weighted dot product for
    scalar fields. Can use derivatives.

    :Definition:

    .. math::
        \int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau) \difd{\tau} \right] q

    :Arguments:
        - ts       : :class:`TimeStepper` instance
        - material : :math:`\Gcal(\tau)`
        - virtual  : :math:`q`
        - state    : :math:`p`
    """
    name = 'dw_volume_dot_w_scalar_th'
    arg_types = ('ts', 'material', 'virtual', 'state')

    function = staticmethod(terms.dw_volume_dot_scalar)

    def get_fargs(self, ts, mats, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        if diff_var is None:
            def iter_kernel():
                for ii, mat in enumerate(mats):
                    val_qp = self.get(state, 'val', step=-ii)
                    mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                    yield ii, (ts.dt * mat, val_qp, vg.bf, vg.bf, vg, 0)
            fargs = iter_kernel

        else:
            val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
            mat = nm.tile(mats[0], (n_el, n_qp, 1, 1))
            fargs = ts.dt * mat, val_qp, vg.bf, vg.bf, vg, 1

        return fargs

class DotSProductVolumeOperatorWETHTerm(ETHTerm):
    r"""
    Fading memory volume :math:`L^2(\Omega)` weighted dot product for
    scalar fields. This term has the same definition as
    dw_volume_dot_w_scalar_th, but assumes an exponential approximation of
    the convolution kernel resulting in much higher efficiency. Can use
    derivatives.

    :Definition:

    .. math::
        \int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau) \difd{\tau} \right] q

    :Arguments:
        - ts         : :class:`TimeStepper` instance
        - material_0 : :math:`\Gcal(0)`
        - material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
        - virtual    : :math:`q`
        - state      : :math:`p`
    """
    name = 'dw_volume_dot_w_scalar_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'virtual', 'state')

    function = staticmethod(terms.dw_volume_dot_scalar)

    def get_fargs(self, ts, mat0, mat1, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _, key = self.get_mapping(state, return_key=True)

        if diff_var is None:
            val_qp = self.get(state, 'val')

            key += tuple(self.arg_names[ii] for ii in [1, 2, 4])
            data = self.get_eth_data(key, state, mat1, val_qp)

            fargs = (ts.dt * mat0, data.history + data.values,
                     vg.bf, vg.bf, vg, 0)

        else:
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            fargs = ts.dt * mat0, aux, vg.bf, vg.bf, vg, 1

        return fargs

class VectorDotGradScalarTerm(Term):
    r"""
    Volume dot product of a vector and a gradient of scalar.
    Can be evaluated.

    :Definition:

    .. math::
        \int_{\Omega} \ul{v} \cdot \nabla p \mbox{ , }
        \int_{\Omega} \ul{u} \cdot \nabla q \\
        \int_{\Omega} c \ul{v} \cdot \nabla p \mbox{ , }
        \int_{\Omega} c \ul{u} \cdot \nabla q \\
        \int_{\Omega} \ul{v} \cdot \ull{M} \cdot \nabla p \mbox{ , }
        \int_{\Omega} \ul{u} \cdot \ull{M} \cdot \nabla q

    :Arguments 1:
        - material : :math:`c` or :math:`\ull{M}` (optional)
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{p}`

    :Arguments 2:
        - material : :math:`c` or :math:`\ull{M}` (optional)
        - state    : :math:`\ul{u}`
        - virtual  : :math:`\ul{q}`

    :Arguments 3:
        - material    : :math:`c` or :math:`\ull{M}` (optional)
        - parameter_v : :math:`\ul{u}`
        - parameter_s : :math:`p`
    """
    name = 'dw_v_dot_grad_s'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'state', 'virtual'),
                 ('opt_material', 'parameter_v', 'parameter_s'))
    modes = ('v_weak', 's_weak', 'eval')

    def check_shapes(self, coef, vvar, svar):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)
        assert_(n_c == dim)
        assert_(svar.n_components == 1)

        if coef is not None:
            assert_((coef.shape[1:] == (n_qp, 1, 1))
                    or (coef.shape[1:] == (n_qp, dim, dim)))
            assert_((coef.shape[0] == 1) or (coef.shape[0] == n_el))

    def get_fargs(self, coef, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)
        if coef is None:
            coef = nm.ones((1, n_qp, 1, 1), dtype=nm.float64)

        if mode == 'weak':
            if self.mode == 'v_weak':
                qp_var, qp_name = svar, 'grad'

            else:
                qp_var, qp_name = vvar, 'val'

            vvg, _ = self.get_mapping(vvar)
            svg, _ = self.get_mapping(svar)

            if diff_var is None:
                val_qp = self.get(qp_var, qp_name)
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return coef, val_qp, vvg.bf, vvg, svg, fmode

        elif mode == 'eval':
            vvg, _ = self.get_mapping(vvar)

            grad = self.get(svar, 'grad')
            val = self.get(vvar, 'val')

            return coef, grad, val, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, coef, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types(self):
        self.function = {
            'v_weak' : terms.dw_v_dot_grad_s_vw,
            's_weak' : terms.dw_v_dot_grad_s_sw,
            'eval' : DotProductVolumeTerm.d_dot,
        }[self.mode]
