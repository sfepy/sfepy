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
        \int_\Omega c p r \mbox{ , } \int_\Omega c \ul{u} \cdot \ul{w} \\
        \int_\Omega \ul{v} \cdot \ull{M} \cdot \ul{u}
        \mbox{ , }
        \int_\Omega \ul{u} \cdot \ull{M} \cdot \ul{w}

    :Arguments 1:
        - material : :math:`c` or :math:`\ull{M}` (optional)
        - virtual  : :math:`q` or :math:`\ul{v}`
        - state    : :math:`p` or :math:`\ul{u}`

    :Arguments 2:
        - material    : :math:`c` or :math:`\ull{M}` (optional)
        - parameter_1 : :math:`p` or :math:`\ul{u}`
        - parameter_2 : :math:`r` or :math:`\ul{w}`
    """
    name = 'dw_volume_dot'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, 'state'),
                   'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1},
                  {'opt_material' : None},
                  {'opt_material' : '1, 1', 'virtual' : ('D', 'state'),
                   'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : 'D, D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    @staticmethod
    def dw_dot(out, mat, val_qp, vgeo, sgeo, fun, fmode):
        status = fun(out, mat, val_qp, vgeo, sgeo, fmode)
        return status

    @staticmethod
    def d_dot(out, mat, val1_qp, val2_qp, geo):
        if mat is None:
            if val1_qp.shape[2] > 1:
                if val2_qp.shape[2] == 1:
                    aux = dot_sequences(val1_qp, geo.normal, mode='ATB')
                    vec = dot_sequences(aux, val2_qp, mode='AB')

                else:
                    vec = dot_sequences(val1_qp, val2_qp, mode='ATB')

            else:
                vec = val1_qp * val2_qp

        elif mat.shape[-1] == 1:
            if val1_qp.shape[2] > 1:
                vec = mat * dot_sequences(val1_qp, val2_qp, mode='ATB')

            else:
                vec = mat * val1_qp * val2_qp

        else:
            aux = dot_sequences(mat, val2_qp, mode='AB')
            vec = dot_sequences(val1_qp, aux, mode='ATB')

        status = geo.integrate(out, vec)

        return status

    def check_shapes(self, mat, virtual, state):
        is_vector_scalar = ((virtual.n_components == 1)
            and (state.n_components == state.dim))\
            or ((virtual.n_components == virtual.dim)
            and (state.n_components == 1))

        assert_((virtual.n_components == state.n_components)
                or is_vector_scalar)

        if mat is not None:
            n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
            assert_((mat.shape[1:] == (n_qp, 1, 1))
                    or ((mat.shape[1:] == (n_qp, dim, dim)) and (n_c == dim)))
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
                if ((self.integration == 'volume')
                    or (virtual.n_components > 1)):
                    fun = terms.dw_volume_dot_vector

                else:
                    fun = terms.dw_surface_s_v_dot_n

            else:
                if ((self.integration == 'volume')
                    or (virtual.n_components == 1)):
                    fun = terms.dw_volume_dot_scalar

                else:
                    fun = terms.dw_surface_v_dot_n_s

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
        \int_\Gamma \ul{v} \cdot \ul{n} p \mbox{ , }
        \int_\Gamma q \ul{n} \cdot \ul{u} \mbox{ , }
        \int_\Gamma p r \mbox{ , } \int_\Gamma \ul{u} \cdot \ul{w}
        \mbox{ , } \int_\Gamma \ul{w} \cdot \ul{n} p \\
        \int_\Gamma c q p \mbox{ , } \int_\Gamma c \ul{v} \cdot \ul{u}
        \mbox{ , }
        \int_\Gamma c p r \mbox{ , } \int_\Gamma c \ul{u} \cdot \ul{w} \\
        \int_\Gamma \ul{v} \cdot \ull{M} \cdot \ul{u}
        \mbox{ , }
        \int_\Gamma \ul{u} \cdot \ull{M} \cdot \ul{w}

    :Arguments 1:
        - material : :math:`c` or :math:`\ull{M}` (optional)
        - virtual  : :math:`q` or :math:`\ul{v}`
        - state    : :math:`p` or :math:`\ul{u}`

    :Arguments 2:
        - material    : :math:`c` or :math:`\ull{M}` (optional)
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
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'virtual' : (1, 'state'), 'state' : 1}
    mode = 'weak'

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
                    yield ii, (ts.dt * mat, val_qp, vg, vg, 0)
            fargs = iter_kernel

        else:
            val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
            mat = nm.tile(mats[0], (n_el, n_qp, 1, 1))
            fargs = ts.dt * mat, val_qp, vg, vg, 1

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

            fargs = (ts.dt * mat0, data.history + data.values, vg, vg, 0)

        else:
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            fargs = ts.dt * mat0, aux, vg, vg, 1

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
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`c` or :math:`\ull{M}` (optional)
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`

    :Arguments 3:
        - material    : :math:`c` or :math:`\ull{M}` (optional)
        - parameter_v : :math:`\ul{u}`
        - parameter_s : :math:`p`
    """
    name = 'dw_v_dot_grad_s'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'state', 'virtual'),
                 ('opt_material', 'parameter_v', 'parameter_s'))
    arg_shapes = [{'opt_material' : '1, 1',
                   'virtual/v_weak' : ('D', None), 'state/v_weak' : 1,
                   'virtual/s_weak' : (1, None), 'state/s_weak' : 'D',
                   'parameter_v' : 'D', 'parameter_s' : 1},
                  {'opt_material' : 'D, D'},
                  {'opt_material' : None}]
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

            return coef, val_qp, vvg, svg, fmode

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

class VectorDotScalarTerm(Term):
    r"""
    Volume dot product of a vector and a scalar.
    Can be evaluated.

    :Definition:

    .. math::
        \int_{\Omega} \ul{v} \cdot \ul{m} p \mbox{ , }
        \int_{\Omega} \ul{u} \cdot \ul{m} q\\

    :Arguments 1:
        - material : :math:`\ul{m}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`\ul{m}`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`

    :Arguments 3:
        - material    : :math:`\ul{m}`
        - parameter_v : :math:`\ul{u}`
        - parameter_s : :math:`p`
    """
    name = 'dw_vm_dot_s'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    arg_shapes = [{'material' : 'D, 1',
                   'virtual/v_weak' : ('D', None), 'state/v_weak' : 1,
                   'virtual/s_weak' : (1, None), 'state/s_weak' : 'D',
                   'parameter_v' : 'D', 'parameter_s' : 1}]
    modes = ('v_weak', 's_weak', 'eval')

    @staticmethod
    def dw_dot(out, mat, val_qp, bfve, bfsc, geo, fmode):

        nel, nqp, dim, nc = mat.shape
        nen = bfve.shape[2]

        status1 = 0
        if fmode in [0, 1, 3]:
            aux = nm.zeros((nel, nqp, dim * nen, nc), dtype=nm.float64)
            status1 = terms.actBfT(aux, bfve, mat)

        if fmode == 0:
            status2 = terms.mulAB_integrate(out, aux, val_qp, geo, 'AB')

        if fmode == 1:
            status2 = terms.mulAB_integrate(out, aux, bfsc, geo, 'AB')

        if fmode == 2:
            aux = (bfsc * dot_sequences(mat, val_qp,
                                        mode='ATB')).transpose((0,1,3,2))
            status2 = geo.integrate(out, nm.ascontiguousarray(aux))

        if fmode == 3:
            status2 = terms.mulAB_integrate(out, bfsc, aux, geo, 'ATBT')

        return status1 and status2

    @staticmethod
    def d_dot(out, mat, val1_qp, val2_qp, geo):
        v1, v2 = (val1_qp, val2_qp) if val1_qp.shape[2] > 1 \
                 else (val2_qp, val1_qp)
        aux = dot_sequences(v1, mat, mode='ATB')
        vec = dot_sequences(aux, v2, mode='AB')
        status = geo.integrate(out, vec)

        return status

    def check_shapes(self, coef, vvar, svar):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)
        assert_(n_c == dim)
        assert_(svar.n_components == 1)
        assert_((coef.shape[1:] == (n_qp, dim, 1))
                or (coef.shape[1:] == (n_qp, 1, dim)))
        assert_((coef.shape[0] == 1) or (coef.shape[0] == n_el))

    def get_fargs(self, coef, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        coef = coef.reshape(coef.shape[:2] + (dim, 1))

        if mode == 'weak':
            apv, vgv = self.get_approximation(vvar)
            aps, vgs = self.get_approximation(svar)

            bfve = apv.get_base('v', 0, self.integral)
            bfsc = aps.get_base('v', 0, self.integral)

            if self.mode == 'v_weak':
                qp_var, geo, fmode = svar, vgv, 0

            else:
                qp_var, geo, fmode = vvar, vgs, 2
                bfve, bfsc = bfsc, bfve

            if diff_var is None:
                val_qp = self.get(qp_var, 'val')

            else:
                val_qp = (nm.array([0], ndmin=4, dtype=nm.float64), 1)
                fmode += 1

            return coef, val_qp, bfve, bfsc, geo, fmode

        elif mode == 'eval':
            vvg, _ = self.get_mapping(vvar)
            vals = self.get(svar, 'val')
            valv = self.get(vvar, 'val')

            return coef, vals, valv, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, coef, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types(self):
        self.function = {
            'v_weak' : self.dw_dot,
            's_weak' : self.dw_dot,
            'eval' : self.d_dot,
        }[self.mode]

class ScalarDotGradIScalarTerm(Term):
    r"""
    Dot product of a scalar and the :math:`i`-th component of gradient of a
    scalar. The index should be given as a 'special_constant' material
    parameter.

    :Definition:

    .. math::
        Z^i = \int_{\Omega} q \nabla_i p

    :Arguments:
        - material : :math:`i`
        - virtual  : :math:`q`
        - state    : :math:`p`
    """
    name = 'dw_s_dot_grad_i_s'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : (1, 'state'), 'state' : 1}

    @staticmethod
    def dw_fun(out, bf, vg, grad, idx, fmode):
        cc = nm.ascontiguousarray
        bft = cc(nm.tile(bf, (out.shape[0], 1, 1, 1)))

        if fmode == 0:
            status = terms.mulAB_integrate(out, bft,
                                           cc(grad[..., idx:idx+1, :]), vg,
                                           mode='ATB')

        else:
            status = terms.mulAB_integrate(out, bft,
                                           cc(vg.bfg[:,:,idx:(idx + 1),:]), vg,
                                           mode='ATB')

        return status

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if mode == 'weak':
            if diff_var is None:
                grad = self.get(state, 'grad')
                fmode = 0

            else:
                grad = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            ap, vg = self.get_approximation(virtual)
            aps, vgs = self.get_approximation(state)

            bf = aps.get_base('v', 0, self.integral)
            idx = int(material[0, 0, 0, 0])

            return bf, vg, grad, idx, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def set_arg_types(self):
        self.function = self.dw_fun
