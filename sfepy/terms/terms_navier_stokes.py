import numpy as nm

from sfepy.linalg import dot_sequences, insert_strided_axis
from sfepy.terms.terms import Term, terms

class DivGradTerm(Term):
    r"""
    Diffusion term.

    :Definition:

    .. math::
        \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u} \mbox{ , }
        \int_{\Omega} \nabla \ul{v} : \nabla \ul{u}

    :Arguments:
        - material: :math:`\nu` (viscosity, optional)
        - virtualparameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'dw_div_grad'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', 'state'),
                   'state' : 'D', 'parameter_1' : 'D', 'parameter_2' : 'D'},
                  {'opt_material' : None}]
    modes = ('weak', 'eval')

    function = staticmethod(terms.term_ns_asm_div_grad)

    def d_div_grad(self, out, grad1, grad2, mat, vg, fmode):
        sh = grad1.shape
        g1 = grad1.reshape((sh[0], sh[1], sh[2] * sh[3]))
        sh = grad2.shape
        g2 = grad2.reshape((sh[0], sh[1], sh[2] * sh[3]))
        aux = mat * dot_sequences(g1[..., None], g2, 'ATB')[..., None]

        if fmode == 2:
            out[:] = aux
            status = 0

        else:
            status = vg.integrate(out, aux, fmode)

        return status

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgs, _ = self.get_mapping(state)
        vgv, _ = self.get_mapping(virtual)

        if mat is None:
            n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
            mat = nm.ones((1, n_qp, 1, 1), dtype=nm.float64)

        if mode == 'weak':
            if diff_var is None:
                grad = self.get(state, 'grad').transpose((0, 1, 3, 2))
                sh = grad.shape
                grad = grad.reshape((sh[0], sh[1], sh[2] * sh[3], 1))
                fmode = 0

            else:
                grad = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return grad, mat, vgv, vgs, fmode

        elif mode == 'eval':
            grad1 = self.get(virtual, 'grad')
            grad2 = self.get(state, 'grad')
            fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

            return grad1, grad2, mat, vgv, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = terms.term_ns_asm_div_grad

        else:
            self.function = self.d_div_grad

class ConvectTerm(Term):
    r"""
    Nonlinear convective term.

    :Definition:

    .. math::
        \int_{\Omega} ((\ul{u} \cdot \nabla) \ul{u}) \cdot \ul{v}

    :Arguments:
        - virtual : :math:`\ul{v}`
        - state   : :math:`\ul{u}`
    """
    name = 'dw_convect'
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D'}

    function = staticmethod(terms.term_ns_asm_convect)

    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
        val_qp = self.get(state, 'val')

        fmode = diff_var is not None

        return grad, val_qp, vg, fmode

class LinearConvectTerm(Term):
    r"""
    Linearized convective term.

    :Definition:

    .. math::
        \int_{\Omega} ((\ul{w} \cdot \nabla) \ul{u}) \cdot \ul{v}

    .. math::
        ((\ul{w} \cdot \nabla) \ul{u})|_{qp}

    :Arguments:
        - virtual   : :math:`\ul{v}`
        - parameter : :math:`\ul{w}`
        - state     : :math:`\ul{u}`
    """
    name = 'dw_lin_convect'
    arg_types = ('virtual', 'parameter', 'state')
    arg_shapes = {'virtual' : ('D', 'state'), 'parameter' : 'D', 'state' : 'D'}

    function = staticmethod(terms.dw_lin_convect)

    def get_fargs(self, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        val_qp = self.get(parameter, 'val')

        if mode == 'weak':
            if diff_var is None:
                grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
                fmode = 0

            else:
                grad = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return grad, val_qp, vg, fmode

        elif mode == 'qp':
            grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
            fmode = 2

            return grad, val_qp, vg, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class LinearConvect2Term(Term):
    r"""
    Linearized convective term with the convection velocity given as a material
    parameter.

    :Definition:

    .. math::
        \int_{\Omega} ((\ul{c} \cdot \nabla) \ul{u}) \cdot \ul{v}

    .. math::
        ((\ul{c} \cdot \nabla) \ul{u})|_{qp}

    :Arguments:
        - material : :math:`\ul{c}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_lin_convect2'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : 'D, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}

    function = staticmethod(terms.dw_lin_convect)

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mode == 'weak':
            if diff_var is None:
                grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
                fmode = 0

            else:
                grad = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return grad, material, vg, fmode

        elif mode == 'qp':
            grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
            fmode = 2

            return grad, material, vg, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class StokesTerm(Term):
    r"""
    Stokes problem coupling term. Corresponds to weak forms of gradient and
    divergence terms. Can be evaluated.

    :Definition:

    .. math::
        \int_{\Omega} p\ \nabla \cdot \ul{v} \mbox{ , }
        \int_{\Omega} q\ \nabla \cdot \ul{u}\\
        \mbox{ or }
        \int_{\Omega} c\ p\ \nabla \cdot \ul{v} \mbox{ , }
        \int_{\Omega} c\ q\ \nabla \cdot \ul{u}

    :Arguments 1:
        - material: :math:`c` (optional)
        - virtual/parameter_v: :math:`\ul{v}`
        - state/parameter_s: :math:`p`

    :Arguments 2:
        - material : :math:`c` (optional)
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
    """
    name = 'dw_stokes'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'state', 'virtual'),
                 ('opt_material', 'parameter_v', 'parameter_s'))
    arg_shapes = [{'opt_material' : '1, 1',
                   'virtual/grad' : ('D', None), 'state/grad' : 1,
                   'virtual/div' : (1, None), 'state/div' : 'D',
                   'parameter_v' : 'D', 'parameter_s' : 1},
                  {'opt_material' : None}]
    modes = ('grad', 'div', 'eval')

    @staticmethod
    def d_eval(out, coef, vec_qp, div, vvg):
        out_qp = coef * vec_qp * div

        status = vvg.integrate(out, out_qp)

        return status

    def get_fargs(self, coef, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'grad':
            qp_var, qp_name = svar, 'val'

        else:
            qp_var, qp_name = vvar, 'div'

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)
        if coef is None:
            coef = nm.ones((1, n_qp, 1, 1), dtype=nm.float64)

        if mode == 'weak':
            vvg, _ = self.get_mapping(vvar)
            svg, _ = self.get_mapping(svar)

            if diff_var is None:
                val_qp = self.get(qp_var, qp_name)
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return coef, val_qp, svg, vvg, fmode

        elif mode == 'eval':
            vvg, _ = self.get_mapping(vvar)

            div = self.get(vvar, 'div')
            vec_qp = self.get(svar, 'val')

            return coef, vec_qp, div, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, coef, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types(self):
        self.function = {
            'grad' : terms.dw_grad,
            'div' : terms.dw_div,
            'eval' : self.d_eval,
        }[self.mode]

class GradTerm(Term):
    r"""
    Evaluate gradient of a scalar or vector field.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\cal{D}} \nabla p \mbox{ or } \int_{\cal{D}} \nabla \ul{u}\\
        \int_{\cal{D}} c \nabla p \mbox{ or } \int_{\cal{D}} c \nabla \ul{u}

    :Arguments:
        - parameter : :math:`p` or :math:`\ul{u}`
    """
    name = 'ev_grad'
    arg_types = ('opt_material', 'parameter')
    arg_shapes = [{'opt_material': '1, 1', 'parameter': 'N'},
                  {'opt_material': None}]
    integration = ('cell', 'facet_extra')

    @staticmethod
    def function(out, mat, grad, vg, fmode):
        if fmode == 2:
            out[:] = grad if mat is None else mat * grad
            status = 0

        else:
            status = vg.integrate(out, grad, fmode) if mat is None else\
                     vg.integrate(out, mat * grad, fmode)

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad', integration=self.act_integration)

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return mat, grad, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, dim, n_c), parameter.dtype

class DivTerm(Term):
    r"""
    Evaluate divergence of a vector field.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
         \int_{\cal{D}} \nabla \cdot \ul{u} \mbox { , }
         \int_{\cal{D}} c \nabla \cdot \ul{u}

    :Arguments:
        - parameter : :math:`\ul{u}`
    """
    name = 'ev_div'
    arg_types = ('opt_material', 'parameter')
    arg_shapes = [{'opt_material': '1, 1', 'parameter': 'D'},
                  {'opt_material': None}]
    integration = ('cell', 'facet_extra')

    @staticmethod
    def function(out, mat, div, vg, fmode):
        if fmode == 2:
            out[:] = div if mat is None else mat * div
            status = 0

        else:
            status = vg.integrate(out, div, fmode) if mat is None else\
                     vg.integrate(out, mat * div, fmode)

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        div = self.get(parameter, 'div', integration=self.act_integration)

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return mat, div, vg, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, 1, 1), parameter.dtype

class DivOperatorTerm(Term):
    r"""
    Weighted divergence term of a test function.

    :Definition:

    .. math::
        \int_{\Omega} \nabla \cdot \ul{v} \mbox { or } \int_{\Omega} c \nabla
        \cdot \ul{v}

    :Arguments:
        - material : :math:`c` (optional)
        - virtual  : :math:`\ul{v}`
    """
    name = 'dw_div'
    arg_types = ('opt_material', 'virtual')
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : ('D', None)},
                  {'opt_material' : None}]

    @staticmethod
    def function(out, mat, vg):
        div_bf = vg.bfg

        n_el, n_qp, dim, n_ep = div_bf.shape
        div_bf = div_bf.reshape((n_el, n_qp, dim * n_ep, 1))
        div_bf = nm.ascontiguousarray(div_bf)

        if mat is not None:
            status = vg.integrate(out, mat * div_bf)
        else:
            status = vg.integrate(out, div_bf)

        return status

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        return mat, vg

class StokesWaveTerm(Term):
    r"""
    Stokes dispersion term with the wave vector :math:`\ul{\kappa}`.

    :Definition:

    .. math::
        \int_{\Omega} (\ul{\kappa} \cdot \ul{v}) (\ul{\kappa} \cdot \ul{u})

    :Arguments:
        - material : :math:`\ul{\kappa}`
        - virtual  : :math:`\ul{v}`
        - statee   : :math:`\ul{u}`
    """
    name = 'dw_stokes_wave'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '.: D',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']

    @staticmethod
    def function(out, out_qp, geo, fmode):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, kappa, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import create_adof_conn, expand_basis

        geo, _ = self.get_mapping(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        ebf = expand_basis(geo.bf, dim)

        aux = nm.einsum('i,...ij->...j', kappa, ebf)[0, :, None, :]
        kebf = insert_strided_axis(aux, 0, n_el)

        if diff_var is None:
            econn = state.field.get_econn('cell', self.region)
            adc = create_adof_conn(nm.arange(state.n_dof, dtype=nm.int32),
                                   econn, n_c, 0)
            vals = state()[adc]
            aux = dot_sequences(kebf, vals[:, None, :, None])
            out_qp = dot_sequences(kebf, aux, 'ATB')
            fmode = 0

        else:
            out_qp = dot_sequences(kebf, kebf, 'ATB')
            fmode = 1

        return out_qp, geo, fmode

class StokesWaveDivTerm(Term):
    r"""
    Stokes dispersion term with the wave vector :math:`\ul{\kappa}` and the
    divergence operator.

    :Definition:

    .. math::
        \int_{\Omega} (\ul{\kappa} \cdot \ul{v}) (\nabla \cdot \ul{u}) \;,
        \int_{\Omega} (\ul{\kappa} \cdot \ul{u}) (\nabla \cdot \ul{v})

    :Arguments 1:
        - material : :math:`\ul{\kappa}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`

    :Arguments 2:
        - material : :math:`\ul{\kappa}`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`\ul{v}`
    """
    name = 'dw_stokes_wave_div'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'))
    arg_shapes = {'material' : '.: D',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['2_3', '2_4', '3_4', '3_8']
    modes = ('kd', 'dk')

    @staticmethod
    def function(out, out_qp, geo, fmode):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, kappa, kvar, dvar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.discrete.variables import create_adof_conn, expand_basis

        geo, _ = self.get_mapping(dvar)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(kvar)

        ebf = expand_basis(geo.bf, dim)

        aux = nm.einsum('i,...ij->...j', kappa, ebf)[0, :, None, :]
        kebf = insert_strided_axis(aux, 0, n_el)

        div_bf = geo.bfg

        div_bf = div_bf.reshape((n_el, n_qp, 1, dim * n_en))
        div_bf = nm.ascontiguousarray(div_bf)

        if diff_var is None:
            avar = dvar if self.mode == 'kd' else kvar
            econn = avar.field.get_econn('cell', self.region)
            adc = create_adof_conn(nm.arange(avar.n_dof, dtype=nm.int32),
                                   econn, n_c, 0)
            vals = avar()[adc]

            if self.mode == 'kd':
                aux = dot_sequences(div_bf, vals[:, None, :, None])
                out_qp = dot_sequences(kebf, aux, 'ATB')

            else:
                aux = dot_sequences(kebf, vals[:, None, :, None])
                out_qp = dot_sequences(div_bf, aux, 'ATB')

            fmode = 0

        else:
            if self.mode == 'kd':
                out_qp = dot_sequences(kebf, div_bf, 'ATB')

            else:
                out_qp = dot_sequences(div_bf, kebf, 'ATB')

            fmode = 1

        return out_qp, geo, fmode

class GradDivStabilizationTerm(Term):
    r"""
    Grad-div stabilization term ( :math:`\gamma` is a global stabilization
    parameter).

    :Definition:

    .. math::
        \gamma \int_{\Omega} (\nabla\cdot\ul{u}) \cdot (\nabla\cdot\ul{v})

    :Arguments:
        - material : :math:`\gamma`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_st_grad_div'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'state' : 'D'}

    function = staticmethod(terms.dw_st_grad_div)

    def get_fargs(self, gamma, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if diff_var is None:
            div = self.get(state, 'div')
            fmode = 0

        else:
            div = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return div, gamma, vg, fmode

from sfepy.terms.terms_diffusion import LaplaceTerm
class PSPGPStabilizationTerm(LaplaceTerm):
    r"""
    PSPG stabilization term, pressure part ( :math:`\tau` is a local
    stabilization parameter), alias to Laplace term dw_laplace.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ \nabla p \cdot \nabla q

    :Arguments:
        - material : :math:`\tau_K`
        - virtual  : :math:`q`
        - state    : :math:`p`
    """
    name = 'dw_st_pspg_p'

class PSPGCStabilizationTerm(Term):
    r"""
    PSPG stabilization term, convective part ( :math:`\tau` is a local
    stabilization parameter).

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ ((\ul{b} \cdot \nabla) \ul{u})
        \cdot \nabla q

    :Arguments:
        - material  : :math:`\tau_K`
        - virtual   : :math:`q`
        - parameter : :math:`\ul{b}`
        - state     : :math:`\ul{u}`
    """
    name = 'dw_st_pspg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : (1, None),
                  'parameter' : 'D', 'state' : 'D'}

    function = staticmethod(terms.dw_st_pspg_c)

    def get_fargs(self, tau, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        svg, _ = self.get_mapping(virtual)
        vvg, _ = self.get_mapping(state)

        val_qp = self.get(parameter, 'val')
        conn = state.field.get_econn(self.act_integration, self.region)

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return val_qp, state(), tau, svg, vvg, conn, fmode

class SUPGPStabilizationTerm(Term):
    r"""
    SUPG stabilization term, pressure part ( :math:`\delta` is a local
    stabilization parameter).

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ \nabla p\cdot ((\ul{b} \cdot
        \nabla) \ul{v})

    :Arguments:
        - material  : :math:`\delta_K`
        - virtual   : :math:`\ul{v}`
        - parameter : :math:`\ul{b}`
        - state     : :math:`p`
    """
    name = 'dw_st_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', None),
                  'parameter' : 'D', 'state' : 1}

    function = staticmethod(terms.dw_st_supg_p)

    def get_fargs(self, delta, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vvg, _ = self.get_mapping(virtual)
        svg, _ = self.get_mapping(state)

        val_qp = self.get(parameter, 'val')

        if diff_var is None:
            grad = self.get(state, 'grad')
            fmode = 0

        else:
            grad = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return val_qp, grad, delta, vvg, svg, fmode

class SUPGCStabilizationTerm(Term):
    r"""
    SUPG stabilization term, convective part ( :math:`\delta` is a local
    stabilization parameter).

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ ((\ul{b} \cdot \nabla)
        \ul{u})\cdot ((\ul{b} \cdot \nabla) \ul{v})

    :Arguments:
        - material  : :math:`\delta_K`
        - virtual   : :math:`\ul{v}`
        - parameter : :math:`\ul{b}`
        - state     : :math:`\ul{u}`
    """
    name = 'dw_st_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'parameter' : 'D', 'state' : 'D'}

    function = staticmethod(terms.dw_st_supg_c)

    def get_fargs(self, delta, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(virtual)

        val_qp = self.get(parameter, 'val')
        conn = virtual.field.get_econn(self.act_integration, self.region)

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return val_qp, state(), delta, vg, conn, fmode
