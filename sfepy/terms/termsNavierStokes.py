import numpy as nm

from sfepy.terms.terms import Term, terms

class DivGradTerm(Term):
    r"""
    :Description:
    Diffusion term.

    :Definition:
    .. math::
        \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}

    :Arguments:
        material : :math:`\nu` (viscosity),
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_div_grad'
    arg_types = ('material', 'virtual', 'state')

    function = staticmethod(terms.term_ns_asm_div_grad)

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if diff_var is None:
            grad = self.get(state, 'grad').transpose((0, 1, 3, 2))
            sh = grad.shape
            grad = grad.reshape((sh[0], sh[1], sh[2] * sh[3], 1))
            fmode = 0

        else:
            grad = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return grad, mat, vg, fmode

class ConvectTerm(Term):
    r"""
    :Description:
    Nonlinear convective term.

    :Definition:
    .. math::
        \int_{\Omega} ((\ul{u} \cdot \nabla) \ul{u}) \cdot \ul{v}

    :Arguments:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{u}`
    """
    name = 'dw_convect'
    arg_types = ('virtual', 'state')

    function = staticmethod(terms.term_ns_asm_convect)

    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
        val_qp = self.get(state, 'val')

        fmode = diff_var is not None

        return grad, val_qp, vg.bf, vg, fmode

class LinearConvectTerm(Term):
    r"""
    :Description:
    Linearized convective term.

    :Definition:
    .. math::
        \int_{\Omega} ((\ul{b} \cdot \nabla) \ul{u}) \cdot \ul{v}

    .. math::
        ((\ul{b} \cdot \nabla) \ul{u})|_{qp}

    :Arguments:
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dw_lin_convect'
    arg_types = ('virtual', 'parameter', 'state')

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

            return grad, val_qp, vg.bf, vg, fmode

        elif mode == 'qp':
            grad = self.get(state, 'grad').transpose((0, 1, 3, 2)).copy()
            fmode = 2

            return grad, val_qp, vg.bf, vg, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

class StokesTerm(Term):
    r"""
    :Description:
    Stokes problem coupling term. Corresponds to weak forms of gradient and
    divergence terms. Can be evaluated.

    :Definition:
    .. math::
        \int_{\Omega} p\ \nabla \cdot \ul{v} \mbox{ , } \int_{\Omega} q\ \nabla
        \cdot \ul{u}

    :Arguments 1:
        virtual : :math:`\ul{v}`,
        state   : :math:`\ul{p}`

    :Arguments 2:
        state   : :math:`\ul{u}`,
        virtual : :math:`\ul{q}`

    :Arguments 3:
        parameter_v : :math:`\ul{u}`,
        parameter_s : :math:`p`
    """
    name = 'dw_stokes'
    arg_types = (('virtual', 'state'),
                 ('state', 'virtual'),
                 ('parameter_v', 'parameter_s'))
    modes = ('grad', 'div', 'eval')

    @staticmethod
    def d_eval(out, coef, vec_qp, div, vvg):
        out_qp = coef * vec_qp * div

        status = vvg.integrate(out, out_qp)

        return status

    def get_fargs(self, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'grad':
            qp_var, qp_name = svar, 'val'

        else:
            qp_var, qp_name = vvar, 'div'

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)
        coef = kwargs.get('material')
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

            return coef, val_qp, svg.bf, vvg, fmode

        elif mode == 'eval':
            vvg, _ = self.get_mapping(vvar)

            div = self.get(vvar, 'div')
            vec_qp = self.get(svar, 'val')

            return coef, vec_qp, div, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, vvar, svar,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(vvar)

        return (n_el, 1, 1, 1), vvar.dtype

    def set_arg_types(self):
        self.function = {
            'grad' : terms.dw_grad,
            'div' : terms.dw_div,
            'eval' : self.d_eval,
        }[self.mode]

class StokesWTerm(StokesTerm):
    r"""
    :Description:
    Stokes problem coupling term weighted by a scalar function. Corresponds to
    weighted weak forms of gradient and divergence terms. Can be evaluated.

    :Definition:
    .. math::
        \int_{\Omega} c\ p\ \nabla \cdot \ul{v} \mbox{ , }
        \int_{\Omega} c\ q\ \nabla \cdot \ul{u}

    :Arguments 1:
        material : :math:`c`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`p`

    :Arguments 2:
        material : :math:`c`,
        state    : :math:`\ul{u}`,
        virtual  : :math:`q`

    :Arguments 3:
        material    : :math:`c`,
        parameter_v : :math:`\ul{u}`,
        parameter_s : :math:`p`
    """
    name = 'dw_stokes_w'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_v', 'parameter_s'))
    modes = ('grad', 'div', 'eval')

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        fargs = StokesTerm(self, virtual, state,
                           mode, term_mode, diff_var,
                           material=material, **kwargs)
        return fargs

class GradQTerm(Term):
    r"""
    :Description:
    Gradient term (weak form) in quadrature points.

    :Definition:
    .. math::
        (\nabla p)|_{qp}

    :Arguments:
        state : :math:`p`
    """
    name = 'dq_grad'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, grad):
        out[:] = grad

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return (self.get(parameter, 'grad'),)

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, n_qp, dim, 1), parameter.dtype

class GradETerm(Term):
    r"""
    :Description:
    Gradient term (weak form) in averaged in elements.

    :Definition:
    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} \nabla p /
        \int_{T_K} 1 \mbox{ or } \int_{T_K} \nabla \ul{w} /
        \int_{T_K} 1

    :Arguments:
        state : :math:`p` or :math:`\ul{w}`
    """
    name = 'de_grad'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, grad, vg, fmode):
        status = vg.integrate(out, grad, fmode)

        return status

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        grad = self.get(parameter, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return grad, vg, fmode

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, dim, 1), parameter.dtype

class DivQTerm(Term):
    r"""
    :Description:
    Divergence term (weak form) in quadrature points.

    :Definition:
    .. math::
        (\nabla \cdot \ul{u})|_{qp}

    :Arguments:
        state : :math:`\ul{u}`
    """
    name = 'dq_div'
    arg_types = ('parameter',)

    function = staticmethod(terms.dq_div_vector)

    @staticmethod
    def function(out, div):
        out[:] = div

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return (self.get(parameter, 'div'),)

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, n_qp, 1, 1), parameter.dtype

class DivETerm(Term):
    r"""
    :Description:
    Evaluate divergence term.

    :Definition:
    .. math::
         \int_{\Omega} \nabla \cdot \ul{u}

    :Arguments:
        parameter : :math:`\ul{u}`
    """
    name = 'de_div'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, div, vg, fmode):
        status = vg.integrate(out, div, fmode)

        return status

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        div = self.get(parameter, 'div')

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return div, vg, fmode

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, 1, 1), parameter.dtype

class GradDivStabilizationTerm(Term):
    r"""
    :Description:
    Grad-div stabilization term ( :math:`\gamma` is a global stabilization
    parameter).

    :Definition:
    .. math::
        \gamma \int_{\Omega} (\nabla\cdot\ul{u}) \cdot (\nabla\cdot\ul{v})

    :Arguments:
        material : :math:`\gamma`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_st_grad_div'
    arg_types = ('material', 'virtual', 'state')

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

from sfepy.terms.termsLaplace import LaplaceTerm
class PSPGPStabilizationTerm(LaplaceTerm):
    r"""
    :Description:
    PSPG stabilization term, pressure part ( :math:`\tau` is a local
    stabilization parameter), alias to Laplace term dw_laplace.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ \nabla p \cdot \nabla q

    :Arguments:
        material : :math:`\tau_K`,
        virtual  : :math:`q`,
        state    : :math:`p`
    """
    name = 'dw_st_pspg_p'

class PSPGCStabilizationTerm(Term):
    r"""
    :Description:
    PSPG stabilization term, convective part ( :math:`\tau` is a local
    stabilization parameter).

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ ((\ul{b} \cdot \nabla) \ul{u})
        \cdot \nabla q

    :Arguments:
        material  : :math:`\tau_K`,
        virtual   : :math:`q`,
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dw_st_pspg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')

    function = staticmethod(terms.dw_st_pspg_c)

    def get_fargs(self, tau, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sap, svg = self.get_approximation(virtual)
        vap, vvg = self.get_approximation(state)

        val_qp = self.get(parameter, 'val')
        conn = vap.get_connectivity(self.region, self.integration)

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return val_qp, state(), tau, svg, vvg, conn, fmode

class SUPGPStabilizationTerm(Term):
    r"""
    :Description:
    SUPG stabilization term, pressure part ( :math:`\delta` is a local
    stabilization parameter).

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ \nabla p\cdot ((\ul{b} \cdot
        \nabla) \ul{v})

    :Arguments:
        material  : :math:`\delta_K`,
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{b}`,
        state     : :math:`p`
    """
    name = 'dw_st_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')

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
    :Description:
    SUPG stabilization term, convective part ( :math:`\delta` is a local
    stabilization parameter).

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ ((\ul{b} \cdot \nabla)
        \ul{u})\cdot ((\ul{b} \cdot \nabla) \ul{v})

    :Arguments:
        material  : :math:`\delta_K`,
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{b}`,
        state     : :math:`\ul{u}`
    """
    name = 'dw_st_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')

    function = staticmethod(terms.dw_st_supg_c)

    def get_fargs(self, delta, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, vg = self.get_approximation(virtual)

        val_qp = self.get(parameter, 'val')
        conn = ap.get_connectivity(self.region, self.integration)

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return val_qp, state(), delta, vg, conn, fmode
