import numpy as nm

from sfepy.terms.terms import Term, terms
from sfepy.base.base import get_default

def grad_as_vector(grad):
    grad = grad.transpose((0, 1, 3, 2))
    sh = grad.shape
    return grad.reshape((sh[0], sh[1], sh[2] * sh[3], 1))

class AdjDivGradTerm(Term):
    r"""
    Gateaux differential of :math:`\Psi(\ul{u}) = \int_{\Omega} \nu\
    \nabla \ul{v} : \nabla \ul{u}` w.r.t. :math:`\ul{u}` in the direction
    :math:`\ul{v}` or adjoint term to `dw_div_grad`.

    :Definition:

    .. math::
        w \delta_{u} \Psi(\ul{u}) \circ \ul{v}

    :Arguments:
        - material_1 : :math:`w` (weight)
        - material_2 : :math:`\nu` (viscosity)
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_adj_div_grad'
    arg_types = ('material_1', 'material_2', 'virtual', 'parameter')
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'virtual' : ('D', None), 'parameter' : 'D'}

    function = staticmethod(terms.term_ns_asm_div_grad)

    def get_fargs(self, mat1, mat2, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgv, _ = self.get_mapping(virtual)
        vgs, _ = self.get_mapping(state)

        if diff_var is None:
            grad = grad_as_vector(self.get(state, 'grad'))
            fmode = 0

        else:
            grad = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return grad, mat1 * mat2, vgv, vgs, fmode

class AdjConvect1Term(Term):
    r"""
    The first adjoint term to nonlinear convective term `dw_convect`.

    :Definition:

    .. math::
        \int_{\Omega} ((\ul{v} \cdot \nabla) \ul{u}) \cdot \ul{w}

    :Arguments:
        - virtual   : :math:`\ul{v}`
        - state     : :math:`\ul{w}`
        - parameter : :math:`\ul{u}`
    """
    name = 'dw_adj_convect1'
    arg_types = ('virtual', 'state', 'parameter' )
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D', 'parameter' : 'D'}

    function = staticmethod(terms.dw_adj_convect1)

    def get_fargs(self, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        val_w = self.get(state, 'val')
        grad_u = self.get(parameter, 'grad') # No transposition here!

        fmode = diff_var is not None

        return val_w, grad_u, vg, fmode

class AdjConvect2Term(Term):
    r"""
    The second adjoint term to nonlinear convective term `dw_convect`.

    :Definition:

    .. math::
        \int_{\Omega} ((\ul{u} \cdot \nabla) \ul{v}) \cdot \ul{w}

    :Arguments:
        - virtual   : :math:`\ul{v}`
        - state     : :math:`\ul{w}`
        - parameter : :math:`\ul{u}`
    """
    name = 'dw_adj_convect2'
    arg_types = ('virtual', 'state', 'parameter' )
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D', 'parameter' : 'D'}

    function = staticmethod(terms.dw_adj_convect2)

    def get_fargs(self, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        val_w = self.get(state, 'val')
        val_u = self.get(parameter, 'val')

        fmode = diff_var is not None

        return val_w, val_u, vg, fmode

class SUPGCAdjStabilizationTerm(Term):
    r"""
    Adjoint term to SUPG stabilization term `dw_st_supg_c`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ [ ((\ul{v} \cdot \nabla)
        \ul{u}) ((\ul{u} \cdot \nabla) \ul{w}) + ((\ul{u} \cdot \nabla)
        \ul{u}) ((\ul{v} \cdot \nabla) \ul{w}) ]

    :Arguments:
        - material  : :math:`\delta_K`
        - virtual   : :math:`\ul{v}`
        - state     : :math:`\ul{w}`
        - parameter : :math:`\ul{u}`
    """
    name = 'dw_st_adj_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'state' : 'D', 'parameter' : 'D'}

    function = staticmethod(terms.dw_st_adj_supg_c)

    def get_fargs(self, mat, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        val_u = self.get(parameter, 'val')
        grad_u = self.get(parameter, 'grad').transpose((0, 1, 3, 2)).copy()
        conn = state.field.get_econn(self.act_integration, self.region)

        fmode = diff_var is not None

        return state(), val_u, grad_u, mat, vg, conn, fmode

class SUPGPAdj1StabilizationTerm(Term):
    r"""
    The first adjoint term to SUPG stabilization term `dw_st_supg_p`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ \nabla p (\ul{v} \cdot
        \nabla \ul{w})

    :Arguments:
        - material  : :math:`\delta_K`
        - virtual   : :math:`\ul{v}`
        - state     : :math:`\ul{w}`
        - parameter : :math:`p`
    """
    name = 'dw_st_adj1_supg_p'
    arg_types = ('material', 'virtual', 'state', 'parameter')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'state' : 'D', 'parameter' : 1}

    function = staticmethod(terms.dw_st_adj1_supg_p)

    def get_fargs(self, mat, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg_w, _ = self.get_mapping(state)

        grad_p = self.get(parameter, 'grad')
        conn_w = state.field.get_econn(self.act_integration, self.region)

        fmode = diff_var is not None

        return state(), grad_p, mat, vg_w, conn_w, fmode

class SUPGPAdj2StabilizationTerm(Term):
    r"""
    The second adjoint term to SUPG stabilization term `dw_st_supg_p`
    as well as adjoint term to PSPG stabilization term `dw_st_pspg_c`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ \nabla r (\ul{v} \cdot \nabla
        \ul{u})

    :Arguments:
        - material  : :math:`\tau_K`
        - virtual   : :math:`\ul{v}`
        - parameter : :math:`\ul{u}`
        - state     : :math:`r`
    """
    name = 'dw_st_adj2_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'state' : 1, 'parameter' : 'D'}

    function = staticmethod(terms.dw_st_adj2_supg_p)

    def get_fargs(self, mat, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg_r, _ = self.get_mapping(state)
        vg_u, _ = self.get_mapping(parameter)

        grad_u = self.get(parameter, 'grad').transpose((0, 1, 3, 2)).copy()
        conn_r = state.field.get_econn(self.act_integration, self.region)

        fmode = diff_var is not None

        return grad_u, state(), mat, vg_u, vg_r, conn_r, fmode

class SDDotTerm(Term):
    r"""
    Sensitivity (shape derivative) of dot product of scalars or vectors.

    :Definition:

    .. math::
        \int_{\Omega} p q (\nabla \cdot \ul{\Vcal}) \mbox{ , }
        \int_{\Omega} (\ul{u} \cdot \ul{w}) (\nabla \cdot \ul{\Vcal})

    :Arguments:
        - parameter_1 : :math:`p` or :math:`\ul{u}`
        - parameter_2 : :math:`q` or :math:`\ul{w}`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_dot'
    arg_types = ('parameter_1', 'parameter_2', 'parameter_mv')
    arg_shapes = [{'parameter_1' : 'D', 'parameter_2' : 'D',
                   'parameter_mv' : 'D'},
                  {'parameter_1' : 1, 'parameter_2' : 1}]

    function = staticmethod(terms.d_sd_volume_dot)

    def get_fargs(self, par1, par2, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par1)

        val1 = self.get(par1, 'val')
        val2 = self.get(par2, 'val')
        div_mv = self.get(par_mv, 'div')

        return val1, val2, div_mv, vg, get_default(term_mode, 1)

    def get_eval_shape(self, par1, par2, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par1)
        return (n_el, 1, 1, 1), par1.dtype

class SDDivTerm(Term):
    r"""
    Sensitivity (shape derivative) of Stokes term `dw_stokes` in 'div' mode.

    Supports the following term modes: 1 (sensitivity) or 0 (original term
    value).

    :Definition:

    .. math::
        \int_{\Omega} p [ (\nabla \cdot \ul{w}) (\nabla \cdot \ul{\Vcal})
        - \pdiff{\Vcal_k}{x_i} \pdiff{w_i}{x_k} ]

    :Arguments:
        - parameter_u : :math:`\ul{u}`
        - parameter_p : :math:`p`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_div'
    arg_types = ('parameter_u', 'parameter_p', 'parameter_mv')
    arg_shapes = {'parameter_u' : 'D', 'parameter_p' : 1,
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_div)

    def get_fargs(self, par_u, par_p, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        div_u = self.get(par_u, 'div')
        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        val_p = self.get(par_p, 'val')
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        return (div_u, grad_u, val_p, div_mv, grad_mv, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, par_u, par_p, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class SDDivGradTerm(Term):
    r"""
    Sensitivity (shape derivative) of diffusion term `dw_div_grad`.

    Supports the following term modes: 1 (sensitivity) or 0 (original term
    value).

    :Definition:

    .. math::
        \int_{\Omega} \hat{I} \nabla \ul{v} : \nabla \ul{u} \mbox{ , }
        \int_{\Omega} \nu \hat{I}  \nabla \ul{v} : \nabla \ul{u}

    .. math::
        \hat{I}_{ijkl} =
            \delta_{ik}\delta_{jl} \nabla \cdot \ul{\Vcal}
          - \delta_{ik}\delta_{js} {\partial \Vcal_l \over \partial x_s}
          - \delta_{is}\delta_{jl} {\partial \Vcal_k \over \partial x_s}

    :Arguments:
        - material  : :math:`\nu` (viscosity, optional)
        - parameter_u : :math:`\ul{u}`
        - parameter_w : :math:`\ul{w}`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_div_grad'
    arg_types = ('opt_material', 'parameter_u', 'parameter_w',
                 'parameter_mv')
    arg_shapes = [{'opt_material' : '1, 1',
                   'parameter_u' : 'D', 'parameter_w' : 'D',
                   'parameter_mv' : 'D'},
                  {'opt_material' : None}]

    function = staticmethod(terms.d_sd_div_grad)

    def get_fargs(self, mat, par_u, par_w, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        grad_w = grad_as_vector(self.get(par_w, 'grad'))
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        if mat is None:
            mat = nm.ones((1, div_mv.shape[1], 1, 1), dtype=nm.float64)

        return (grad_u, grad_w, div_mv, grad_mv, mat, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, mat, par_u, par_w, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class SDConvectTerm(Term):
    r"""
    Sensitivity (shape derivative) of convective term `dw_convect`.

    Supports the following term modes: 1 (sensitivity) or 0 (original term
    value).

    :Definition:

    .. math::
        \int_{\Omega} [ u_k \pdiff{u_i}{x_k} w_i (\nabla \cdot \Vcal)
        - u_k \pdiff{\Vcal_j}{x_k} \pdiff{u_i}{x_j} w_i ]

    :Arguments:
        - parameter_u : :math:`\ul{u}`
        - parameter_w : :math:`\ul{w}`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_convect'
    arg_types = ('parameter_u', 'parameter_w', 'parameter_mv')
    arg_shapes = {'parameter_u' : 'D', 'parameter_w' : 'D',
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_convect)

    def get_fargs(self, par_u, par_w, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        val_u = self.get(par_u, 'val')
        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        val_w = self.get(par_w, 'val')
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        return (val_u, grad_u, val_w, div_mv, grad_mv, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, par_u, par_w, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class NSOFMinGradTerm(Term):
    name = 'd_of_ns_min_grad'
    arg_types = ('material_1', 'material_2', 'parameter')
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'parameter' : 1}

    function = staticmethod(terms.d_of_nsMinGrad)

    def get_fargs(self, weight, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        grad = grad_as_vector(self.get(parameter, 'grad'))

        return grad, weight * mat, vg

    def get_eval_shape(self, weight, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        return (1, 1, 1, 1), parameter.dtype

class NSOFSurfMinDPressTerm(Term):
    r"""
    Sensitivity of :math:`\Psi(p)`.

    :Definition:

    .. math::
        \delta \Psi(p) = \delta \left( \int_{\Gamma_{in}}p -
        \int_{\Gamma_{out}}bpress \right)

    :Arguments:
        - material_1 : :math:`w` (weight)
        - material_2 : :math:`bpress` (given pressure)
        - parameter  : :math:`p`
    """
    name = 'ev_of_ns_surf_min_d_press'
    arg_types = ('material_1', 'material_2', 'parameter')
    arg_shapes = {'material_1' : 1, 'material_2' : 1,
                  'parameter' : 1}
    integration = 'facet'

    function = staticmethod(terms.d_of_nsSurfMinDPress)

    def get_fargs(self, weight, bpress, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        val_p = self.get(parameter, 'val')

        return val_p, weight, bpress, sg, 0

    def get_eval_shape(self, weight, bpress, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        return (1, 1, 1, 1), parameter.dtype

class NSOFSurfMinDPressDiffTerm(NSOFSurfMinDPressTerm):
    r"""
    Gateaux differential of :math:`\Psi(p)` w.r.t. :math:`p` in the
    direction :math:`q`.

    :Definition:

    .. math::
        w \delta_{p} \Psi(p) \circ q

    :Arguments:
        - material : :math:`w` (weight)
        - virtual  : :math:`q`
    """
    name = 'dw_of_ns_surf_min_d_press_diff'
    arg_types = ('material', 'virtual')
    arg_shapes = {'material' : 1, 'virtual' : (1, None)}

    def get_fargs(self, weight, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        aux = nm.array([0], ndmin=4, dtype=nm.float64)

        return aux, weight, 0.0, sg, 1

class SDGradDivStabilizationTerm(Term):
    r"""
    Sensitivity (shape derivative) of stabilization term `dw_st_grad_div`.

    :Definition:

    .. math::
        \gamma \int_{\Omega} [ (\nabla \cdot \ul{u}) (\nabla \cdot \ul{w})
        (\nabla \cdot \ul{\Vcal})
        - \pdiff{u_i}{x_k} \pdiff{\Vcal_k}{x_i} (\nabla \cdot \ul{w})
        - (\nabla \cdot \ul{u}) \pdiff{w_i}{x_k} \pdiff{\Vcal_k}{x_i} ]

    :Arguments:
        - material    : :math:`\gamma`
        - parameter_u : :math:`\ul{u}`
        - parameter_w : :math:`\ul{w}`
        - parameter_mv : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'ev_sd_st_grad_div'
    arg_types = ('material', 'parameter_u', 'parameter_w',
                 'parameter_mv')
    arg_shapes = {'material' : '1, 1',
                  'parameter_u' : 'D', 'parameter_w' : 'D',
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_st_grad_div)

    def get_fargs(self, mat, par_u, par_w, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        div_u = self.get(par_u, 'div')
        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        div_w = self.get(par_w, 'div')
        grad_w = grad_as_vector(self.get(par_w, 'grad'))
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        return (div_u, grad_u, div_w, grad_w, div_mv, grad_mv, mat, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, mat, par_u, par_w, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class SDSUPGCStabilizationTerm(Term):
    r"""
    Sensitivity (shape derivative) of stabilization term `dw_st_supg_c`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ [ (\ul{b} \cdot \nabla u_k)
        (\ul{b} \cdot \nabla w_k) (\nabla \cdot \Vcal) -
        (\ul{b} \cdot \nabla \Vcal_i) \pdiff{u_k}{x_i}
        (\ul{b} \cdot \nabla w_k) - (\ul{u} \cdot \nabla u_k)
        (\ul{b} \cdot \nabla \Vcal_i) \pdiff{w_k}{x_i} ]

    :Arguments:
        - material    : :math:`\delta_K`
        - parameter_b : :math:`\ul{b}`
        - parameter_u : :math:`\ul{u}`
        - parameter_w : :math:`\ul{w}`
        - parameter_mv : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'ev_sd_st_supg_c'
    arg_types = ('material', 'parameter_b', 'parameter_u', 'parameter_w',
                'parameter_mv')
    arg_shapes = {'material' : '1, 1',
                  'parameter_b' : 'D', 'parameter_u' : 'D', 'parameter_w' : 'D',
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_st_supg_c)

    def get_fargs(self, mat, par_b, par_u, par_w, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        val_b = self.get(par_b, 'val')
        grad_u = self.get(par_u, 'grad').transpose((0, 1, 3, 2)).copy()
        grad_w = self.get(par_w, 'grad').transpose((0, 1, 3, 2)).copy()
        div_mv = self.get(par_mv, 'div')
        grad_mv = self.get(par_mv, 'grad').transpose((0, 1, 3, 2)).copy()

        return (val_b, grad_u, grad_w, div_mv, grad_mv, mat, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, mat, par_b, par_u, par_w, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class SDPSPGCStabilizationTerm(Term):
    r"""
    Sensitivity (shape derivative) of stabilization terms `dw_st_supg_p` or
    `dw_st_pspg_c`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\
        [ \pdiff{r}{x_i} (\ul{b} \cdot \nabla u_i) (\nabla \cdot \Vcal) -
        \pdiff{r}{x_k} \pdiff{\Vcal_k}{x_i} (\ul{b} \cdot \nabla u_i)
        - \pdiff{r}{x_k} (\ul{b} \cdot \nabla \Vcal_k) \pdiff{u_i}{x_k} ]

    :Arguments:
        - material    : :math:`\delta_K`
        - parameter_b : :math:`\ul{b}`
        - parameter_u : :math:`\ul{u}`
        - parameter_r : :math:`r`
        - parameter_mv : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'ev_sd_st_pspg_c'
    arg_types = ('material', 'parameter_b', 'parameter_u', 'parameter_r',
                'parameter_mv')
    arg_shapes = {'material' : '1, 1',
                  'parameter_b' : 'D', 'parameter_u' : 'D', 'parameter_r' : 1,
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_st_pspg_c)

    def get_fargs(self, mat, par_b, par_u, par_r, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        val_b = self.get(par_b, 'val')
        grad_u = self.get(par_u, 'grad').transpose((0, 1, 3, 2)).copy()
        grad_r = self.get(par_r, 'grad')
        div_mv = self.get(par_mv, 'div')
        grad_mv = self.get(par_mv, 'grad').transpose((0, 1, 3, 2)).copy()

        return (val_b, grad_u, grad_r, div_mv, grad_mv, mat, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, mat, par_b, par_u, par_r, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class SDPSPGPStabilizationTerm(Term):
    r"""
    Sensitivity (shape derivative) of stabilization term `dw_st_pspg_p`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ [ (\nabla r \cdot \nabla p)
        (\nabla \cdot \Vcal) - \pdiff{r}{x_k} (\nabla \Vcal_k \cdot \nabla p) -
        (\nabla r \cdot \nabla \Vcal_k) \pdiff{p}{x_k} ]

    :Arguments:
        - material    : :math:`\tau_K`
        - parameter_r : :math:`r`
        - parameter_p : :math:`p`
        - parameter_mv : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'ev_sd_st_pspg_p'
    arg_types = ('material', 'parameter_r', 'parameter_p',
                'parameter_mv')
    arg_shapes = {'material' : '1, 1',
                  'parameter_r' : 1, 'parameter_p' : 1,
                  'parameter_mv' : 'D'}

    function = staticmethod(terms.d_sd_st_pspg_p)

    def get_fargs(self, mat, par_r, par_p, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_p)

        grad_r = self.get(par_r, 'grad')
        grad_p = self.get(par_p, 'grad')
        div_mv = self.get(par_mv, 'div')
        grad_mv = self.get(par_mv, 'grad').transpose((0, 1, 3, 2)).copy()

        return (grad_r, grad_p, div_mv, grad_mv, mat, vg,
                get_default(term_mode, 1))

    def get_eval_shape(self, mat, par_r, par_p, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_p)
        return (n_el, 1, 1, 1), par_p.dtype
