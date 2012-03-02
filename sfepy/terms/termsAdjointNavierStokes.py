import numpy as nm

from sfepy.terms.terms import Term, terms

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

    function = staticmethod(terms.term_ns_asm_div_grad)

    def get_fargs(self, mat1, mat2, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if diff_var is None:
            grad = grad_as_vector(self.get(state, 'grad'))
            fmode = 0

        else:
            grad = nm.array([0], ndmin=4, dtype=nm.float64)
            fmode = 1

        return grad, mat1 * mat2, vg, fmode

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

    function = staticmethod(terms.dw_adj_convect1)

    def get_fargs(self, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        val_w = self.get(state, 'val')
        grad_u = self.get(parameter, 'grad') # No transposition here!

        fmode = diff_var is not None

        return val_w, grad_u, vg.bf, vg, fmode

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

    function = staticmethod(terms.dw_adj_convect2)

    def get_fargs(self, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        val_w = self.get(state, 'val')
        val_u = self.get(parameter, 'val')

        fmode = diff_var is not None

        return val_w, val_u, vg.bf, vg, fmode

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

    function = staticmethod(terms.dw_st_adj_supg_c)

    def get_fargs(self, mat, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, vg = self.get_approximation(state)

        val_u = self.get(parameter, 'val')
        grad_u = self.get(parameter, 'grad').transpose((0, 1, 3, 2)).copy()
        conn = ap.get_connectivity(self.region, self.integration)

        fmode = diff_var is not None

        return state(), val_u, grad_u, mat, vg.bf, vg, conn, fmode

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

    function = staticmethod(terms.dw_st_adj1_supg_p)

    def get_fargs(self, mat, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap_w, vg_w = self.get_approximation(state)

        grad_p = self.get(parameter, 'grad')
        conn_w = ap_w.get_connectivity(self.region, self.integration)

        fmode = diff_var is not None

        return state(), grad_p, mat, vg_w.bf, vg_w, conn_w, fmode

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

    function = staticmethod(terms.dw_st_adj2_supg_p)

    def get_fargs(self, mat, virtual, parameter, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap_r, vg_r = self.get_approximation(state)
        vg_u, _ = self.get_mapping(parameter)

        grad_u = self.get(parameter, 'grad').transpose((0, 1, 3, 2)).copy()
        conn_r = ap_r.get_connectivity(self.region, self.integration)

        fmode = diff_var is not None

        return grad_u, state(), mat, vg_u.bf, vg_u, vg_r, conn_r, fmode

class TestPQTerm(Term):
    name = 'd_sd_test_pq'
    arg_types = ('parameter_p', 'parameter_q', 'parameter_mesh_velocity',
                'mode')

    function = staticmethod(terms.d_sd_testPQ)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par_p, par_q, par_mv, mode = self.get_args( **kwargs )
        ap, vg = self.get_approximation(par_p)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_p = par_p()
        vec_q = par_q()
##         print par_p, par_q, par_mv, char_fun, mode
##         pause()
        bf = ap.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_p, 0, vec_q, 0, vec_mv, 0,
                                    bf, vg, ap.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDDivTerm(Term):
    r"""
    Sensitivity (shape derivative) of Stokes term `dw_stokes` in 'div' mode.

    Supports the following term modes: 1 (sensitivity) or 0 (original term
    value).

    :Definition:

    .. math::
        \int_{\Omega_D} p [ (\nabla \cdot \ul{w}) (\nabla \cdot \ul{\Vcal})
        - \pdiff{\Vcal_k}{x_i} \pdiff{w_i}{x_k} ]

    :Arguments:
        - parameter_u : :math:`\ul{u}`
        - parameter_p : :math:`p`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
    """
    name = 'd_sd_div'
    arg_types = ('parameter_u', 'parameter_p', 'parameter_mesh_velocity')

    function = staticmethod(terms.d_sd_div)

    def get_fargs(self, par_u, par_p, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        div_u = self.get(par_u, 'div')
        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        val_p = self.get(par_p, 'val')
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        return div_u, grad_u, val_p, div_mv, grad_mv, vg, term_mode

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
        w \nu \int_{\Omega_D} [ \pdiff{u_i}{x_k} \pdiff{w_i}{x_k}
        (\nabla \cdot \ul{\Vcal})
        - \pdiff{\Vcal_j}{x_k} \pdiff{u_i}{x_j} \pdiff{w_i}{x_k}
        - \pdiff{u_i}{x_k} \pdiff{\Vcal_l}{x_k} \pdiff{w_i}{x_k} ]

    :Arguments:
        - material_1  : :math:`w` (weight)
        - material_2  : :math:`\nu` (viscosity)
        - parameter_u : :math:`\ul{u}`
        - parameter_w : :math:`\ul{w}`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
    """
    name = 'd_sd_div_grad'
    arg_types = ('material_1', 'material_2', 'parameter_u', 'parameter_w',
                 'parameter_mesh_velocity')

    function = staticmethod(terms.d_sd_div_grad)

    def get_fargs(self, mat1, mat2, par_u, par_w, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        grad_w = grad_as_vector(self.get(par_w, 'grad'))
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        return grad_u, grad_w, div_mv, grad_mv, mat1 * mat2, vg, term_mode

    def get_eval_shape(self, mat1, mat2, par_u, par_w, par_mv,
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
        \int_{\Omega_D} [ u_k \pdiff{u_i}{x_k} w_i (\nabla \cdot \Vcal)
        - u_k \pdiff{\Vcal_j}{x_k} \pdiff{u_i}{x_j} w_i ]

    :Arguments:
        - parameter_u : :math:`\ul{u}`
        - parameter_w : :math:`\ul{w}`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
    """
    name = 'd_sd_convect'
    arg_types = ('parameter_u', 'parameter_w', 'parameter_mesh_velocity')

    function = staticmethod(terms.d_sd_convect)

    def get_fargs(self, par_u, par_w, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(par_u)

        val_u = self.get(par_u, 'val')
        grad_u = grad_as_vector(self.get(par_u, 'grad'))
        val_w = self.get(par_w, 'val')
        div_mv = self.get(par_mv, 'div')
        grad_mv = grad_as_vector(self.get(par_mv, 'grad'))

        return val_u, grad_u, val_w, div_mv, grad_mv, vg, term_mode

    def get_eval_shape(self, par_u, par_w, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)
        return (n_el, 1, 1, 1), par_u.dtype

class NSOFMinGradTerm(Term):
    name = 'd_of_ns_min_grad'
    arg_types = ('material_1', 'material_2', 'parameter')

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
    name = 'd_of_ns_surf_min_d_press'
    arg_types = ('material_1', 'material_2', 'parameter')
    integration = 'surface'

    function = staticmethod(terms.d_of_nsSurfMinDPress)

    def get_fargs(self, weight, bpress, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        val_p = self.get(parameter, 'val')

        return val_p, weight, bpress, sg.bf, sg, 0

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

    def get_fargs(self, weight, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        aux = nm.array([0], ndmin=4, dtype=nm.float64)

        return aux, weight, 0.0, sg.bf, sg, 1

class SDGradDivStabilizationTerm(Term):
    r"""
    Sensitivity of stabilization term `dw_st_grad_div`.

    :Definition:

    .. math::
        \gamma \int_{\Omega_D} [ (\nabla \cdot \ul{u}) (\nabla \cdot \ul{w})
        (\nabla \cdot \ul{\Vcal})
        - \pdiff{u_i}{x_k} \pdiff{\Vcal_k}{x_i} (\nabla \cdot \ul{w})
        - (\nabla \cdot \ul{u}) \pdiff{w_i}{x_k} \pdiff{\Vcal_k}{x_i} ]

    :Arguments:
        - material    : :math:`\gamma`
        - parameter_w : :math:`\ul{w}`
        - parameter_u : :math:`\ul{u}`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_grad_div'
    arg_types = ('material', 'parameter_w', 'parameter_u',
                'parameter_mesh_velocity', 'mode')

    function = staticmethod(terms.d_sd_st_grad_div)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_w, par_u, par_mv, mode = self.get_args( **kwargs )
        ap_u, vg_u = self.get_approximation(par_u)
        ap_mv, vg_mv = self.get_approximation(par_mv)
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape(self.integral)

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_w = par_w()
        vec_u = par_u()

        try:
            coef = float( coef )
        except TypeError:
            coef = 1.0

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_w, 0,
                                    vec_mv, 0, coef, vg_u, vg_mv,
                                    ap_u.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDSUPGCStabilizationTerm(Term):
    r"""
    Sensitivity of stabilization term `dw_st_supg_c`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ [ (\ul{u} \cdot \nabla)u_k
        (\ul{u} \cdot \nabla)w_k (\nabla \cdot \Vcal) -
        (\ul{u} \cdot \nabla)\Vcal_i \pdiff{u_k}{x_i}
        (\ul{u} \cdot \nabla)w_k - (\ul{u} \cdot \nabla)u_k
        (\ul{u} \cdot \nabla)\Vcal_i \pdiff{w_k}{x_i} ]

    :Arguments:
        - material    : :math:`\delta_K`
        - parameter_w : :math:`\ul{w}`
        - parameter_b : :math:`\ul{b}`
        - parameter_u : :math:`\ul{u}`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_supg_c'
    arg_types = ('material', 'parameter_w', 'parameter_b', 'parameter_u',
                'parameter_mesh_velocity', 'mode')

    function = staticmethod(terms.d_sd_st_supg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_w, par_b, par_u, par_mv, mode = self.get_args( **kwargs )
        ap_u, vg_u = self.get_approximation(par_u)
        ap_mv, vg_mv = self.get_approximation(par_mv)
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape(self.integral)

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_w = par_w()
        vec_b = par_b()
        vec_u = par_u()

        bf = ap_u.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_b, 0,
                                    vec_w, 0, vec_mv, 0,
                                    bf, coef, vg_u, vg_mv,
                                    ap_u.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDPSPGCStabilizationTerm(Term):
    r"""
    Sensitivity of stabilization terms `dw_st_supg_p` or `dw_st_pspg_c`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\
        [ \pdiff{p}{x_i} (\ul{u} \cdot \nabla)w_i (\nabla \cdot \Vcal) -
        \pdiff{p}{x_k} \pdiff{\Vcal_k}{x_i} (\ul{u} \cdot \nabla)w_i
        - \pdiff{p}{x_k} (\ul{u} \cdot \nabla)\Vcal_k  \pdiff{w_i}{x_k} ]

    :Arguments:
        - material    : :math:`\delta_K`
        - parameter_r : :math:`r`
        - parameter_b : :math:`\ul{b}`
        - parameter_u : :math:`\ul{u}`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_pspg_c'
    arg_types = ('material', 'parameter_r', 'parameter_b', 'parameter_u',
                'parameter_mesh_velocity', 'mode')

    function = staticmethod(terms.d_sd_st_pspg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_r, par_b, par_u, par_mv, mode = self.get_args( **kwargs )
        ap_u, vg_u = self.get_approximation(par_u)
        ap_r, vg_r = self.get_approximation(par_r)
        ap_mv, vg_mv = self.get_approximation(par_mv)
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape(self.integral)

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_r = par_r()
        vec_b = par_b()
        vec_u = par_u()

        bf = ap_u.get_base('v', 0, self.integral)
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_b, 0,
                                    vec_r, 0, vec_mv, 0,
                                    bf, coef, vg_u, vg_r, vg_mv,
                                    ap_u.econn, ap_r.econn, ap_mv.econn,
                                    chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDPSPGPStabilizationTerm(Term):
    r"""
    Sensitivity of stabilization term `dw_st_pspg_p`.

    :Definition:

    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ [ (\nabla p \cdot \nabla q)
        (\nabla \cdot \Vcal)
        - \pdiff{p}{x_k} (\nabla \Vcal_k \cdot \nabla q)
        - (\nabla p \cdot \nabla \Vcal_k) \pdiff{q}{x_k} ]

    :Arguments:
        - material    : :math:`\tau_K`
        - parameter_r : :math:`r`
        - parameter_p : :math:`p`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
        - mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_pspg_p'
    arg_types = ('material', 'parameter_r', 'parameter_p',
                'parameter_mesh_velocity', 'mode')

    function = staticmethod(terms.d_sd_st_pspg_p)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_r, par_p, par_mv, mode = self.get_args( **kwargs )
        ap_p, vg_p = self.get_approximation(par_p)
        ap_mv, vg_mv = self.get_approximation(par_mv)
        n_el, n_qp, dim, n_ep = ap_p.get_v_data_shape(self.integral)
        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_r = par_r()
        vec_p = par_p()

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_p, 0, vec_r, 0,
                                    vec_mv, 0, coef, vg_p, vg_mv,
                                    ap_p.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status
