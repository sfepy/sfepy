from sfepy.terms.terms import *
from sfepy.terms.termsNavierStokes import DivGradTerm

class AdjDivGradTerm(DivGradTerm):
    r"""
    :Description:
    Gateaux differential of :math:`\Psi(\ul{u}) = \int_{\Omega} \nu\
    \nabla \ul{v} : \nabla \ul{u}` w.r.t. :math:`\ul{u}` in the direction
    :math:`\ul{v}` or adjoint term to `dw_div_grad`.

    :Definition:
    .. math::
        w \delta_{u} \Psi(\ul{u}) \circ \ul{v}

    :Arguments:
        material_1 : :math:`w` (weight),
        material_2 : :math:`\nu` (viscosity),
        virtual    : :math:`\ul{v}`,
        state      : :math:`\ul{u}`
    """
    name = 'dw_adj_div_grad'
    arg_types = ('material_1', 'material_2', 'virtual', 'parameter')
    geometry = [(Volume, 'virtual')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        weight, viscosity, virtual, parameter = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        else:
            raise StopIteration

        vec = parameter()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, weight * viscosity,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

class AdjConvect1Term(Term):
    r"""
    :Description:
    The first adjoint term to nonlinear convective term `dw_convect`.

    :Definition:
    .. math::
        \int_{\Omega} ((\ul{v} \cdot \nabla) \ul{u}) \cdot \ul{w}

    :Arguments:
        virtual   : :math:`\ul{v}`,
        state     : :math:`\ul{w}`,
        parameter : :math:`\ul{u}`
    """
    name = 'dw_adj_convect1'
    arg_types = ('virtual', 'state', 'parameter' )
    geometry = [(Volume, 'virtual')]

    function = staticmethod(terms.dw_adj_convect1)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state, parameter = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec_w = state()
        vec_u = parameter()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_w, 0, vec_u, 0,
                                    bf, vg, ap.econn, chunk, mode )
            yield out, chunk, status

class AdjConvect2Term(AdjConvect1Term):
    r"""
    :Description:
    The second adjoint term to nonlinear convective term `dw_convect`.

    :Definition:
    .. math::
        \int_{\Omega} ((\ul{u} \cdot \nabla) \ul{v}) \cdot \ul{w}

    :Arguments:
        virtual   : :math:`\ul{v}`,
        state     : :math:`\ul{w}`,
        parameter : :math:`\ul{u}`
    """
    name = 'dw_adj_convect2'
    arg_types = ('virtual', 'state', 'parameter' )
    geometry = [(Volume, 'virtual')]

    function = staticmethod(terms.dw_adj_convect2)

class AdjSUPGCtabilizationTerm(Term):
    r"""
    :Description:
    Adjoint term to SUPG stabilization term `dw_st_supg_c`.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ [ ((\ul{v} \cdot \nabla)
        \ul{u}) ((\ul{u} \cdot \nabla) \ul{w}) + ((\ul{u} \cdot \nabla)
        \ul{u}) ((\ul{v} \cdot \nabla) \ul{w}) ]

    :Arguments:
        material  : :math:`\delta_K`,
        virtual   : :math:`\ul{v}`,
        state     : :math:`\ul{w}`,
        parameter : :math:`\ul{u}`
    """
    name = 'dw_st_adj_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    function = staticmethod(terms.dw_st_adj_supg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, virtual, par, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec_w = state()
        vec_u = par()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_w, 0,
                                    coef, bf, vg, ap.econn,
                                    chunk, mode )
            yield out, chunk, status

class SUPGPAdj1StabilizationTerm(Term):
    r"""
    :Description:
    The first adjoint term to SUPG stabilization term `dw_st_supg_p`.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ \nabla p (\ul{v} \cdot
        \nabla \ul{w})

    :Arguments:
        material  : :math:`\delta_K`,
        virtual   : :math:`\ul{v}`,
        state     : :math:`\ul{w}`,
        parameter : :math:`p`
    """
    name = 'dw_st_adj1_supg_p'
    arg_types = ('material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'parameter')]

    function = staticmethod(terms.dw_st_adj1_supg_p)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, virtual, state, par = self.get_args( **kwargs )
        cg = self.get_current_group()
        apr, vgr = virtual.get_approximation( cg, 'Volume' )
        app, vgp = par.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = apr.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim * n_ep, dim * n_ep )
            mode = 1
        else:
            raise StopIteration

        vec_w = state()
        vec_p = par()
        bf = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_p, 0, vec_w, 0,
                                    coef, bf, vgr, vgp,
                                    apr.econn, app.econn, chunk, mode )
            yield out, chunk, status

class SUPGPAdj2StabilizationTerm(Term):
    r"""
    :Description:
    The second adjoint term to SUPG stabilization term `dw_st_supg_p`
    as well as adjoint term to PSPG stabilization term `dw_st_pspg_c`.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ \nabla r (\ul{v} \cdot \nabla
        \ul{u})

    :Arguments:
        material  : :math:`\tau_K`,
        virtual   : :math:`\ul{v}`,
        parameter : :math:`\ul{u}`,
        state     : :math:`r`
    """
    name = 'dw_st_adj2_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    function = staticmethod(terms.dw_st_adj2_supg_p)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, virtual, par, state = self.get_args( **kwargs )
        cg = self.get_current_group()
        apr, vgr = virtual.get_approximation( cg, 'Volume' )
        apc, vgc = state.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_epr, 1 )
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape( self.integral_name )[3]
            shape = (chunk_size, 1, dim * n_epr, n_epc )
            mode = 1
        else:
            raise StopIteration

        vec_u = par()
        vec_r = state()
        bf = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_r, 0,
                                    coef, bf, vgr, vgc,
                                    apr.econn, apc.econn, chunk, mode )
            yield out, chunk, status

class TestPQTerm(Term):
    name = 'd_sd_test_pq'
    arg_types = ('parameter_p', 'parameter_q', 'parameter_mesh_velocity',
                'mode')
    geometry = [(Volume, 'parameter_p')]

    function = staticmethod(terms.d_sd_testPQ)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par_p, par_q, par_mv, mode = self.get_args( **kwargs )
        ap, vg = par_p.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_p = par_p()
        vec_q = par_q()
##         print par_p, par_q, par_mv, char_fun, mode
##         pause()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_p, 0, vec_q, 0, vec_mv, 0,
                                    bf, vg, ap.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDDivTerm(Term):
    r"""
    :Description:
    Sensitivity of Stokes term `dw_stokes` in 'div' mode.

    :Definition:
    .. math::
        \int_{\Omega_D} p [ (\nabla \cdot \ul{w}) (\nabla \cdot \ul{\Vcal})
        - \pdiff{\Vcal_k}{x_i} \pdiff{w_i}{x_k} ]

    :Arguments:
        parameter_u : :math:`\ul{u}`,
        parameter_p : :math:`p`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_div'
    arg_types = ('parameter_u', 'parameter_p', 'parameter_mesh_velocity',
                'mode')
    geometry = [(Volume, 'parameter_u'), (Volume, 'parameter_p'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_div)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par_u, par_p, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_u, vg_u = par_u.get_approximation( cg, 'Volume' )
        ap_p, vg_p = par_p.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_p.get_v_data_shape( self.integral_name )

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_p = par_p()
        vec_u = par_u()

        bf = ap_p.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_p, 0, vec_mv, 0,
                                    bf, vg_u, vg_p, vg_mv,
                                    ap_u.econn, ap_p.econn, ap_mv.econn,
                                    chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDDivGradTerm(Term):
    r"""
    :Description:
    Sensitivity of diffusion term `dw_div_grad`.

    :Definition:
    .. math::
        w \nu \int_{\Omega_D} [ \pdiff{u_i}{x_k} \pdiff{w_i}{x_k}
        (\nabla \cdot \ul{\Vcal})
        - \pdiff{\Vcal_j}{x_k} \pdiff{u_i}{x_j} \pdiff{w_i}{x_k}
        - \pdiff{u_i}{x_k} \pdiff{\Vcal_l}{x_k} \pdiff{w_i}{x_k} ]

    :Arguments:
        material_1  : :math:`w` (weight),
        material_2  : :math:`\nu` (viscosity),
        parameter_u : :math:`\ul{u}`,
        parameter_w : :math:`\ul{w}`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_div_grad'
    arg_types = ('material_1', 'material_2', 'parameter_u', 'parameter_w',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_div_grad)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        u_field: u_name ... fluid velocity
        w_field: w_name ... adjoint velocity
        mesh_v_field: mv_name ... mesh velsocity
        """
        weight, viscosity, par_u, par_w, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_u, vg_u = par_u.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape( self.integral_name )

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_w = par_w()
        vec_u = par_u()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_w, 0, vec_mv, 0,
                                    weight * viscosity[0,0,0,0], vg_u, vg_mv,
                                    ap_u.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDConvectTerm(Term):
    r"""
    :Description:
    Sensitivity of convective term `dw_convect`.

    :Definition:
    .. math::
        \int_{\Omega_D} [ u_k \pdiff{u_i}{x_k} w_i (\nabla \cdot \Vcal)
        - u_k \pdiff{\Vcal_j}{x_k} \pdiff{u_i}{x_j} w_i ]

    :Arguments:
        parameter_u : :math:`\ul{u}`,
        parameter_w : :math:`\ul{w}`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_convect'
    arg_types = ('parameter_u', 'parameter_w',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'), (Volume, 'parameter_w'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_convect)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par_u, par_w, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_u, vg_u = par_u.get_approximation( cg, 'Volume' )
        ap_w, vg_w = par_w.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape( self.integral_name )

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_w = par_w()
        vec_u = par_u()

        bf_u = ap_u.get_base( 'v', 0, self.integral_name )
        bf_w = ap_w.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_w, 0, vec_mv, 0,
                                    bf_u, bf_w, vg_u, vg_w, vg_mv,
                                    ap_u.econn, ap_w.econn, ap_mv.econn,
                                    chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class NSOFMinGrad1Term(Term):
    name = 'd_of_ns_min_grad1'
    arg_types = ('material_1', 'material_2', 'parameter')
    geometry = [(Volume, 'parameter')]

    function = staticmethod(terms.d_of_nsMinGrad)
    
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        parameter: fluid velocity
        """
        weight, viscosity, par = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        shape = (1,)

        vec = par()
        for out, chunk in self.char_fun( chunk_size, shape,
                                        zero = True, set_shape = False ):
            status = self.function( out, vec, 0, viscosity * weight,
                                    vg, ap.econn, chunk )
            yield out, chunk, status

class NSOFMinGrad2Term(NSOFMinGrad1Term):
    name = 'd_of_ns_min_grad2'
    arg_types = ('material_1', 'material_2', 'parameter')
    geometry = [(Volume, 'parameter')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        parameter: fluid velocity
        uses 1.0 instead of material.viscosity
        """
        weight, material, par = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        shape = (1,)

        vec = par()
        for out, chunk in self.char_fun( chunk_size, shape,
                                        zero = True, set_shape = False ):
            status = self.function( out, vec, 0, weight, vg, ap.econn, chunk )
            yield out, chunk, status

class NSOFSurfMinDPressTerm(Term):
    r"""
    :Description:
    Sensitivity of :math:`\Psi(p)`.

    :Definition:
    .. math::
        \delta \Psi(p) = \delta \left( \int_{\Gamma_{in}}p -
        \int_{\Gamma_{out}}bpress \right)

    :Arguments:
        material_1 : :math:`w` (weight),
        material_2 : :math:`bpress` (given pressure),
        parameter  : :math:`p`,
    """
    name = 'd_of_ns_surf_min_d_press'
    arg_types = ('material_1', 'material_2', 'parameter')
    geometry = [(Surface, 'parameter')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface',
                      function=terms.d_of_nsSurfMinDPress, **kwargs)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        Integrates over surface.
        
        material_1: weight
        material_2: given constant pressure on outlet
        parameter: fluid pressure
        """
        weight, press, par = self.get_args( **kwargs )
        ap, sg = par.get_approximation( self.get_current_group(), 'Surface' )
        shape = (1,)

        sd = ap.surface_data[self.region.name]

        vec = par()
##         sg.str( sys.stdout, 0 )
##         print self.char_fun.region
##         pause()
        bf = ap.get_base( sd.face_type, 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape,
                                        zero = True, set_shape = False ):
            lchunk = self.char_fun.get_local_chunk()
            status = self.function( out, vec, 0, weight, press,
                                    bf, sg, sd.econn, lchunk, 0 )
            yield out, lchunk, status

class NSOFSurfMinDPressDiffTerm(NSOFSurfMinDPressTerm):
    r"""
    :Description:
    Gateaux differential of :math:`\Psi(p)` w.r.t. :math:`p` in the
    direction :math:`q`.

    :Definition:
    .. math::
        w \delta_{p} \Psi(p) \circ q

    :Arguments:
        material : :math:`w` (weight),
        virtual  : :math:`q`,
    """
    name = 'dw_of_ns_surf_min_d_press_diff'
    arg_types = ('material', 'virtual')
    geometry = [(Surface, 'virtual')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        Integrates over surface.
        
        material_1: weight
        virtual: fluid pressure-like
        """
        weight, virtual = self.get_args( **kwargs )
        ap, sg = virtual.get_approximation( self.get_current_group(), 'Surface' )
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape( self.integral_name,
                                                     self.region.name )
        shape = (chunk_size, 1, n_fp, 1 )

        sd = ap.surface_data[self.region.name]

        aux = nm.array( [], dtype = nm.float64 )
##         print sd.econn, sd.econn.shape
        bf = ap.get_base( sd.face_type, 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
##             print chunk, chunk.shape
##             print lchunk, lchunk.shape
##             pause()
            status = self.function( out, aux, 0, weight, 0.0,
                                    bf, sg, sd.econn, lchunk, 1 )
            yield out, lchunk, status

class SDGradDivStabilizationTerm(Term):
    r"""
    :Description:
    Sensitivity of stabilization term `dw_st_grad_div`.

    :Definition:
    .. math::
        \gamma \int_{\Omega_D} [ (\nabla \cdot \ul{u}) (\nabla \cdot \ul{w})
        (\nabla \cdot \ul{\Vcal})
        - \pdiff{u_i}{x_k} \pdiff{\Vcal_k}{x_i} (\nabla \cdot \ul{w})
        - (\nabla \cdot \ul{u}) \pdiff{w_i}{x_k} \pdiff{\Vcal_k}{x_i} ]

    :Arguments:
        material    : :math:`\gamma`,
        parameter_w : :math:`\ul{w}`,
        parameter_u : :math:`\ul{u}`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_grad_div'
    arg_types = ('material', 'parameter_w', 'parameter_u',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_st_grad_div)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_w, par_u, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_u, vg_u = par_u.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape( self.integral_name )

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
    :Description:
    Sensitivity of stabilization term `dw_st_supg_c`.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ [ (\ul{u} \cdot \nabla)u_k
        (\ul{u} \cdot \nabla)w_k (\nabla \cdot \Vcal) -
        (\ul{u} \cdot \nabla)\Vcal_i \pdiff{u_k}{x_i}
        (\ul{u} \cdot \nabla)w_k - (\ul{u} \cdot \nabla)u_k
        (\ul{u} \cdot \nabla)\Vcal_i \pdiff{w_k}{x_i} ]

    :Arguments:
        material    : :math:`\delta_K`,
        parameter_w : :math:`\ul{w}`,
        parameter_b : :math:`\ul{b}`,
        parameter_u : :math:`\ul{u}`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_supg_c'
    arg_types = ('material', 'parameter_w', 'parameter_b', 'parameter_u',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_st_supg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_w, par_b, par_u, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_u, vg_u = par_u.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape( self.integral_name )

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_w = par_w()
        vec_b = par_b()
        vec_u = par_u()

        bf = ap_u.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_b, 0,
                                    vec_w, 0, vec_mv, 0,
                                    bf, coef, vg_u, vg_mv,
                                    ap_u.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDPSPGCStabilizationTerm(Term):
    r"""
    :Description:
    Sensitivity of stabilization terms `dw_st_supg_p` or `dw_st_pspg_c`.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \delta_K\
        [ \pdiff{p}{x_i} (\ul{u} \cdot \nabla)w_i (\nabla \cdot \Vcal) -
        \pdiff{p}{x_k} \pdiff{\Vcal_k}{x_i} (\ul{u} \cdot \nabla)w_i
        - \pdiff{p}{x_k} (\ul{u} \cdot \nabla)\Vcal_k  \pdiff{w_i}{x_k} ]

    :Arguments:
        material    : :math:`\delta_K`,
        parameter_r : :math:`r`,
        parameter_b : :math:`\ul{b}`,
        parameter_u : :math:`\ul{u}`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_pspg_c'
    arg_types = ('material', 'parameter_r', 'parameter_b', 'parameter_u',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_r'),
                (Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_st_pspg_c)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_r, par_b, par_u, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_u, vg_u = par_u.get_approximation( cg, 'Volume' )
        ap_r, vg_r = par_r.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_u.get_v_data_shape( self.integral_name )

        if not chunk_size:
            chunk_size = n_el
        shape = (chunk_size, 1, 1, 1)

        vec_mv = par_mv()
        vec_r = par_r()
        vec_b = par_b()
        vec_u = par_u()

        bf = ap_u.get_base( 'v', 0, self.integral_name )
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
    :Description:
    Sensitivity of stabilization term `dw_st_pspg_p`.

    :Definition:
    .. math::
        \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ [ (\nabla p \cdot \nabla q)
        (\nabla \cdot \Vcal)
        - \pdiff{p}{x_k} (\nabla \Vcal_k \cdot \nabla q)
        - (\nabla p \cdot \nabla \Vcal_k) \pdiff{q}{x_k} ]

    :Arguments:
        material    : :math:`\tau_K`,
        parameter_r : :math:`r`,
        parameter_p : :math:`p`,
        parameter_mesh_velocity : :math:`\ul{\Vcal}`,
        mode        : 1 (sensitivity) or 0 (original term value)
    """
    name = 'd_sd_st_pspg_p'
    arg_types = ('material', 'parameter_r', 'parameter_p',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_p'),
                (Volume, 'parameter_mesh_velocity')]

    function = staticmethod(terms.d_sd_st_pspg_p)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        coef, par_r, par_p, par_mv, mode = self.get_args( **kwargs )
        cg = self.get_current_group()
        ap_p, vg_p = par_p.get_approximation( cg, 'Volume' )
        ap_mv, vg_mv = par_mv.get_approximation( cg, 'Volume' )
        n_el, n_qp, dim, n_ep = ap_p.get_v_data_shape( self.integral_name )
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
