from sfepy.terms.terms import *
from sfepy.terms.termsNavierStokes import DivGradTerm
from sfepy.terms.utils import fix_scalar_in_el

class AdjDivGrad1Term( DivGradTerm ):
    """Uses material.viscosity as viscosity."""
    name = 'dw_adj_div_grad1'
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
            status = self.function( out, vec, 0,
                                    weight * viscosity,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

class AdjDivGrad2Term( DivGradTerm ):
    """Uses 1.0 as viscosity."""
    name = 'dw_adj_div_grad2'
    arg_types = ('material_1', 'material_2', 'virtual', 'parameter')
    geometry = [(Volume, 'virtual')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        weight, material, virtual, parameter = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if diff_var is None:
            shape = (chunk_size, 1, dim * n_ep, 1 )
            mode = 0
        else:
            raise StopIteration

        vec = parameter()
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, weight,
                                    vg, ap.econn, chunk, mode )
            yield out, chunk, status

class AdjConvect1Term( Term ):
    name = 'dw_adj_convect1'
    arg_types = ('virtual', 'state', 'parameter' )
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_adj_convect1 )
        
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

class AdjConvect2Term( AdjConvect1Term ):
    name = 'dw_adj_convect2'
    arg_types = ('virtual', 'state', 'parameter' )
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_adj_convect2 )

class AdjSUPGCtabilizationTerm( Term ):
    name = 'dw_st_adj_supg_c'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_adj_supg_c )
        
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

        coef_in_el = fix_scalar_in_el( coef, n_el, nm.float64 )

        vec_w = state()
        vec_u = par()
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_w, 0,
                                    coef_in_el, bf, vg, ap.econn,
                                    chunk, mode )
            yield out, chunk, status

class SUPGPAdj1StabilizationTerm( Term ):
    name = 'dw_st_adj1_supg_p'
    arg_types = ('material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual'), (Volume, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_adj1_supg_p )
        
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

        coef_in_el = fix_scalar_in_el( coef, n_el, nm.float64 )
            
        vec_w = state()
        vec_p = par()
        bf = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_p, 0, vec_w, 0,
                                    coef_in_el, bf, vgr, vgp,
                                    apr.econn, app.econn, chunk, mode )
            yield out, chunk, status

class SUPGPAdj2StabilizationTerm( Term ):
    name = 'dw_st_adj2_supg_p'
    arg_types = ('material', 'virtual', 'parameter', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_st_adj2_supg_p )
        
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

        coef_in_el = fix_scalar_in_el( coef, n_el, nm.float64 )
            
        vec_u = par()
        vec_r = state()
        bf = apr.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_r, 0,
                                    coef_in_el, bf, vgr, vgc,
                                    apr.econn, apc.econn, chunk, mode )
            yield out, chunk, status

class TestPQTerm( Term ):
    name = 'd_sd_test_pq'
    arg_types = ('parameter_p', 'parameter_q', 'parameter_mesh_velocity',
                'mode')
    geometry = [(Volume, 'parameter_p')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_test_pq )

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

class SDDivTerm( Term ):
    name = 'd_sd_div'
    arg_types = ('parameter_u', 'parameter_p', 'parameter_mesh_velocity',
                'mode')
    geometry = [(Volume, 'parameter_u'), (Volume, 'parameter_p'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_div )

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

class SDDivGradTerm( Term ):
    name = 'd_sd_div_grad'
    arg_types = ('material_1', 'material_2', 'parameter_u', 'parameter_w',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_div_grad )

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
                                    weight * viscosity, vg_u, vg_mv,
                                    ap_u.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDConvectTerm( Term ):
    name = 'd_sd_convect'
    arg_types = ('parameter_u', 'parameter_w',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'), (Volume, 'parameter_w'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_convect )

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

class NSOFMinGrad1Term( Term ):
    name = 'd_of_ns_min_grad1'
    arg_types = ('material_1', 'material_2', 'parameter')
    geometry = [(Volume, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_of_nsMinGrad )
    
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

class NSOFMinGrad2Term( NSOFMinGrad1Term ):
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

class NSOFSurfMinDPressTerm( Term ):
    name = 'd_of_ns_surf_min_d_press'
    arg_types = ('material_1', 'material_2', 'parameter')
    geometry = [(Surface, 'parameter')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign,
                       terms.d_of_nsSurfMinDPress )
        self.dof_conn_type = 'surface'

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

class NSOFSurfMinDPressDiffTerm( NSOFSurfMinDPressTerm ):
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

        shape = (chunk_size, 1, sg.n_fp, 1 )

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

class SDGradDivStabilizationTerm( Term ):
    name = 'd_sd_st_grad_div'
    arg_types = ('material', 'parameter_w', 'parameter_u',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_st_grad_div )
        
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

class SDSUPGCStabilizationTerm( Term ):
    name = 'd_sd_st_supg_c'
    arg_types = ('material', 'parameter_w', 'parameter_b', 'parameter_u',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_st_supg_c )
        
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

        coef_in_el = fix_scalar_in_el( coef, n_el, nm.float64 )
        bf = ap_u.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_b, 0,
                                    vec_w, 0, vec_mv, 0,
                                    bf, coef_in_el, vg_u, vg_mv,
                                    ap_u.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDPSPGCStabilizationTerm( Term ):
    name = 'd_sd_st_pspg_c'
    arg_types = ('material', 'parameter_r', 'parameter_b', 'parameter_u',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_r'),
                (Volume, 'parameter_u'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_st_pspg_c)
        
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

        coef_in_el = fix_scalar_in_el( coef, n_el, nm.float64 )
        bf = ap_u.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_u, 0, vec_b, 0,
                                    vec_r, 0, vec_mv, 0,
                                    bf, coef_in_el, vg_u, vg_r, vg_mv,
                                    ap_u.econn, ap_r.econn, ap_mv.econn,
                                    chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status

class SDPSPGPStabilizationTerm( Term ):
    name = 'd_sd_st_pspg_p'
    arg_types = ('material', 'parameter_r', 'parameter_p',
                'parameter_mesh_velocity', 'mode')
    geometry = [(Volume, 'parameter_p'),
                (Volume, 'parameter_mesh_velocity')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.d_sd_st_pspg_p)
        
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

        coef_in_el = fix_scalar_in_el( coef, n_el, nm.float64 )
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec_p, 0, vec_r, 0,
                                    vec_mv, 0, coef_in_el, vg_p, vg_mv,
                                    ap_p.econn, ap_mv.econn, chunk, mode )
            out1 = nm.sum( nm.squeeze( out ) )
            yield out1, chunk, status
