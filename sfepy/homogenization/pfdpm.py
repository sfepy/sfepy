from sfepy.base.base import *
from sfepy.fem.evaluate import eval_term_op
from sfepy.fem.meshio import HDF5MeshIO
from sfepy.solvers.ts import TimeStepper
from sfepy.base.la import MatrixAction, eig
from sfepy.base.plotutils import spy
from sfepy.homogenization.utils import build_op_pi, iter_sym
from sfepy.homogenization.coefficients import Coefficients
from sfepy.fem.history import History
from sfepy.solvers import Solver

##
# 02.05.2007, c
def debug_data( name, data_template, data, eps = 1e-10 ):

    diff = nm.abs( data_template - data )
    iis = nm.where( diff >= eps )
    ok = not len( iis[0] )

    if ok:
        output( 'data: %s ... ok' % name )
    else:
        output( 'data: %s ... failed' % name )
        output( 'data template:' )
        print data_template
        output( 'data:' )
        print data
        output( 'indices of failed entries:' )
        print iis
        output( 'absolute diff:' )
        print diff
        output( 'relative (w.r.t. template) diff:' )
        print diff / nm.abs( data_template )
        pause()

##
# c: 16.07.2007, r: 27.02.2008
def volumes( problem, variables_all, coef_eqs, kwargs ):
    region_names = kwargs['regions']
    coef_term = coef_eqs[kwargs['term']]

    variables = select_by_names( variables_all, kwargs['variables'] )
    problem.set_variables( variables )

    volumes = {}
    for region_name in region_names:
        val = eval_term_op( None, coef_term % region_name, problem )
        volumes[region_name] = nm.asarray( val, dtype = nm.float64 )

    return volumes, region_names[-1]
                  
##
# c: 09.03.2007, r: 02.04.2008
def coef_e( problem, corrs_rs, pis, variables_all, coef_eqs, kwargs ):
    """Elastic coefficient."""
    coef_term = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2
    coef = nm.zeros( (sym, sym), dtype = nm.float64 )

    indx = corrs_rs.di.indx['uc']
    for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
        omega1 = corrs_rs.states_rs[irr,icr][indx]
        pi1 = pis[irr,icr] + omega1
        problem.variables[var_names[0]].data_from_data( pi1 )
            
        for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
            omega2 = corrs_rs.states_rs[irc,icc][indx]
            pi2 = pis[irc,icc] + omega2
            problem.variables[var_names[1]].data_from_data( pi2 )

            val = eval_term_op( None, coef_term, problem, call_mode = 'd_eval' )

            coef[ir,ic] = val
    return coef

##
# c: 05.03.2008, r: 05.03.2008
def coef_biot_e( problem, corrs_rs, variables_all, coef_eqs, kwargs ):
    """Elastic Biot coefficient."""
    coef_term_a, coef_term_b = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2

    coef = eval_term_op( None, coef_term_a, problem, shape = (sym,),
                         mode = 'const' )
    print coef, coef.shape

    one_var = problem.variables[var_names[0]]
    one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
    one_var.data_from_data( one )

    indx = corrs_rs.di.indx['uc']
    for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
        omega = corrs_rs.states_rs[ir,ic][indx]
        problem.variables[var_names[1]].data_from_data( omega )
        val = eval_term_op( None, coef_term_b, problem, call_mode = 'd_eval' )
        print ii, ir, ic, val
        coef[ii] += val
    return coef

##
# c: 07.03.2008, r: 23.06.2008
def coef_biot_h( problem, filenames_rs, variables_all, coef_eqs, kwargs ):
    """Fading memory Biot coefficient."""
    coef_term = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )
##     problem.set_equations( {'eq': 'dw_biot_grad.i1.Y( m.alpha, vc, p_y_one )'} )
##     problem.time_update( conf_ebc = {}, conf_epbc = {} )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2

    ts = TimeStepper( *HDF5MeshIO( filenames_rs[0,0] ).read_time_stepper() )
    coef = nm.zeros( (ts.n_step, sym), dtype = nm.float64 )

    vv = problem.variables
    one = nm.ones( (vv[var_names[0]].field.n_nod,), dtype = nm.float64 )
    vv[var_names[0]].data_from_data( one )

    one_m = one[:vv[var_names[1]].field.n_nod]
    vv[var_names[1]].data_from_data( one_m )

    for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
        io = HDF5MeshIO( filenames_rs[ir,ic] )
        for step, time in ts:
            data = io.read_data( step )
            vv[var_names[2]].data_from_data( data['u'].data )
            vv[var_names[3]].data_from_data( data['dp'].data )
            val = eval_term_op( None, coef_term, problem, call_mode = 'd_eval' )

##             vv = eval_term_op( None, 'd_biot_div.i1.Y( m.alpha, p_y_one, Pi1 )',
##                              problem )
##             print vv
##             v = eval_term_op( None, 'dw_biot_grad.i1.Y( m.alpha, vc, p_y_one )',
##                             problem )
##             print nm.sum( data['u'].data * v ) - vv
##             debug()

##             print '>>>>>>>', step, ii, ir, ic, val
            coef[step,ii] = val
    return coef

##
# c: 16.05.2008, r: 23.06.2008
def coef_biot_h2( problem, filenames_rs, filename_pressure, corrss_rs, corrs_pressure,
                variables_all, coef_eqs, kwargs ):
    """Fading memory Biot coefficient, alternative form."""
    coef_term = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2

    ts = TimeStepper( *HDF5MeshIO( filenames_rs[0,0] ).read_time_stepper() )
    coef = nm.zeros( (ts.n_step, sym), dtype = nm.float64 )

    vv = problem.variables

##     one = nm.ones( (vv['p_y_one'].field.n_nod,), dtype = nm.float64 )
##     vv['p_y_one'].data_from_data( one )

##     one_m = one[:vv['p_ym_one'].field.n_nod]
##     vv['p_ym_one'].data_from_data( one_m )

    p0 = HDF5MeshIO( filename_pressure ).read_data( 0 )['p'].data
    vv[var_names[1]].data_from_data( p0 )
##     indx_u = corrs_pressure.di.indx['uc']
##     indx_p = corrs_pressure.di.indx['pc']
    for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
        io = HDF5MeshIO( filenames_rs[ir,ic] )
        for step, time in ts:
            data_rs = io.read_data( step )
            vv[var_names[0]].data_from_data( data_rs['p'].data )
##             vv['pp3'].data_from_data( data_rs['dp'].data )
##             vv['Pi1'].data_from_data( data_rs['u'].data )
            val = eval_term_op( None, coef_term, problem, call_mode = 'd_eval' )

##             val1 = eval_term_op( None, 'd_biot_div.i1.Y( m.alpha, p_y_one, Pi1 )',
##                              problem )
##             print val1
##             omega_p = corrs_pressure.state[indx_u]
##             vv['Pi2'].data_from_data( omega_p )
            
##             pp = corrs_pressure.state[indx_p]
##             vv['pp2'].data_from_data( pp )
##             val2 = eval_term_op( None, """d_lin_elastic.i2.Y( m.D, Pi2, Pi1 )
##             - d_biot_div.i1.Ym( m.alpha, pp2, Pi1 )""",
##                              problem )
##             print val2

##             val3 = eval_term_op( None,
##                                """d_volume_wdot.i1.Ym( m.imu, p_ym_one, pp3 )""",
##                                problem )
##             print val3
##             val4 = eval_term_op( None,
##                                """- d_volume_wdot.i1.Ym( m.imu, pp2, pp3 )
##                                - d_biot_div.i1.Ym( m.alpha, pp3, Pi2 )""",
##                                problem )
##             print val4

##             val5 = eval_term_op( None,
##                                """d_biot_div.i1.Ym( m.alpha, pp3, Pi2 )""",
##                                problem )
##             print val5
## ##             val5 = eval_term_op( None,
## ##                                """d_biot_grad.i1.Ym( m.alpha, Pi2, pp1 )""",
## ##                                problem )
## ##             print val5
##             val6 = eval_term_op( None,
##                                """d_lin_elastic.i2.Y( m.D, Pi2, Pi1 )""",
##                                problem )
##             print val6

##             val7 = eval_term_op( None, """d_volume_wdot.i1.Ym( m.imu, pp2, pp3 )""",
##                              problem )
##             print val7

##             val8 = eval_term_op( None,
##                                """- d_diffusion.i1.Ym( m.K, pp1, pp2 )
##                                - d_biot_div.i1.Ym( m.alpha, pp2, Pi1 )""",
##                                problem )
##             print val8
##             val9 = eval_term_op( None,
##                                """d_diffusion.i1.Ym( m.K, pp1, pp2 )""",
##                                problem )
##             print val9
##             print val2 + val4
##             print val2 - val6 - val8
##             print val1 + val3
##             print val

##             debug()

##             print '>>>>>>>', step, ii, ir, ic, val
            coef[step,ii] = val
    return coef

##
# c: 05.03.2008, r: 17.03.2008
def coef_biot_mir( problem, corrs_pressure, variables_all, coef_eqs, kwargs ):
    """Instantaneous reciprocal Biot modulus."""
    coef_term_a, coef_term_b = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    coef = eval_term_op( None, coef_term_a, problem, shape = (1,),
                         mode = 'const' )
    print coef, coef.shape

    one_var = problem.variables[var_names[0]]
    one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
    one_var.data_from_data( one )

    one_var = problem.variables[var_names[1]]
    one_m = one[:one_var.field.n_nod]
    one_var.data_from_data( one_m )

    omega = corrs_pressure.state[corrs_pressure.di.indx['uc']]
    problem.variables[var_names[2]].data_from_data( omega )

    pp = corrs_pressure.state[corrs_pressure.di.indx['pc']]
    problem.variables[var_names[3]].data_from_data( pp )

    val = eval_term_op( None, coef_term_b, problem, call_mode = 'd_eval' )
    print val
    coef += val

    return coef

##
# c: 07.03.2008, r: 23.06.2008
def coef_biot_mhr( problem, filename_pressure, variables_all, coef_eqs, kwargs ):
    """Fading memory reciprocal Biot modulus."""
    coef_term = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    io = HDF5MeshIO( filename_pressure )
    ts = TimeStepper( *io.read_time_stepper() )
    coef = nm.zeros( (ts.n_step, 1), dtype = nm.float64 )

    vv = problem.variables

    one_var = vv[var_names[0]]
    one = nm.ones( (one_var.field.n_nod,), dtype = nm.float64 )
    one_var.data_from_data( one )

    one_var = vv[var_names[1]]
    one_m = one[:one_var.field.n_nod]
    one_var.data_from_data( one_m )

    for step, time in ts:
        data = io.read_data( step )
        vv[var_names[2]].data_from_data( data['u'].data )
        vv[var_names[3]].data_from_data( data['dp'].data )
        val = eval_term_op( None, coef_term, problem, call_mode = 'd_eval' )
        coef[step] = val
    return coef

##
# c: 13.03.2007, r: 23.06.2008
def coef_b_bar( problem, states_rs, filenames_rs, coef_term = None ):
    """Instantaneous viscosity coefficient."""
    dim = problem.domain.mesh.dim
    coef = nm.zeros( (dim**2, dim**2), dtype = nm.float64 )
    
    if coef_term is None:
        coef_term = 'd_div.i1.Y3( pp1, Pi1 )'

    get_state = problem.variables.get_state_part_view
    for irr in range( dim ):
        for icr in range( dim ):
            pc = get_state( states_rs[irr,icr], 'pc' )
            problem.variables['pp1'].data_from_data( pc )

            for irc in range( dim ):
                for icc in range( dim ):
                    io = HDF5MeshIO( filenames_rs[irc,icc] )
                    omega = io.read_data( 0 )['state_u'].data

                    problem.variables['Pi1'].data_from_data( omega )
                    
                    val = eval_term_op( None, coef_term, problem )

                    coef[dim*irr+icr,dim*irc+icc] = 0.5 * val
    return coef

##
# c: 12.03.2007, r: 23.06.2008
def coef_h( problem, corrs_rs, filenames_rs, variables_all, coef_eqs, kwargs ):
    """Viscous fading memory coefficient."""
    coef_term = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2
    ts = TimeStepper( *HDF5MeshIO( filenames_rs[0,0] ).read_time_stepper() )
    coef = nm.zeros( (ts.n_step, sym, sym), dtype = nm.float64 )

    indx = corrs_rs.di.indx['pc']
    for ir, (irr, icr) in enumerate( iter_sym( dim ) ):
        io = HDF5MeshIO( filenames_rs[irr,icr] )
        for step, time in ts:
            dpc = io.read_data( step )['dp'].data
            problem.variables['pp1'].data_from_data( dpc )

            for ic, (irc, icc) in enumerate( iter_sym( dim ) ):
                pc = corrs_rs.states_rs[irc,icc][indx]
                problem.variables['pp2'].data_from_data( pc )

                val = eval_term_op( None, coef_term, problem,
                                    call_mode = 'd_eval' )

                coef[step,ir,ic] = val
    return coef

##
# c: 12.04.2007, r: 15.04.2008
def coef_p( problem, corrs_alpha, pis, variables_all, coef_eqs, kwargs ):
    """Biot-like PBar coefficient."""
    coef_term_a, coef_term_b = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2
    coef = nm.zeros( (sym,), dtype = nm.float64 )

    indx = corrs_alpha.di.indx['uc']
    omega = corrs_alpha.state[indx]
    indx = corrs_alpha.di.indx['pc']
    pc = corrs_alpha.state[indx]

    problem.variables['pp1'].data_from_data( pc )
    problem.variables['Pi1'].data_from_data( omega )
    for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
        problem.variables['Pi2'].data_from_data( pis[ir,ic] )

        val1 = eval_term_op( None, coef_term_a, problem, call_mode = 'd_eval' )
        if ir == ic:
            val2 = eval_term_op( None, coef_term_b, problem )
        else:
            val2 = 0.0
        coef[ii] = val2 - val1
    return coef

##
# c: 12.04.2007, r: 23.06.2008
def coef_r( problem, filename_alpha, pis, variables_all, coef_eqs, kwargs ):
    """Biot-like RTilde coefficient."""
    coef_term_a, coef_term_b = coef_eqs[kwargs['term']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    sym = (dim + 1) * dim / 2
    io = HDF5MeshIO( filename_alpha )
    ts = TimeStepper( *io.read_time_stepper() )
    coef = nm.zeros( (ts.n_step, sym), dtype = nm.float64 )

    for step, time in ts:
        data = io.read_data( step )
        omega = data['u'].data
        pc = data['dp'].data

        problem.variables['pp1'].data_from_data( pc )
        problem.variables['Pi1'].data_from_data( omega )
        for ii, (ir, ic) in enumerate( iter_sym( dim ) ):
            problem.variables['Pi2'].data_from_data( pis[ir,ic] )

            val1 = eval_term_op( None, coef_term_a, problem,
                                 call_mode = 'd_eval' )
            if ir == ic:
                val2 = eval_term_op( None, coef_term_b, problem )
            else:
                val2 = 0.0
            coef[step,ii] = val2 - val1
    return coef

##
# c: 24.04.2007, r: 14.04.2008
def coef_g_star( problem, corrs_alpha, variables_all, coef_eqs, kwargs ):
    """Barenblatt coefficient for incompressible interface Y_3."""
    coef_term = coef_eqs[kwargs['term']]
    region_name = kwargs['region']

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    coef = nm.zeros( (1,), dtype = nm.float64 )

    indx = corrs_alpha.di.indx['uc']
    omega = corrs_alpha.state[indx]
    problem.variables[var_names[0]].data_from_data( omega )

    ct = coef_term % region_name
    coef[0] = eval_term_op( None, ct, problem )

    return coef

##
# c: 02.05.2007, r: 13.02.2008
def eval_boundary_diff_vel_grad( problem, uc, pc, equation, region_name,
                             pi = None ):
    set_state = problem.variables.set_state_part
    get_state = problem.variables.get_state_part_view

    problem.set_equations( {'eq_2' : equation} )

    state = problem.create_state_vector()
    set_state( state, uc, 'uc' )
    set_state( state, pc, 'pc' )
    if pi is not None:
        problem.variables['Pi'].data_from_data( pi )

#    problem.time_update( conf_ebc = {} )
    problem.time_update( conf_ebc = {}, conf_epbc = {} )
    
    problem.apply_ebc( state )
    aux = problem.get_evaluator().eval_residual( state )[0]
    pc = get_state( aux, 'pc', True )
    pc = problem.variables.make_full_vec( pc, 'pc', 0 )

    field = problem.variables['pc'].field

    reg = problem.domain.regions[region_name]
    nods = reg.get_field_nodes( field, merge = True )
    val = pc[nods].sum()

#    assert nm.all( problem.variables.di.ptr == problem.variables.adi.ptr )
    problem.time_update() # Restore EBC.

    return val

##
# c: 24.04.2007, r: 23.06.2008
def coef_g_plus( problem, filename_alpha, ebcs_all, epbcs_all, variables_all,
               coef_eqs, kwargs ):
    """Standard Barenblatt coefficient."""
    coef_term_a, coef_term_b = coef_eqs[kwargs['term']]
    region_name = kwargs['region']
    equation = coef_eqs[kwargs['aux_eq']]

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )
    problem.set_variables( variables )

    io = HDF5MeshIO( filename_alpha )
    ts = TimeStepper( *io.read_time_stepper() )
    coef = nm.zeros( (ts.n_step, 1), dtype = nm.float64 )

    for step, time in ts:
        data = io.read_data( step )
        omega, pc = data['u'].data, data['p'].data

        problem.variables['Pi1'].data_from_data( omega )
        ct = coef_term_a % region_name
        val1 = eval_term_op( None, ct, problem )
        val2 = eval_boundary_diff_vel_grad( problem, omega, pc, equation,
                                        'EBC' + region_name )

        problem.variables['pp1'].data_from_data( pc )
        val3 = 0
        val3 = eval_term_op( None, coef_term_b % region_name, problem )

        print '->>>', val1, val2, val3
        coef[step,0] = val1 + val2
        
    return coef

##
# c: 24.01.2008, r: 14.04.2008
def coef_g_bar( problem, ebcs_all, epbcs_all, variables_all,
              coef_eqs, kwargs ):
    """
    Note:

    solving "dw_hdpm_d.i1.Y3( m.K, qc, pc ) = 0" solve, in fact
    "C p^{\infty} = \hat{C} \hat{\pi}" with the result "\hat{p^{\infty}}",
    where the rhs comes from E(P)BC.
    - it is preferable to computing directly by
    "\hat{p^{\infty}} = \hat{C^-1 \strip(\hat{C} \hat{\pi})}", as it checks
    explicitly the rezidual.
    """
    
    mtx_term = coef_eqs[kwargs['aux_eq']]
    region_name = kwargs['region']

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names )

    problem.set_variables( variables )
    problem.set_equations( {'eq': mtx_term} )

    ebc_names = kwargs['ebcs']
    ebc = select_by_names( ebcs_all, ebc_names )
    epbc = select_by_names( epbcs_all, kwargs['epbcs'] )

    problem.time_update( conf_ebc = {}, conf_epbc = {} )
    dummy_f = problem.create_state_vector()
    mtx_cf = eval_term_op( dummy_f, mtx_term, problem,
                        dw_mode = 'matrix', tangent_matrix = problem.mtx_a )
    
    coef = nm.zeros( (1,), dtype = nm.float64 )

    problem.time_update( conf_ebc = ebc, conf_epbc = epbc, create_matrix = True )
    vec_pi_inf_f = problem.solve()
    vec_g_bar = mtx_cf * vec_pi_inf_f

    reg = problem.domain.regions[region_name]
    field = problem.variables[var_names[1]].field
    nods = reg.get_field_nodes( field, merge = True )
    coef[0] = vec_g_bar[nods].sum()

    return coef

##         vec_pi_hat = problem.create_state_vector()
##         problem.apply_ebc( vec_pi_hat )
## ##         problem.save_ebc( 'ebc_%d.vtk' % alpha,
## ##                          force = False, default = -1.0 )

##         mtx_c = eval_term_op( vec_pi_hat, mtx_term, problem,
##                            dw_mode = 'matrix', tangent_matrix = None )
##         print mtx_c.shape
##         vec_gf = - mtx_cf * vec_pi_hat
##         vec_g = problem.variables.strip_state_vector( vec_gf )
##         print vec_gf.shape
##         print vec_g.shape
        
##         ##
##         # \pi_{\infty}.
##         ls = Solver.any_from_conf( problem.get_solver_conf( kwargs['solver'] ),
##                                  mtx = mtx_c )
##         vec_pi_inf = ls( vec_g )
##         vec_pi_inf_f = problem.variables.make_full_vec( vec_pi_inf )

        
#        vec_pi_hat = problem.create_state_vector()
#        problem.apply_ebc( vec_pi_hat )
#        print problem.variables.has_ebc( vec_pi_hat )
#        problem.variables[var_names[2]].data_from_data( vec_pi_hat )

##         vec_pi_inf2 = problem.variables.strip_state_vector( vec_pi_inf_f2,
##                                                         follow_epbc = False )

##
# c: 24.04.2007, r: 23.06.2008
def coef_q_hat( problem, states_rs, filenames_rs, pis, region_name,
              coef_term = None ):
    """Biot-like coefficient for impermeable-like interface."""
    dim = problem.domain.mesh.dim
    coef = nm.zeros( (dim, dim), dtype = nm.float64 )
    
    if coef_term is None:
        coef_term = 'd_surface_integrate.isurf.%s( Pi1 )' % region_name

    equation = """dw_hdpm_d.i1.Y3( m.K, qc, pc ) + dw_div.i1.Y3( qc, uc )
                  + dw_div.i1.Y3( qc, Pi ) = 0"""

#    print filenames_rs
    get_state = problem.variables.get_state_part_view
    for ir in range( dim ):
        for ic in range( dim ):
            data = HDF5MeshIO( filenames_rs[ir,ic] ).read_data( 0 )
            omega = data['state_u'].data
            omega_bar = get_state( states_rs[ir,ic], 'uc' )

            problem.variables['Pi1'].data_from_data( omega + omega_bar )

            val1 = eval_term_op( None, coef_term, problem )

            pc = data['state_p'].data
            pc_bar = get_state( states_rs[ir,ic], 'pc' )
            val2 = eval_boundary_diff_vel_grad( problem, omega, pc, equation,
                                            'EBC' + region_name, pis[ir,ic] )
            val3 = eval_boundary_diff_vel_grad( problem, omega_bar, pc_bar, equation,
                                            'EBC' + region_name, pis[ir,ic] )

            print val1, val2, val3
            coef[ir,ic] = val1 + val2 + val3
    return coef

##
# c: 23.04.2007, r: 13.06.2008
def solve_permeability( problem, ebcs_all, epbcs_all, variables_all,
                       coef_eqs, kwargs ):
    """variable names:
         auxiliary problem: ppc, pqc
         coefficient: pp1, pp2"""
    region_name = kwargs['region']
    if kwargs.has_key( 'variant' ):
        variant = kwargs['variant']
        fld_replace = ('field', variant)
        reg_replace = ('region', region_name)
    else:
        variant = None
        fld_replace, reg_replace = None, None

    ##
    # Auxiliary permeability problem.
    auxpb = problem.copy( share = ['domain', 'conf', 'fields', 'materials'] )
    
    dim = auxpb.domain.mesh.dim

    equations = {
        'eq' : coef_eqs[kwargs['aux_eq']] % ((region_name, ) * 2)
    }

    ebc = select_by_names( ebcs_all, kwargs['ebcs'], replace = reg_replace )
    epbc = select_by_names( epbcs_all, kwargs['epbcs'] )

    var_names = kwargs['variables']
    variables = select_by_names( variables_all, var_names,
                               replace = fld_replace )

    auxpb.set_variables( variables )
    auxpb.set_equations( equations )
        
    auxpb.time_update( conf_ebc = ebc, conf_epbc = epbc )
##     auxpb.save_ebc( '_ebc.vtk', force = False, default = -1.0 )
##     pause()

    states_eta = []
    for ir in range( dim ):
        state = auxpb.solve( ir = ir )
##         auxpb.save_state( 'perm_%d_%s.vtk' % (ir, region_name), state )
        assert_( auxpb.variables.has_ebc( state ) )
        states_eta.append( state.copy() )
##     pause()

    ##
    # Permeability.
    coef_term = coef_eqs[kwargs['term']] % region_name
    output( coef_term )
    
    coef = nm.zeros( (dim, dim), dtype = nm.float64 )

    pressure = auxpb.variables[var_names[1]]
    coor = pressure.field.get_coor()
##     print coor
##     print coor.shape
##     pause()

    for ir in range( dim ):
        pp1 = coor[:,ir] + states_eta[ir]
        auxpb.variables[var_names[2]].data_from_data( pp1 )
        for ic in range( dim ):
            pp2 = coor[:,ic] + states_eta[ic]
            auxpb.variables[var_names[3]].data_from_data( pp2 )

            val = eval_term_op( None, coef_term, auxpb, call_mode = 'd_eval' )
            coef[ir,ic] = val

    return coef

##
# c: 29.08.2007, r: 06.05.2008
def get_matrix_parts( problem, ebc, epbc, variables_all,
                    conf_equations, mode = 'explicit' ):
    variables = select_by_names( variables_all, ['uc', 'vc', 'pc', 'qc'] )
    problem.set_variables( variables )
    indx = problem.variables.get_indx

    matrices = {}
    for key, mtx_term in conf_equations.iteritems():
        ks = key.split( ',' )
        mtx_name, var_names = ks[0], ks[1:]
        output( mtx_name, var_names )

        problem.set_equations( {'eq': mtx_term} )
        problem.time_update( conf_ebc = ebc, conf_epbc = epbc )

        ir = indx( var_names[0], stripped = True, allow_dual = True )
        ic = indx( var_names[1], stripped = True, allow_dual = True )

        dummy = problem.create_state_vector()
        mtx = eval_term_op( dummy, mtx_term, problem, dw_mode = 'matrix' )
##         print mtx
##         print ir, ic
##         if ir == ic:
##             aux = mtx[ir,ic].toarray()
##             ee = nm.sort( eig( aux )[0] )
##             print nla.det( aux )
##             print nm.prod( ee )
##             print ee
##         pause()
##         import pylab as p
##         import sfepy.base.plotutils as plu
##         plu.spy( mtx, eps = 1e-16 )
##         p.show()
##         debug()
        matrices[mtx_name] = mtx[ir,ic]


    mtx_a = matrices['mtx_a']
    mtx_bt = matrices['mtx_bt']
    output( 'full A size: %.3f MB' % (8.0 * nm.prod( mtx_a.shape ) / 1e6) )
    output( 'full B size: %.3f MB' % (8.0 * nm.prod( mtx_bt.shape ) / 1e6) )
##    debug()
##     mtx_bt = matrices['mtx_bt'].toarray()
##     mtx_b = matrices['mtx_b'].toarray()
##     assert nm.allclose( mtx_bt.T,  mtx_b )

    ls = Solver.any_from_conf( problem.ls_conf, presolve = True, mtx = mtx_a )
    if mode == 'explicit':
        tt = time.clock()
        mtx_aibt = nm.zeros( mtx_bt.shape, dtype = mtx_bt.dtype )
        for ic in xrange( mtx_bt.shape[1] ):
            mtx_aibt[:,ic] = ls( mtx_bt[:,ic].toarray().squeeze() )
        output( 'mtx_aibt: %.2f s' % (time.clock() - tt) )
        action_aibt = MatrixAction.from_array( mtx_aibt )
    else:
        ##
        # c: 30.08.2007, r: 13.02.2008
        def fun_aibt( vec ):
            # Fix me for sparse mtx_bt...
            rhs = sc.dot( mtx_bt, vec )
            out = ls( rhs )
            return out
        action_aibt = MatrixAction.from_function( fun_aibt,
                                                (mtx_a.shape[0],
                                                 mtx_bt.shape[1]),
                                                nm.float64 )
    matrices['action_aibt'] = action_aibt

    return matrices
    
##
# c: 11.04.2007, r: 09.04.2008
def solve_pressure_eigenproblem( matrices, eig_problem, n_eigs = 0, check = False ):
    """G = B*AI*BT or B*AI*BT+D"""

    def get_slice( n_eigs, nn ):
        if n_eigs > 0:
            ii = slice( 0, n_eigs )
        elif n_eigs < 0:
            ii = slice( nn + n_eigs, nn )
        else:
            ii = slice( 0, 0 )
        return ii

    ms = matrices
    mtx_c, mtx_b, action_aibt = ms['mtx_c'], ms['mtx_b'], ms['action_aibt']
    mtx_g = mtx_b * action_aibt.to_array() # mtx_b must be sparse!
    if eig_problem == 'B*AI*BT+D':
        mtx_g += ms['mtx_d'].toarray()

    ms['mtx_g'] = mtx_g
    output( mtx_c.shape, mtx_g.shape )

    eigs, mtx_q = eig( mtx_c.toarray(), mtx_g, method = 'eig.sgscipy' )

    if check:
        ee = nm.diag( sc.dot( mtx_q.T * mtx_c, mtx_q ) ).squeeze()
        oo = nm.diag( sc.dot( sc.dot( mtx_q.T,  mtx_g ), mtx_q ) ).squeeze()
##         print ee
##         print eigs
##         print oo
        try:
            assert_( nm.allclose( ee, eigs ) )
            assert_( nm.allclose( oo, nm.ones_like( eigs ) ) )
        except ValueError:
            debug()

    nn = mtx_c.shape[0]
    if isinstance( n_eigs, tuple ):
        output( 'required number of eigenvalues: (%d, %d)' % n_eigs )
        if sum( n_eigs ) < nn:
            ii0 = get_slice( n_eigs[0], nn )
            ii1 = get_slice( -n_eigs[1], nn )
            eigs = nm.concatenate( (eigs[ii0], eigs[ii1] ) )
            mtx_q = nm.concatenate( (mtx_q[:,ii0], mtx_q[:,ii1]), 1 ) 
    else:
        output( 'required number of eigenvalues: %d' % n_eigs )
        if (n_eigs != 0) and (abs( n_eigs ) < nn):
            ii = get_slice( n_eigs, nn )
            eigs = eigs[ii]
            mtx_q = mtx_q[:,ii]

##     from sfepy.base.plotutils import pylab, iplot
##     pylab.semilogy( eigs )
##     pylab.figure( 2 )
##     iplot( eigs )
##     pylab.show()
##     debug()

    out = Struct( eigs = eigs, mtx_q = mtx_q, matrices = matrices )
    return out

def make_t_correctors_via_evp( problem, p_evp, state0, ts, file_conf,
                               alpha = None,
                               post_process_hook = None, file_per_var = False ):
    format = '====== time %%e (step %%%dd of %%%dd) =====' % ((ts.n_digit,) * 2)

    dump_filename = file_conf

    get_state = problem.variables.get_state_part_view
    make_full_vec = problem.variables.make_full_vec
    to_output = problem.variables.state_to_output

    nr, nc = p_evp.mtx_q.shape

    vec_g = None
    if alpha is not None:
        # Nonzero EBC case.
        auxpb = problem.copy( share = ['domain', 'conf', 'fields', 'materials'] )
        variables = select_by_names( auxpb.conf.variables, ['pc', 'qc'] )
        auxpb.set_variables( variables )
        auxpb.set_equations( {'eq' : 'dw_diffusion.i1.Y3( m.K, qc, pc )'}  )
        auxpb.time_update( conf_ebc = p_evp.ebc, conf_epbc = p_evp.epbc )

        state = auxpb.create_state_vector()
        auxpb.apply_ebc( state )

        vec_g = eval_term_op( state, "- dw_diffusion.i1.Y3( m.K, qc, pc )",
                           auxpb, dw_mode = 'vector' )
#        print vec_g
        if nm.allclose( vec_g, 0.0 ):
            vec_g = None
        else:
            one = nm.ones( (nc,), dtype = nm.float64 )

    if vec_g is not None:
        output( 'nonzero pressure EBC' )

    ##
    # follow_epbc = False -> R1 = - R2 as required. ? for other correctors?
    sstate0 = problem.variables.strip_state_vector( state0, follow_epbc = False )
    vec_p0 = get_state( sstate0, 'pc', True )
##     print state0
##     print vec_p0
##     print vec_p0.min(), vec_p0.max(), nla.norm( vec_p0 )
##     debug()

    # xi0 = Q^{-1} p(0) = Q^T G p(0)
    vec_xi0 = sc.dot( p_evp.mtx_q.T,
                     sc.dot( p_evp.matrices['mtx_g'],
                             vec_p0[:,nm.newaxis] ) ).squeeze()
    action_aibt = p_evp.matrices['action_aibt']

    e_e_qg = 0.0
    iee_e_qg = 0.0
    for step, time in ts:
        output( format % (time, step + 1, ts.n_step) )

        e_e = nm.exp( - p_evp.eigs * time )
        e_e_qp = e_e * vec_xi0 # exp(-Et) Q^{-1} p(0)

        if vec_g is not None:
            Qg = sc.dot( p_evp.mtx_q.T, vec_g )
            e_e_qg = e_e * Qg
            iee_e_qg = ((one - e_e) / p_evp.eigs) * Qg

        vec_p = sc.dot( p_evp.mtx_q, e_e_qp + iee_e_qg )
        vec_dp = - sc.dot( p_evp.mtx_q, (p_evp.eigs * e_e_qp - e_e_qg) )
        vec_u = action_aibt( vec_dp )
##         bbb = sc.dot( vec_dp.T, - p_evp.matrices['mtx_c'] * vec_p0 )

        vec_u = make_full_vec( vec_u, 'uc', None )
        vec_p = make_full_vec( vec_p, 'pc', None )
        vec_dp = make_full_vec( vec_dp, 'pc', 0.0 ) # time derivative of constant!
##         aaa = sc.dot( vec_xi0.T, p_evp.eigs * (p_evp.eigs * e_e_qp) )
##         print aaa
##         print bbb
        out = {'u' : Struct( name = 'dump', mode = 'nodes', data = vec_u,
                             dofs = None, var_name = 'uc' ),
               'p' : Struct( name = 'dump', mode = 'nodes', data = vec_p,
                             dofs = None, var_name = 'pc' ),
               'dp' : Struct( name = 'dump', mode = 'nodes', data = vec_dp,
                              dofs = None, var_name = 'pc' )}
        problem.save_state( dump_filename, out = out, file_per_var = False,
                            ts = ts )
        # For visualization...
        out = {}
        extend = not file_per_var
        out.update( to_output( vec_u, var_info = {'uc' : (True, 'uc')},
                               extend = extend ) )
        out.update( to_output( vec_p, var_info = {'pc' : (True, 'pc')},
                               extend = extend ) )
        out.update( to_output( vec_dp, var_info = {'pc' : (True, 'dpc')},
                               extend = extend ) )
        if post_process_hook is not None:
            out = post_process_hook( out, problem,
                                     {'u' : vec_u, 'p' : vec_p, 'dp' : vec_dp},
                                     extend = extend )
        problem.save_state( dump_filename, out = out,
                            file_per_var = file_per_var, ts = ts )

def verify_t_correctors( problem, ebc, epbc,
                         conf_equations, initial_state, dump_filename ):
    variables = select_by_names( problem.conf.variables,
                               ['uc', 'vc', 'pc', 'qc', 'pp1'] )
    problem.set_variables( variables )
    problem.set_equations( conf_equations )

    io = HDF5MeshIO( dump_filename )
    ts = TimeStepper( *io.read_time_stepper() )

    get_state = problem.variables.get_state_part_view
    p0 = get_state( initial_state, 'pc' )

    format = '====== time %%e (step %%%dd of %%%dd) =====' % ((ts.n_digit,) * 2)
    vv = problem.variables
    ok = True
    for step, time in ts:
        output( format % (time, step + 1, ts.n_step) )

        data = io.read_data( step )
        if step == 0:
            assert_( nm.allclose( data['p'].data, p0 ) )

        problem.time_update( conf_ebc = ebc, conf_epbc = epbc )

        state0 = problem.create_state_vector()
        state0[vv.di.indx['uc']] = data['u'].data
        state0[vv.di.indx['pc']] = data['p'].data
        vv['pp1'].data_from_data( data['dp'].data )

        state = problem.solve( state0 = state0, ts = ts )
        err = nla.norm( state - state0 )
        print state.min(), state.max()
        print state0.min(), state0.max()
        print '>>>>>', err

        ok = ok and (err < 1e-15)
        problem.advance( ts )

    return ok

##
# c: 24.01.2008, r: 11.04.2008
def get_ebc_gamma( ebc_in, ebc_names, alpha ):
    replace = ('region', {'fixed_p_zero' : 'EBCGammaY%d' % (3 - alpha),
                          'fixed_p_one' : 'EBCGammaY%d' % (alpha)})
    ebc = select_by_names( ebc_in, ebc_names,
                         replace = replace,
                         simple = False )
    return ebc

##
# c: 13.02.2008, r: 27.02.2008
def create_pis( problem, variables_all, var_name ):
    variables = select_by_names( variables_all, [var_name] )
    problem.set_variables( variables )

    dim = problem.domain.mesh.dim
    pis = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = build_op_pi( var_name, problem, ir, ic )
            pis[ir,ic] = pi
    return pis


##
# c: 14.04.2008, r: 14.04.2008
def get_alpha_equations( conf_equations, alpha ):
    equations = {}
    for key, val in conf_equations.iteritems():
        if key == 'eq_1':
            val = val % alpha
        equations[key] = val
    return equations

def save_steady_correctors( base_name, state, problem,
                            post_process_hook, file_per_var ):
    get_state = problem.variables.get_state_part_view

    problem.save_state( base_name + '.vtk', state,
                        post_process_hook = post_process_hook,
                        file_per_var = file_per_var )

    out = {'u' : Struct( name = 'dump', mode = 'nodes',
                         data = get_state( state, 'uc' ),
                         dofs = None, var_name = 'uc' ),
           'p' : Struct( name = 'dump', mode = 'nodes',
                         data = get_state( state, 'pc' ),
                         dofs = None, var_name = 'pc' )}
    problem.save_state( base_name + '.h5', out = out,
                        file_per_var = False )

    
def  solve_steady_correctors_alpha( problem, ebcs_all, epbcs_all, variables_all,
                                    conf_equations, ebc_names, epbc_names,
                                    file_conf, alpha = 0,
                                    post_process_hook = None,
                                    file_per_var = False ):
    """
    variant .. int -> alpha: 0    : both
                             1, 2 : 1 or 2
    """
    ebc = select_by_names( ebcs_all, ebc_names )
    epbc = select_by_names( epbcs_all, epbc_names )
    variables = select_by_names( variables_all, ['uc', 'vc', 'pc', 'qc'] )

    problem.set_variables( variables )
    eqs = get_alpha_equations( problem.conf.equations_steady_alpha, alpha )
    problem.set_equations( eqs )
    problem.time_update( conf_ebc = ebc, conf_epbc = epbc )
##     problem.save_ebc( ofn_trunk + '_ebc_%d.vtk' % alpha,
##                      force = False, default = -1.0 )
##     problem.save_ebc( ofn_trunk + '_ebc2_%d.vtk' % alpha )

    state = problem.solve()
    assert_( problem.variables.has_ebc( state ) )

    save_steady_correctors( file_conf % alpha, state, problem,
                            post_process_hook, file_per_var )

    return Struct( name = 'Steady alpha correctors',
                   state = state,
                   alpha = alpha,
                   di = problem.variables.di )

def  solve_steady_correctors_rs( problem, ebcs_all, epbcs_all, variables_all,
                                 conf_equations, ebc_names, epbc_names, pis,
                                 file_conf, post_process_hook = None,
                                 file_per_var = False ):
    dim = problem.domain.mesh.dim

    ebc = select_by_names( ebcs_all, ebc_names )
    epbc = select_by_names( epbcs_all, epbc_names )
    variables = select_by_names( variables_all,
                                 ['uc', 'vc', 'pc', 'qc', 'Pi'] )

    problem.set_variables( variables )
    problem.set_equations( conf_equations )

    problem.time_update( conf_ebc = ebc, conf_epbc = epbc )

    states_rs = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = pis[ir,ic]
            problem.variables['Pi'].data_from_data( pi )

            state = problem.create_state_vector()
            problem.apply_ebc( state )
            state = problem.solve()
            assert_( problem.variables.has_ebc( state ) )
            states_rs[ir,ic] = state

            save_steady_correctors( file_conf % (ir, ic), state, problem,
                                    post_process_hook, file_per_var )

    return Struct( name = 'Steady RS correctors',
                   states_rs = states_rs,
                   di = problem.variables.di )

def  solve_steady_correctors_pressure( problem, ebcs_all, epbcs_all,
                                       variables_all,
                                       conf_equations, ebc_names, epbc_names,
                                       file_conf, post_process_hook = None,
                                       file_per_var = False ):
    dim = problem.domain.mesh.dim

    ebc = select_by_names( ebcs_all, ebc_names )
    epbc = select_by_names( epbcs_all, epbc_names )
    variables = select_by_names( variables_all,
                                 ['uc', 'vc', 'pc', 'qc',
                                  'p_y_one', 'p_ym_one'] )

    problem.set_variables( variables )
    problem.set_equations( conf_equations )

    problem.time_update( conf_ebc = ebc, conf_epbc = epbc )

    vv = problem.variables

    one = nm.ones( (vv['p_y_one'].field.n_nod,), dtype = nm.float64 )
    vv['p_y_one'].data_from_data( one )

    one_m = one[:vv['p_ym_one'].field.n_nod]
    vv['p_ym_one'].data_from_data( one_m )

    state = problem.create_state_vector()
    problem.apply_ebc( state )
    state = problem.solve()
    assert_( problem.variables.has_ebc( state ) )

    save_steady_correctors( file_conf, state, problem,
                            post_process_hook, file_per_var )

    return Struct( name = 'Steady pressure correctors',
                   state = state,
                   di = problem.variables.di )
