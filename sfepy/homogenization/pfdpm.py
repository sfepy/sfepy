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
# c: 24.01.2008, r: 11.04.2008
def get_ebc_gamma( ebc_in, ebc_names, alpha ):
    replace = ('region', {'fixed_p_zero' : 'EBCGammaY%d' % (3 - alpha),
                          'fixed_p_one' : 'EBCGammaY%d' % (alpha)})
    ebc = select_by_names( ebc_in, ebc_names,
                         replace = replace,
                         simple = False )
    return ebc

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

