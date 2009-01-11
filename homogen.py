#!/usr/bin/env python
# 11.07.2006, c 
import os.path as op
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.evaluate import eval_term_op
import sfepy.homogenization.pfdpm as pfdpm
from sfepy.homogenization.prolong import ProlongatorToQP
from sfepy.base.conf import get_standard_keywords
from sfepy.solvers.generic import solve_evolutionary_op, solve_stationary_op

##
# 19.09.2006, c
# 21.09.2006
# 11.06.2007
def post_process( state, problem ):

    dvel1 = eval_term_op( state, 'de_hdpm_dvel.Omega( m.K1, p1 )', problem,
                        new_geometries = False )
    dvel2 = eval_term_op( state, 'de_hdpm_dvel.Omega( m.K2, p2 )', problem,
                        new_geometries = False )
    p1grad = eval_term_op( state, 'de_hdpm_pgrad.Omega( p1 )', problem,
                         new_geometries = False )
    p2grad = eval_term_op( state, 'de_hdpm_pgrad.Omega( p2 )', problem,
                         new_geometries = False )
    strain = eval_term_op( state, 'de_sdcc_strain.Omega( u )', problem,
                         new_geometries = False )

    out = problem.state_to_output( state )
    out['dvel1'] = Struct( name = 'output_data',
                           mode = 'cell', data = dvel1,
                           dof_types = None )
    out['dvel2'] = Struct( name = 'output_data',
                           mode = 'cell', data = dvel2,
                           dof_types = None )
    out['p1grad'] = Struct( name = 'output_data',
                            mode = 'cell', data = p1grad,
                            dof_types = None )
    out['p2grad'] = Struct( name = 'output_data',
                            mode = 'cell', data = p2grad,
                            dof_types = None )
    out['strain'] = Struct( name = 'output_data',
                            mode = 'cell', data = strain,
                            dof_types = None )
    out['p2-p1'] = copy( out['p1'] )
    out['p2-p1'].data = out['p2'].data - out['p1'].data
    out.update( post_process_force_load( 'force', 'u', 'v', state, problem ) )

    return out

##
# 21.09.2006, c
# 11.06.2007
def post_process_force_load( force_name, var, test, state, problem ):

    force = eval_term_op( state,
                        'dw_volume_lvf.Omega( %s, %s )' % (force_name, test),
                        problem, new_geometries = False )
    aux = problem.create_state_vector()
    problem.update_vec( aux, -force )
    out1 = problem.state_to_output( aux )
    out = {}
    out[force_name] = Struct( name = 'output_data',
                             mode = 'vertex',
                             data = problem.materials[force_name].data,
                             dof_types = (0, 1) )
    out[force_name+'2'] = out1[var]

    return out

##
# 19.09.2006, c
# 20.09.2006
# 22.09.2006
# 25.09.2006
# 26.09.2006
# 27.09.2006
# 01.12.2006
# 05.12.2006
# 07.06.2007
# 08.06.2007
def solve_evolutionary( conf, options ):
    if conf.nls.problem == 'linear':
        print 'time-dependent term caches require final residual evaluation!'
        raise NotImplementedError

    problem = ProblemDefinition.from_conf( conf )
    problem.equations.set_cache_mode( conf.fe.cache_override )

##     for cc in problem.equations.caches.values():
##         print cc
##     pause()

    td = {'state_history' : []}
    nls_context = Struct( problem = problem, data = td )

    ts = TimeStepper.from_conf( conf.ts )
    n_digit = int( nm.log10( ts.n_step - 1 ) + 1 )
    format = '====== time %%e (step %%%dd of %%%dd) =====' % (n_digit, n_digit)
#    print n_digit, format
#    pause()

    state = problem.create_state_vector()

    if options.output_filename_trunk:
        ofn_trunk = options.output_filename_trunk
    else:
        ofn_trunk = get_trunk( conf.filename_mesh ) + '_out'

    if conf.save_format == 'hdf5':
        output_filename = 'output/' + ofn_trunk + '.h5'
        write_hdf5( output_filename, mesh = problem.domain.mesh, ts = ts )
    else:
        suffix = '.%%0%dd.vtk' % n_digit
#        print suffix

    is_save = nm.linspace( 0, ts.n_step - 1, conf.save_steps ).astype( int )
    is_save = nm.unique1d( is_save )
    ii = 0

    state0 = state.copy()

    for step, time in ts:
        output( format % (time, step + 1, ts.n_step) )

        td.update( {'ts' : ts} )
        problem.time_update( ts )
        problem.apply_ebc( state )

        if step == 0:
            problem.save_ebc( get_trunk( conf.filename_mesh ) + '_ebc.vtk' )

        outs = ['u_0', 'p1_0', 'p2_0']
        ins = ['u', 'p1', 'p2']
        problem.variables.non_state_data_from_state( outs, state0, ins )

        ##
        # Step 0 has zero out-of-balance and no history terms...
        if step == 1:
            if conf.nls.matrix == 'external':
                problem.mtx_a, status = tangent_matrix( state, problem.mtx_a,
                                                      nls_context )
                if status != 0:
                    raise RuntimeError, 'tangent matrix computation failed'

        if step == 2:
            conf.nls.check = 0
##        problem.regions.print_names()
##        problem.equations.print_names()
##        pause()
        state = newton( conf.nls, state, problem.mtx_a, nls_context )

        if is_save[ii] == step:
            out = post_process( state, problem )
##             problem.set_equations( conf.equations,
##                                   cache_override = conf.fe.cache_override )

            if conf.save_format == 'hdf5':
                write_hdf5( output_filename, ts = ts, out = out )
            else:
                fd = open( 'output/' + ofn_trunk + suffix % step, 'w' )
                write_vtk( fd, problem.domain.mesh, out )
                fd.close()

            ii += 1

        state0 = state.copy()
        problem.advance( ts )

        for cc in problem.equations.caches.values():
            print cc.name, cc.mem_sizes

##
# 19.09.2006, c
# 20.09.2006
# 13.06.2007
def solve_steady( conf, options ):
    ts = TimeStepper.from_conf( conf.ts )
    ts.set_step( None, conf.steady_state_n_time )
    n_digit = int( nm.log10( ts.n_step - 1 ) + 1 )
    format = '====== time %%e (step %%%dd of %%%dd) =====' % (n_digit, n_digit)

    output( format % (ts.time, ts.step + 1, ts.n_step) )

    problem = ProblemDefinition.from_conf( conf, init_variables = False )

    ##
    # Pressure subproblem.
    p_keys = ['qh1', 'ph1', 'qh2', 'ph2']
    conf_variables = {}
    for key in p_keys:
        conf_variables[key] = conf.variables[key]

    print conf_variables

    problem.set_variables( conf_variables )
    problem.set_equations( conf.equations_steady_p )

    state = problem.create_state_vector()
    print state.shape

    problem.time_update( ts )
    problem.apply_ebc( state )

    problem.save_ebc( get_trunk( conf.filename_mesh ) + '_pb' + '_ebc.vtk' )

    nls_context = Struct( problem = problem, data = {} )
    state = newton( conf.nls, state, problem.mtx_a, nls_context )

    pb1 = problem.variables.get_state_part_view( state, 'ph1' )
    pb2 = problem.variables.get_state_part_view( state, 'ph2' )

    out1 = problem.state_to_output( state )

    ##
    # Displacement subproblem.
    u_keys = ['vh', 'uh', 'pb1', 'pb2']
    conf_variables = {}
    for key in u_keys:
        conf_variables[key] = conf.variables[key]

    print conf_variables

    problem.set_variables( conf_variables )
    problem.set_equations( conf.equations_steady_u )

    state = problem.create_state_vector()
    print state.shape

    problem.update_bc( ts, conf.ebc, conf.epbc, conf.funmod )
    problem.apply_ebc( state )

    problem.save_ebc( get_trunk( conf.filename_mesh ) + '_ub' + '_ebc.vtk' )

    problem.variables['pb1'].data_from_data( pb1, slice( 0, len( pb1 ) ) )
    problem.variables['pb2'].data_from_data( pb2, slice( 0, len( pb1 ) ) )

    nls_context = Struct( problem = problem, data = {'ts': ts} )
    state = newton( conf.nls, state, problem.mtx_a, nls_context )

    ub = problem.variables.get_state_part_view( state, 'uh' )

    out = problem.state_to_output( state )
    out['pb1'] = out1['ph1']
    out['pb2'] = out1['ph2']
    out.update( post_process_force_load( 'force_steady', 'uh', 'vh',
                                      state, problem ) )

    ##
    # Save complete steady state solution.
    if options.output_filename_trunk:
        ofn_trunk = options.output_filename_trunk
    else:
        ofn_trunk = get_trunk( conf.filename_mesh ) + '_out'

    fd = open( 'output/' + ofn_trunk + '_steady.vtk', 'w' )
    write_vtk( fd, problem.domain.mesh, out )
    fd.close()

    return problem, ub, pb1, pb2

##
# 20.09.2006, c
# 13.06.2007
def solve_impact( conf, options, problem, ub, pb1, pb2 ):

    ts = TimeStepper.from_conf( conf.ts )
    ts.set_step( None, conf.steady_state_n_time )
    
    problem.set_variables( conf.variables )
    problem.set_equations( conf.equations_impact )
    
    state = problem.create_state_vector()
    print state.shape

    problem.time_update( ts )
    problem.apply_ebc( state, force_values = 0.0 )

    problem.save_ebc( get_trunk( conf.filename_mesh ) + '_pb' + '_ebc.vtk' )

    problem.variables['ub'].data_from_data( ub, slice( 0, len( ub ) ) )
    problem.variables['pb1'].data_from_data( pb1, slice( 0, len( pb1 ) ) )
    problem.variables['pb2'].data_from_data( pb2, slice( 0, len( pb1 ) ) )

    nls_context = Struct( problem = problem, data = {'ts': ts} )
    state = newton( conf.nls, state, problem.mtx_a, nls_context )

    out = problem.state_to_output( state )
    out.update( post_process_force_load( 'force_impact', 'uh', 'vh',
                                      state, problem ) )
    
    ##
    # Save impact solution.
    if options.output_filename_trunk:
        ofn_trunk = options.output_filename_trunk
    else:
        ofn_trunk = get_trunk( conf.filename_mesh ) + '_out'

    fd = open( 'output/' + ofn_trunk + '_impact.vtk', 'w' )
    write_vtk( fd, problem.domain.mesh, out )
    fd.close()

    return state


##
# 28.08.2007, c
def solve_steady_correctors( problem, conf, options, ofn_trunk ):

    coefs = pfdpm.Coefficients()
        
    ##
    # Homogenized permeability.
    coefs.K1 = pfdpm.solve_permeability( problem, conf, 'Y1' )
    print coefs.K1
    coefs.K2 = pfdpm.solve_permeability( problem, conf, 'Y2' )
    print coefs.K2

##     if options.debug_micro:
##         mat = problem.materials['debug_macro'].get_data
##         pfdpm.debug_data( 'C1', mat( 'C1' ), coefs.K1 )
##         pfdpm.debug_data( 'C2', mat( 'C2' ), coefs.K2 )

    problem.set_variables( conf.variables )

    ##
    # Volume fractions based on 'uc' geometry.
    aux = pfdpm.volumes( problem, ['Y1', 'Y2', 'Y3', 'Y'] )
    coefs.volume_y1, coefs.volume_y2, coefs.volume_y3, coefs.volume_y = aux
    vv = coefs.volume_y1 + coefs.volume_y2 + coefs.volume_y3
    assert_( abs( vv - coefs.volume_y ) < 1e-15 )
    coefs.VF = nm.array( [coefs.volume_y1 / coefs.volume_y,
                          coefs.volume_y2 / coefs.volume_y], dtype = nm.float64 )
    ##
    # Solve steady_alpha.
    aux = pfdpm.solve_steady_correctors( problem,
                                       conf.equations_steady_alpha,
                                       ofn_trunk, conf, options.variant )
    corrs_alpha, alphas = aux 
##     vec_p = problem.variables.get_state_part_view( corrs_alpha[1], 'pc' )
##     print vec_p.min(), vec_p.max(), nla.norm( vec_p )
##     pause()
    ##
    # Solve steady_rs. 
    corrs_rs, pis = pfdpm.solve_steady_correctors( problem,
                                                conf.equations_steady_rs,
                                                ofn_trunk, conf )

    ##
    # Homogenized coefficients.
    coefs.GStar = pfdpm.coef_g_star( problem, alphas, corrs_alpha )
    
    coefs.E = pfdpm.coef_e( problem, corrs_rs, pis )
    ii = [0, 3, 1, 2]
    indx = nm.meshgrid( ii, ii )
    print coefs.E[indx]


    sol = Struct( corrs_alpha = corrs_alpha,
                  corrs_rs = corrs_rs,
                  alphas = alphas, pis = pis, indx = indx )
##     print coefs
##     print sol
##     pause()
    
    return sol, coefs

##
# 28.08.2007, c
# 29.08.2007
# 30.08.2007
# 01.09.2007
# 03.09.2007
def solve_time_variant_correctors( sol_s, problem, conf, options, ofn_trunk ):
    import time as _time

    coefs = pfdpm.Coefficients()

    ##
    # Eigensolutions for time-variant correctors.
    output( 'eigensolutions for time-variant...' )
    problem.set_equations( conf.equations_eig_time_variant )

    if conf.options.grid_mode == 'single_grid':
        if conf.options.eig_solver == 'dense':
            p_evp = pfdpm.solve_pressure_eigenproblem( problem,
                                                    conf.options.n_eigs )
        else:
            raise NotImplementedError
    else:
        cconf = copy( conf )
        cconf.filename_mesh = conf.options.filename_coarse_mesh
        cpb = ProblemDefinition.from_conf( cconf,
                                          init_equations = False )
        cpb.set_equations( conf.equations_eig_time_variant )
        if conf.options.eig_solver == 'dense':
            p_evp = pfdpm.solve_pressure_eigenproblem( cpb,
                                                    conf.options.n_eigs,
                                                    ret_matrices = False )
        else:
            raise NotImplementedError

        aux = pfdpm.get_matrix_parts( problem, conf.options.AIBTAction )
        p_evp.action_aibt = aux[2]

        prolong = ProlongatorToQP( problem, cpb, conf.options, 'qc', 'pc',
                                   debug = False )

        print p_evp.mtx_q.shape
        p_evp.mtx_q = prolong( p_evp.mtx_q )
        print p_evp.mtx_q.shape

    print p_evp.eigs
#    print p_evp.mtx_q
#    pause()

    ##
    # Solve time-variant correctors via eigensolutions.
    dim = problem.domain.mesh.dim

    ts = TimeStepper.from_conf( conf.ts )
    coefs.dt = ts.dt
    coefs.n_step = ts.n_step
    coefs.times = ts.times
    
    output( 'time-variant rs correctors via eigensolutions...' )
    filenames_rs = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            filename = ofn_trunk + '_t%d%d.h5' % (ir,ic)
            filenames_rs[ir,ic] = filename
            pfdpm.make_t_correctors_via_evp( problem, p_evp,
                                         -sol_s.corrs_rs[ir,ic], ts,
                                         filename, force_ebc = 0.0 )

    output( 'time-variant alpha correctors via eigensolutions...' )
    n_alpha = len( sol_s.alphas )
    filenames_alpha = nm.zeros( (n_alpha,), dtype = nm.object )
    for ii, alpha in enumerate( sol_s.alphas ):
        filename = ofn_trunk + '_t%d.h5' % alpha
        filenames_alpha[ii] = filename
        problem.time_update( conf_ebc = getattr( conf, 'ebc%d' % alpha ) )
        pfdpm.make_t_correctors_via_evp( problem, p_evp, sol_s.corrs_alpha[ii], ts,
                                     filename )
    ##
    # Solve time-variant correctors.
##     problem.set_equations( conf.equations_time_variant )


    # Alpha correctors.
##     pfdpm.solve_time_variant_correctors( problem, state_alpha, ts,
##                                       ofn_trunk + '_t_alpha.h5' )

##     # RS correctors.
##     filenames_rs = nm.zeros( (dim, dim), dtype = nm.object )
##     for ir in range( dim ):
##         for ic in range( dim ):
##             filename = ofn_trunk + '_t%d%d.h5' % (ir,ic)
##             filenames_rs[ir,ic] = filename
##             pfdpm.solve_time_variant_correctors( problem, -corrs_rs[ir,ic], ts,
##                                               filename, force_ebc = 0.0 )

    coefs.BBar = pfdpm.coef_b_bar( problem, sol_s.corrs_rs, filenames_rs )
    print coefs.BBar[solS.indx]

    coefs.PBar = pfdpm.coef_rp( problem, sol_s.alphas, sol_s.corrs_alpha,
                               filenames_alpha, sol_s.pis, 0, 'P' )
    print coefs.PBar

    coefs.QHat = nm.zeros( (n_alpha, dim, dim), dtype = nm.float64 )
    for ii, alpha in enumerate( sol_s.alphas ):
        rname = 'Gamma%d' % alpha
        coefs.QHat[ii] = pfdpm.coef_q_hat( problem, sol_s.corrs_rs, filenames_rs,
                                         sol_s.pis, rname )
    print coefs.QHat

    coefs.Hs = []
    coefs.RTildes = []
    coefs.GPluses = []
    for step, time in ts:
        print step
        tt = _time.clock()
        coef_h = pfdpm.coef_h( problem, sol_s.corrs_rs, filenames_rs, step )
        print _time.clock() - tt

#        print coef_h[sol_s.indx]
        coefs.Hs.append( coef_h )

        coef_r_tilde = pfdpm.coef_rp( problem, sol_s.alphas, sol_s.corrs_alpha,
                                   filenames_alpha, sol_s.pis, step, 'R' )
#        print coef_r_tilde
        coefs.RTildes.append( coef_r_tilde )

        coef_g_plus = pfdpm.coef_g_plus( problem, sol_s.alphas, filenames_alpha,
                                     step )
        print coef_g_plus
        coefs.GPluses.append( coef_g_plus )
#        pause()

    if options.debug_micro:
        mat = problem.materials['debug_macro'].get_data
        problem.materials.set_current_group( 0 )
        
        ip = [0, 3, 1]
        indx = nm.meshgrid( ip, ip )

        pfdpm.debug_data( 'E', mat( 'E' ), coefs.E[indx] )
        pfdpm.debug_data( 'BBar', mat( 'BBar' ), coefs.BBar[indx] )

        for step, time in ts:
            output( 'step', step )
            pfdpm.debug_data( 'H', mat( 'H' )[step], coefs.Hs[step][indx] )

        for ii, alpha in enumerate( sol_s.alphas ):
            pfdpm.debug_data( 'GStar%d' % alpha, mat( 'GStar' ),
                             coefs.GStar[ii] )
            pfdpm.debug_data( 'PBar%d' % alpha, mat( 'PBar%d' % alpha ),
                             coefs.PBar[ii].flat[ip] )
            pfdpm.debug_data( 'QHat%d' % alpha, mat( 'QHat%d' % alpha ),
                             coefs.QHat[ii].flat[ip] )

            for step, time in ts:
                output( 'step', step )
                sgn = (alpha == 1) and 1. or -1.
                pfdpm.debug_data( 'RTilde%d' % alpha,
                                 sgn * mat( 'RTilde1' )[step].squeeze(),
                                 coefs.RTildes[step][ii].flat[ip] )


        pause()

    sol = Struct( filenames_alpha = filenames_alpha,
                  filenames_rs = filenames_rs )
    return sol, coefs

##
# 12.02.2007, c
# 19.02.2007
# 28.02.2007
# 01.03.2007
# 02.03.2007
# 09.03.2007
# 12.03.2007
# 11.04.2007
# 12.04.2007
# 23.04.2007
# 24.04.2007
# 26.04.2007
# 02.05.2007
# 10.05.2007
# 15.06.2007
# 09.07.2007
# 16.07.2007
# 17.07.2007
# 28.08.2007
# 01.09.2007
# 03.09.2007
def solve_micro( conf, options ):

    if options.output_filename_trunk:
        ofn_trunk = options.output_filename_trunk
    else:
        ofn_trunk = get_trunk( conf.filename_mesh ) + '_out'
    ofn_trunk = 'output/' + ofn_trunk

    problem = ProblemDefinition.from_conf( conf,
                                          init_variables = False,
                                          init_equations = False )

    sol_s, coefs_s = solve_steady_correctors( problem, conf, options, ofn_trunk )

    if conf.options.check_eigs is not None:
        def pn( p_eigs, n_eigs_total ):
            if p_eigs == 100:
                return n_eigs_total
            elif p_eigs > 0:
                return nm.clip( int(0.01 * p_eigs * n_eigs_total), 0, n_eigs_total )
            else:
                return nm.clip( int(0.01 * p_eigs * n_eigs_total),
                                -n_eigs_total, -0 )
            
        indx = problem.variables.adi.indx
        ip = indx['pc']
        n_eigs_total = ip.stop - ip.start
        for ii, p_eigs in enumerate( conf.options.check_eigs ):
            if isinstance( p_eigs, tuple ):
                n_eigs = [pn( p_eigs[0], n_eigs_total ), pn( p_eigs[1], n_eigs_total )]
                if (sum( p_eigs ) == 100) and (sum( n_eigs ) < n_eigs_total):
                    n_eigs[1] = n_eigs_total - n_eigs[0]
            else:
                n_eigs = pn( p_eigs, n_eigs_total ).tolist()

            conf.options.n_eigs = n_eigs
#            conf.options.n_eigs = n_eigs
            print 'eig test %d: %s' % (ii, n_eigs)
#            pause()
            sol_t, coefs_t = solve_time_variant_correctors( sol_s, problem,
                                                       conf, options,
                                                       ofn_trunk + ('_%d' % ii) )
            coefs = coefs_s + coefs_t
            if isinstance( n_eigs, tuple ):
                fname = 'coefs_%d_%d_%03.2f_%03.2f.h5' % (n_eigs + p_eigs)
            else:
                fname = 'coefs_%d_%03.2f.h5' % (n_eigs, p_eigs)
            coefs.to_file_hdf5( fname )
    else:
        sol_t, coefs_t = solve_time_variant_correctors( sol_s, problem, conf, options,
                                                   ofn_trunk )

    return sol_s + sol_t, coefs_s + coefs_t

def get_evp( key, cache_evp, problem, conf, equivalence = None ):
    """
    Calls pfdpm.solve_pressure_eigenproblem() if the same or equivalent EVP
    was not already solved - in that case returns the cached EVP.
    """
    evp = None
    if equivalence is None:
        if cache_evp.has_key( key ):
            evp = cache_evp[key]
    else:
        for key2 in equivalence[key]:
            if key2 in cache_evp:
                evp = cache_evp[key2]
                cache_evp[key] = evp
                break
                
    if evp is None:
        ebc = pfdpm.select_by_names( conf.ebcs, conf.ebc_sets[key] )
        epbc = pfdpm.select_by_names( conf.epbcs, conf.epbc_sets[key] )

        solve = pfdpm.get_matrix_parts
        matrices = solve( problem, ebc, epbc, conf.variables,
                          conf.equations_time )
        solve = pfdpm.solve_pressure_eigenproblem
        evp = solve( matrices, conf.options.eig_problem,
                     conf.options.n_eigs,
                     conf.options.check.get( 'diagonalization', False ) )
        evp.ebc, evp.epbc = ebc, epbc
        cache_evp[key] = evp
    else:
        # Just create/restore equation mappings.
        variables = select_by_names( conf.variables, ['uc', 'vc', 'pc', 'qc'] )
        problem.set_variables( variables )
        problem.set_equations( conf.equations_time )
        problem.time_update( conf_ebc = evp.ebc, conf_epbc = evp.epbc )

    return evp

def build_evp_equivalence( equivs_in ):
    if equivs_in is None:
        return None

    out = {}
    for equiv in equivs_in:
        for item in equiv:
            if item in out:
                out[item].add( equiv )
            else:
                out[item] = set( equiv )
    return out

def compute_micro_cefficients( conf, options, ret_all = False ):
    """Dependencies must be listed in a correct order."""
    aux = get_default_attr( conf, 'evp_equivalence', None )
    evp_equivalence = build_evp_equivalence( aux )

    opts = conf.options
    if hasattr( opts, 'post_process_hook' ) and opts.post_process_hook is not None:
        # User postprocessing.
        post_process_hook = getattr( conf.funmod, opts.post_process_hook )
    else:
        post_process_hook = None
    file_per_var = get_default_attr( opts, 'file_per_var', True )

    problem = ProblemDefinition.from_conf( conf,
                                           init_variables = False,
                                           init_equations = False )

    dependencies = {}
    cache_evp = {}
    
    coefs = pfdpm.Coefficients()
    coefs.filename = conf._filename
    
    coef_info = getattr( conf, opts.coef_info )
    for coef_name, cargs in coef_info.iteritems():
        output( 'computing %s...' % coef_name )
        requires = cargs.get( 'requires', [] )
        for req in requires:
            if dependencies.has_key( req ) and (dependencies[req] is not None):
                continue

            output( 'computing dependency %s...' % req )

            if req == 'pis':
                dependencies['pis'] = pfdpm.create_pis( problem,
                                                       conf.variables, 'uc' )

            elif req == 'corrs_rs':
                ##
                # Solve steady_rs.
                solve = pfdpm.solve_steady_correctors_rs
                aux = solve( problem, conf.ebcs, conf.epbcs, conf.variables,
                             conf.equations_steady_rs,
                             conf.ebc_sets['corrs_rs'], conf.epbc_sets['corrs_rs'],
                             dependencies['pis'], opts.file_conf[req],
                             post_process_hook, file_per_var )
                dependencies['corrs_rs'] = aux

            elif req == 'corrs_time_rs':
                evp = get_evp( req, cache_evp, problem, conf,
                              equivalence = evp_equivalence )

                dim = problem.domain.mesh.dim
                ts = problem.get_time_solver().ts
                corrs_rs = dependencies['corrs_rs']

                output( 'time-variant rs correctors via eigensolutions...' )
                filenames_rs = nm.zeros( (dim, dim), dtype = nm.object )
                solve = pfdpm.make_t_correctors_via_evp
                for ir in range( dim ):
                    for ic in range( dim ):
                        filename =( opts.file_conf[req] % (ir,ic)) + '.h5'
                        filenames_rs[ir,ic] = filename
                        solve( problem, evp, -corrs_rs.states_rs[ir,ic], ts,
                               filename, post_process_hook = post_process_hook,
                               file_per_var = file_per_var )
                output( '...done' )

                if opts.check.get( 'time_correctors', False ):
                    output( 'verifying correctors %s...' % req )
                    verify = pfdpm.verify_t_correctors
                    ok = True
                    for ir in range( dim ):
                        for ic in range( dim ):
                            oo = verify( problem, evp.ebc, evp.epbc,
                                         conf.equations_time_debug,
                                         -corrs_rs.states_rs[ir,ic],
                                         filenames_rs[ir,ic] )
                            ok = ok and oo
                    output( '...done, ok: %s' % ok )


                dependencies['corrs_time_rs'] = filenames_rs

            elif req in ['corrs_alpha1', 'corrs_alpha2']:
                ##
                # Solve steady_alpha.
                alpha = int( req[-1] )
                solve = pfdpm.solve_steady_correctors_alpha
                aux = solve( problem, conf.ebcs, conf.epbcs, conf.variables,
                             conf.equations_steady_alpha,
                             conf.ebc_sets[req],
                             conf.epbc_sets[req],
                             opts.file_conf[req],
                             alpha, post_process_hook, file_per_var )
                dependencies[req] = aux

            elif req in ['corrs_time_alpha1', 'corrs_time_alpha2']:
                ts = problem.get_time_solver().ts
                alpha = int( req[-1] )
                corrs_alpha = dependencies['corrs_alpha%d' % alpha]

                output( 'time-variant alpha correctors via eigensolutions...' )
                solve = pfdpm.make_t_correctors_via_evp

                evp = get_evp( req, cache_evp, problem, conf,
                              equivalence = evp_equivalence )
                filename = (opts.file_conf[req] % alpha) + '.h5'
                solve( problem, evp, corrs_alpha.state, ts, filename,
                       alpha = alpha, post_process_hook = post_process_hook,
                       file_per_var = file_per_var )
                output( '...done' )

                if opts.check.get( 'time_correctors', False ):
                    output( 'verifying correctors %s...' % req )
                    verify = pfdpm.verify_t_correctors
                    ok = verify( problem, evp.ebc, evp.epbc,
                                 conf.equations_time_debug,
                                 corrs_alpha.state,
                                 filename )
                    output( '...done, ok: %s' % ok )
                dependencies[req] = filename

            elif req == 'corrs_pressure':
                ##
                # Solve steady_pressure.
                solve = pfdpm.solve_steady_correctors_pressure
                aux = solve( problem, conf.ebcs, conf.epbcs, conf.variables,
                             conf.equations_steady_pressure,
                             conf.ebc_sets['corrs_pressure'],
                             conf.epbc_sets['corrs_pressure'],
                             opts.file_conf[req], post_process_hook, file_per_var )
                dependencies['corrs_pressure'] = aux

            elif req == 'corrs_time_pressure':
                evp = get_evp( req, cache_evp, problem, conf,
                              equivalence = evp_equivalence )

                ts = problem.get_time_solver().ts
                corrs_pressure = dependencies['corrs_pressure']

                output( 'time-variant pressure correctors via eigensolutions...' )
                solve = pfdpm.make_t_correctors_via_evp
                filename = opts.file_conf[req] + '.h5'
                solve( problem, evp, corrs_pressure.state, ts,
                       filename, post_process_hook = post_process_hook,
                       file_per_var = file_per_var )

                output( '...done' )

                if opts.check.get( 'time_correctors', False ):
                    output( 'verifying correctors %s...' % req )
                    ok = pfdpm.verify_t_correctors( problem, evp.ebc, evp.epbc,
                                                  conf.equations_time_debug,
                                                  corrs_pressure.state,
                                                  filename )
                    output( '...done, ok: %s' % ok )

                dependencies['corrs_time_pressure'] = filename

            else:
                print 'unknown dependency: %s' % req
                raise ValueError

            output( '...done' )
                
        if coef_name == 'times':
            val = problem.get_time_solver().ts.times

        elif coef_name in ['C1', 'C2', 'C']:
            ##
            # Permeabilities related to channels.
            val = pfdpm.solve_permeability( problem,
                                           conf.ebcs, conf.epbcs,
                                           conf.variables,
                                           conf.equations_coefs, cargs )

        elif coef_name == 'VF':
            ##
            # Volume fractions based on 'uc' geometry.
            aux, volume_all = pfdpm.volumes( problem, conf.variables,
                                             conf.equations_coefs, cargs )
            vf = {}
            vv = 0.0
            for key, val in aux.iteritems():
                setattr( coefs, 'volume%s' % key, val )
                if key != volume_all:
                    vf[key] =  nm.array( val / aux[volume_all],
                                         dtype = nm.float64 )
                    vv += val
            assert_( abs( vv - coefs.volumeY ) < 1e-14 )
            val = vf

        elif coef_name == 'E':
            val = pfdpm.coef_e( problem,
                               dependencies['corrs_rs'],
                               dependencies['pis'],
                               conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['GStar1', 'GStar2']:
            alpha = req[-1]
            val = pfdpm.coef_g_star( problem, dependencies['corrs_alpha' + alpha],
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['GBar1', 'GBar2']:
            val = pfdpm.coef_g_bar( problem,
                                  conf.ebcs, conf.epbcs,
                                  conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['GPlus1', 'GPlus2']:
            alpha = req[-1]
            val = pfdpm.coef_g_plus( problem,
                                   dependencies['corrs_time_alpha' + alpha],
                                   conf.ebcs, conf.epbcs,
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotE':
            val = pfdpm.coef_biot_e( problem,
                                   dependencies['corrs_rs'],
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotMIR':
            val = pfdpm.coef_biot_mir( problem,
                                     dependencies['corrs_pressure'],
                                     conf.variables, conf.equations_coefs,
                                     cargs )
        elif coef_name == 'H':
            val = pfdpm.coef_h( problem,
                               dependencies['corrs_rs'],
                               dependencies['corrs_time_rs'],
                               conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotH':
            val = pfdpm.coef_biot_h( problem,
                                   dependencies['corrs_time_rs'],
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotH2':
            val = pfdpm.coef_biot_h2( problem,
                                    dependencies['corrs_time_rs'],
                                    dependencies['corrs_time_pressure'],
                                    dependencies['corrs_rs'],
                                    dependencies['corrs_pressure'],
                                    conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotMHR':
            val = pfdpm.coef_biot_mhr( problem,
                                     dependencies['corrs_time_pressure'],
                                     conf.variables, conf.equations_coefs,
                                     cargs )

        elif coef_name in ['P1', 'P2']:
            alpha = req[-1]
            val = pfdpm.coef_p( problem,
                               dependencies['corrs_alpha' + alpha],
                               dependencies['pis'],
                               conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['R1', 'R2']:
            alpha = req[-1]
            val = pfdpm.coef_r( problem,
                               dependencies['corrs_time_alpha' + alpha],
                               dependencies['pis'],
                               conf.variables, conf.equations_coefs, cargs )

        else:
            print 'unknown coefficient: %s' % coef_name
            raise ValueError

        setattr( coefs, coef_name, val )
        output( '...done' )

    prec = nm.get_printoptions()[ 'precision']
    if hasattr( opts, 'print_digits' ):
        nm.set_printoptions( precision = opts.print_digits )
    print coefs
    nm.set_printoptions( precision = prec )
##     pause()

    if ret_all:
        return coefs, dependencies
    else:
        return coefs

##
# c: 13.06.2008, r: 22.06.2008
def verify_steady_solution( conf, options ):
    if not hasattr( conf, 'equations_steady' ):
        output( 'set "equations_steady" in the input!' )
        return False

    ok = True
    pb = ProblemDefinition.from_conf( conf )
    opts = conf.options
    if hasattr( opts, 'post_process_hook' ) and opts.post_process_hook is not None:
        # User postprocessing.
        pph = getattr( conf.funmod, opts.post_process_hook )
    else:
        pph = None
    state_t, data_t = solve_evolutionary_op( pb, options, post_process_hook = pph )

    pb.set_equations( conf.equations_steady )
    tsc = pb.ts_conf
    ts = pb.get_default_ts( tsc.t0, tsc.t1, tsc.dt, tsc.n_step, tsc.n_step - 1 )
    # To reassemble matrix with new equations.
    pb.set_linear( False )
    state_s, data_s = solve_stationary_op( pb, options, ts = ts,
                                       post_process_hook = pph )

    err = nla.norm( state_s - state_t )
    eps = get_default_attr( conf.options, 'steady_state_error', 1e-8 )
    if err > eps:
        ok = False

    output( 'error: %.2e' % err )

    aux = eval_term_op( state_s, 'dw_volume_wdot.i1.Omega( m.BiotMHI, q, p )',
                      pb, dw_mode = 'vector' )
    print nla.norm( aux )
    aux = eval_term_op( state_t, 'dw_volume_wdot.i1.Omega( m.BiotMHI, q, p )',
                      pb, dw_mode = 'vector' )
    print nla.norm( aux )

    import pylab
    pylab.plot( state_s )
    pylab.plot( state_t )
    pylab.show()

    return ok

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'steady' :
    'steady state solution',
    'impact' :
    'impact loading solution',
    'micro' :
    'microproblem solution',
    'debug_micro' :
    'debug microproblem solution (requires -m) [default: %default]',
    'variant' :
    'Y_1 or Y_2 variant of alpha correctors (requires -m) [default: both]',
    'recovery' :
    'microproblem solution recovery using macroscopic problem solution'
    ' [default: %default]',
    'new_micro' :
    'new microproblem solution',
    'verify_steady' :
    'verify steady state',
}

##
# c: 11.07.2006, r: 18.06.2008
def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "-s", "--steady",
                       action = "store_true", dest = "steady",
                       default = False, help = help['steady'] )
    parser.add_option( "-i", "--impact",
                       action = "store_true", dest = "impact",
                       default = False, help = help['impact'] )
    parser.add_option( "-m", "--micro",
                       action = "store_true", dest = "micro",
                       default = False, help = help['micro'] )
    parser.add_option( "-r", "--recovery",
                       action = "store_true", dest = "recovery",
                       default = False, help = help['recovery'] )
    parser.add_option( "", "--debug-micro",
                       action = "store_true", dest = "debug_micro",
                       default = False, help = help['debug_micro'] )
    parser.add_option( "", "--variant", type = "int", metavar = '1 or 2',
                       action = "store", dest = "variant",
                       default = 0, help = help['variant'] )
    parser.add_option( "-n", "--new-micro",
                       action = "store_true", dest = "new_micro",
                       default = False, help = help['new_micro'] )
    parser.add_option( "-v", "--verify-steady",
                       action = "store_true", dest = "verify_steady",
                       default = False, help = help['verify_steady'] )
    (options, args) = parser.parse_args()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    if options.impact: options.steady = True

    if ((options.debug_micro and not options.micro)
        or (options.variant and not options.micro)) :
        output( 'assuming microproblem solution' )
        options.micro = True
        options.steady = options.impact = False
        

    required, other = get_standard_keywords()

    if options.steady:
        required.remove( 'save_steps' )
        required.remove( 'equations' )
        required += ['equations_steady_u', 'equations_steady_p',
                     'steady_state_n_time']
    elif options.micro:
        required.remove( 'save_steps' )
        required.remove( 'equations' )
        other.remove( 'epbc' )
        other.remove( 'options' )
        required += ['epbc', 'equations_steady_rs', 'equations_steady_alpha',
                     'equations_time_variant', 'options']
    if options.recovery:
        required.remove( 'save_steps' )
        required += ['macro_solution', 'micro_correctors']

    if options.impact:
        required += ['equations_impact']

    if options.new_micro:
        required.remove( 'equations' )

    conf = ProblemConf.from_file( filename_in, required, other )
##     print conf
##     pause()

    if options.steady:
        out = solve_steady( conf, options )
        if options.impact:
            solve_impact( conf, options, *out )
    elif options.micro:
        micro_solution, coefs = solve_micro( conf, options )
        coefs.to_file_hdf5( 'coefs.h5' )
        print coefs
        pause()
        c2 = pfdpm.Coefficients.from_file_hdf5( 'coefs.h5' )
        print c2
    elif options.recovery:
        macro_solution = conf.macro_solution
        micro_correctors = conf.micro_correctors
        solve_recovery( conf, options, macro_solution, micro_correctors )
    elif options.new_micro:
        coefs = compute_micro_cefficients( conf, options )

        coefs.to_file_hdf5( 'coefs.h5' )
        coefs.to_file_txt( 'coefs.txt',
                           conf.options.tex_names, conf.options.float_format )
    elif options.verify_steady:
        ok = verify_steady_solution( conf, options )
        if not ok:
            output( 'failed!' )
        else:
            output( 'ok!' )
    else:
        solve_evolutionary( conf, options )

    

if __name__ == '__main__':
    main()

    ##
    # Test state vector.
##     aux = nm.arange( len( state ) / 2, dtype = nm.float64 )
##     state[0::2] = aux
##     state[1::2] = -aux

##     state = nm.arange( len( state ), dtype = nm.float64 )
##     print state
##     pause()

##         vec_u0 = problem.variables.get_state_part_view( state0, 'u' )
##         print vec_u0[var.eq_map.eq_ebc]
##         vec_u = problem.variables.get_state_part_view( state, 'u' )
##         print vec_u[var.eq_map.eq_ebc]
##         pause()

##             out2 = read_data_hdf5( 'pokus.h5', step )
##             for key in out.keys():
##                 print out[key]
##                 print out2[key]
##                 pause()

def debug_dw_surface_ltr( problem, conf ):
##     problem.update_bc( ts, {}, {}, conf.funmod )
    problem.materials.time_update( ts, conf.funmod, problem.domain )
    val = eval_term_op( None, 'dw_surface_ltr.Gamma2( one, vc )', problem )
    trac = problem.variables.get_state_part_view( val, 'uc' )
    print trac.shape
    print problem.variables.di
    print nm.sum( trac )
    print len( nm.where( trac != 0.0 )[0] )
    tr = nm.c_[trac[0::2], trac[1::2]]
    print tr
    from matplotlib.mlab import save
    save( 'tr2', tr, fmt=' % .7e' )
    import pylab
    pylab.plot( trac[0::2] )
    pylab.plot( trac[1::2] )
    pylab.show()
    pause()
