from sfepy.base.base import *

import sfepy.base.ioutils as io
from sfepy.fem import ProblemDefinition
from sfepy.base.conf import get_standard_keywords

##
# c: 03.07.2007, r: 27.02.2008
def save_only( conf, save_names, problem = None ):
    """Save information available prior to setting equations and
    solving them."""
    if problem is None:
        problem = ProblemDefinition.from_conf( conf, init_variables = False )

    if save_names.regions is not None:
        problem.save_regions( save_names.regions )

    if save_names.field_meshes is not None:
        problem.save_field_meshes( save_names.field_meshes )

    if save_names.region_field_meshes is not None:
        problem.save_region_field_meshes( save_names.region_field_meshes )

    if save_names.ebc is not None:
        if not hasattr( problem, 'variables' ):
            problem.set_variables( conf.variables )
        try:
            ts = TimeStepper.from_conf( conf.ts )
            ts.set_step( 0 )
        except:
            ts = None
        try:
            problem.variables.equation_mapping( conf.ebcs, conf.epbcs,
                                                problem.domain.regions, ts,
                                                conf.funmod )
        except Exception, e:
            output( 'cannot make equation mapping!' )
            output( 'reason: %s' % e )
        else:
            problem.save_ebc( save_names.ebc )

##
# 20.03.2007, c
# 30.03.2007
# 28.05.2007
# 03.07.2007
# 18.07.2007
# 02.10.2007
# 03.10.2007
def solve_stationary( conf, data = None, save_names = None, nls_status = None ):

    if data is None:
        # Term-dependent data.
        data = {}
    problem = ProblemDefinition.from_conf( conf )

    problem.time_update( None )

    if save_names is not None:
        save_only( conf, save_names, problem = problem )

    state = problem.solve( nls_status = nls_status )

    return problem, state, data


def prepare_save_data( ts, conf ):
    try:
        save_steps = conf.options.save_steps
    except:
        save_steps = -1
    if save_steps == -1:
        save_steps = ts.n_step

    is_save = nm.linspace( 0, ts.n_step - 1, save_steps ).astype( nm.int32 )
    is_save = nm.unique1d( is_save )

    return ts.suffix, is_save

def time_step_function( ts, state0, problem, data ):
    problem.time_update( ts )

    if ts.step == 0:
        problem.apply_ebc( state0 )
        state = state0.copy()
        problem.init_time( ts )

        if problem.equations.caches:
            # Initialize caches.
            ev = problem.get_evaluator( ts = ts, **data )
            try:
                vec_r = ev.eval_residual( state, is_full = True )
            except ValueError:
                output( 'initial residual evaluation failed, giving up...' )
                raise
            else:
                err = nla.norm( vec_r )
                output( 'initial residual: %e' % err )

        else:
            # Just initialize data of state variables.
            problem.variables.data_from_state( state )

        if problem.is_linear():
            # Assemble linear system matrix for all
            # time steps.
            ev = problem.get_evaluator( ts = ts, mtx = problem.mtx_a, **data )
            try:
                mtx_a = ev.eval_tangent_matrix( state, is_full = True )
            except ValueError:
                output( 'matrix evaluation failed, giving up...' )
                raise
        else:
            mtx_a = None

        # Initialize solvers (and possibly presolve the matrix).
        problem.init_solvers( ts = ts, mtx = mtx_a, **data )
        # Initialize variables with history.
        problem.init_variables( state0 )

    else:
        state = problem.solve( state0 = state0, ts = ts, **data )

    return state

def solve_evolutionary_op( problem,
                           save_results = True, return_history = False,
                           post_process_hook = None, step_hook = None ):
    """TODO  return_history"""
    
    data = {}
    time_solver = problem.get_time_solver( step_fun = time_step_function,
                                           step_args = (problem, data) )

    suffix, is_save = prepare_save_data( time_solver.ts,
					 problem.conf )

    state0 = problem.create_state_vector()
    problem.setup_ic()
    problem.apply_ic( state0 )
    
    ii = 0
    for ts, state in time_solver( state0 ):

        if step_hook is not None:
            step_hook( problem, ts, state )

        if save_results and (is_save[ii] == ts.step):
            filename = problem.get_output_name( suffix = suffix % ts.step )
            problem.save_state( filename, state,
                                post_process_hook = post_process_hook )
            ii += 1

        problem.advance( ts )
    return state, data

def solve_stationary_op( problem, save_results = True, ts = None,
                         post_process_hook = None ):
    data = {}
    problem.time_update( ts )
    state = problem.solve()

    if save_results:
        problem.save_state( problem.get_output_name(), state,
                            post_process_hook = post_process_hook )

    return state, data
    
def solve_direct( conf, options, problem = None, step_hook = None,
                  post_process_hook = None ):
    """Generic (simple) problem solver."""
    if problem is None:
	problem = ProblemDefinition.from_conf( conf )
	if options.output_filename_trunk:
	    ofn_trunk = options.output_filename_trunk
	    problem.ofn_trunk = ofn_trunk
	if options.output_format:
	    problem.output_format = options.output_format
    ofn_trunk = problem.ofn_trunk
    
    save_names = Struct( ebc = None, regions = None, field_meshes = None,
                         region_field_meshes = None )
    if options.save_ebc:
        save_names.ebc = ofn_trunk + '_ebc.vtk'
    if options.save_regions:
        save_names.regions = ofn_trunk + '_region'
    if options.save_field_meshes:
        save_names.field_meshes = ofn_trunk + '_field'
    if options.save_region_field_meshes:
        save_names.region_field_meshes = ofn_trunk + '_region_field'

    is_extra_save = False
    for name, val in save_names.to_dict().iteritems():
        if val is not None:
            is_extra_save = True
            break
    if is_extra_save:
        save_only( conf, save_names )

    if options.solve_not:
        return None, None, None
            
    if hasattr( conf.options, 'ts' ):
        ##
        # Time-dependent problem.
        out = solve_evolutionary_op( problem, options,
				     post_process_hook = post_process_hook,
                                     step_hook = step_hook )
    else:
        ##
        # Stationary problem.
        out = solve_stationary_op( problem, options,
				   post_process_hook = post_process_hook )

    state, data = out
    return problem, state, data
