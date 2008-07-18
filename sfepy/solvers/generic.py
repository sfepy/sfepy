from sfepy.base.base import *

import sfepy.base.ioutils as io
from sfepy.fem.problemDef import ProblemDefinition
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


##
# c: 06.02.2008, r: 09.07.2008
def prepare_save_data( ts, conf, options ):
    if options.output_file_name_trunk:
        ofn_trunk = options.output_file_name_trunk
    else:
        ofn_trunk = io.get_trunk( conf.file_name_mesh ) + '_out'

    try:
        save_steps = conf.options.save_steps
    except:
        save_steps = -1
    if save_steps == -1:
        save_steps = ts.n_step

    is_save = nm.linspace( 0, ts.n_step - 1, save_steps ).astype( nm.int32 )
    is_save = nm.unique1d( is_save )

    return ofn_trunk, ts.suffix, is_save

def time_step_function( ts, state0, problem, data ):
    problem.time_update( ts )

    if ts.step == 0:
        problem.apply_ebc( state0 )
        state = state0.copy()
        problem.init_time( ts )

        if problem.equations.caches:
            # Initialize caches.
            ev = problem.get_evaluator( ts = ts, **data )
            vec_r, ret = ev.eval_residual( state )
            if ret == 0: # OK.
                err = nla.norm( vec_r )
                output( 'initial residual: %e' % err )
            else:
                output( 'initial residual evaluation failed, giving up...' )
                raise ValueError

        if problem.is_linear():
            # Assemble linear system matrix for all
            # time steps.
            ev = problem.get_evaluator( ts = ts, mtx = problem.mtx_a, **data )
            mtx_a, ret = ev.eval_tangent_matrix( state )
            if ret != 0:
                output( 'matrix evaluation failed, giving up...' )
                raise ValueError
        else:
            mtx_a = None

        # Initialize solvers (and possibly presolve the matrix).
        problem.init_solvers( ts = ts, mtx = mtx_a, **data )
        # Initialize variables with history.
        problem.init_variables( state0 )

    else:
        state = problem.solve( state0 = state0, ts = ts, **data )

    problem.advance( ts )

    return state

def solve_evolutionary_op( problem, options,
                         save_results = True, return_history = False,
                         post_process_hook = None ):
    """TODO  return_history"""
    
    data = {}
    time_solver = problem.get_time_solver( step_fun = time_step_function,
                                        step_args = (problem, data) )

    ofn_trunk, suffix, is_save = prepare_save_data( time_solver.ts,
                                                problem.conf, options )

    state0 = problem.create_state_vector()
    ii = 0
    for step, time, state in time_solver( state0 ):

        if save_results and (is_save[ii] == step):
            problem.save_state( ofn_trunk + suffix % step, state,
                               post_process_hook = post_process_hook )

            ii += 1
    return state, data

##
# c: 13.06.2008, r: 13.06.2008
def solve_stationary_op( problem, options, save_results = True, ts = None,
                       post_process_hook = None ):
    data = {}
    problem.time_update( ts )
    state = problem.solve()

    if save_results:
        if options.output_file_name_trunk:
            ofn_trunk = options.output_file_name_trunk
        else:
            ofn_trunk = io.get_trunk( problem.conf.file_name_mesh ) + '_out'
        problem.save_state( ofn_trunk + '.vtk', state,
                           post_process_hook = post_process_hook )

    return state, data
    
##
# c: 12.01.2007, r: 22.06.2008
def solve_direct( conf, options ):
    """Generic (simple) problem solver."""
    if options.output_file_name_trunk:
        ofn_trunk = options.output_file_name_trunk
    else:
        ofn_trunk = io.get_trunk( conf.file_name_mesh )

    opts = conf.options
    if hasattr( opts, 'post_process_hook' ) and opts.post_process_hook is not None:
        # User postprocessing.
        post_process_hook = getattr( conf.funmod, opts.post_process_hook )
    else:
        post_process_hook = None

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
            
    pb = ProblemDefinition.from_conf( conf )
    if hasattr( conf.options, 'ts' ):
        ##
        # Time-dependent problem.
        state, data = solve_evolutionary_op( pb, options,
                                           post_process_hook = post_process_hook )
    else:
        ##
        # Stationary problem.
        state, data = solve_stationary_op( pb, options,
                                         post_process_hook = post_process_hook )

    return pb, state, data
