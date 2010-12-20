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
        problem = ProblemDefinition.from_conf(conf, init_equations=False)

    if save_names.regions is not None:
        problem.save_regions( save_names.regions )

    if save_names.regions_as_groups is not None:
        problem.save_regions_as_groups( save_names.regions_as_groups )

    if save_names.field_meshes is not None:
        problem.save_field_meshes( save_names.field_meshes )

    if save_names.region_field_meshes is not None:
        problem.save_region_field_meshes( save_names.region_field_meshes )

    if save_names.ebc is not None:
        problem.save_ebc( save_names.ebc )

def solve_stationary(conf, save_names=None, nls_status=None):

    problem = ProblemDefinition.from_conf( conf )

    problem.time_update( None )

    if save_names is not None:
        save_only( conf, save_names, problem = problem )

    state = problem.solve( nls_status = nls_status )

    return problem, state


def prepare_save_data( ts, conf ):
    try:
        save_steps = conf.options.save_steps
    except:
        save_steps = -1
    if save_steps == -1:
        save_steps = ts.n_step

    is_save = nm.linspace( 0, ts.n_step - 1, save_steps ).astype( nm.int32 )
    is_save = nm.unique( is_save )

    return ts.suffix, is_save

def prepare_matrix(problem, state):
    """
    Pre-assemble tangent system matrix.
    """
    problem.update_materials()

    ev = problem.get_evaluator()
    try:
        mtx = ev.eval_tangent_matrix(state(), is_full=True)

    except ValueError:
        output('matrix evaluation failed, giving up...')
        raise

    return mtx

def time_step_function(ts, state0, problem, nls_status=None):
    problem.time_update( ts )

    if ts.step == 0:
        state0.apply_ebc()
        state = state0.copy(deep=True)

        if not ts.is_quasistatic:
            problem.init_time( ts )

            if problem.equations.caches:
                # Initialize caches.
                ev = problem.get_evaluator()
                try:
                    vec_r = ev.eval_residual(state(), is_full=True)
                except ValueError:
                    output( 'initial residual evaluation failed, giving up...' )
                    raise
                else:
                    err = nla.norm( vec_r )
                    output( 'initial residual: %e' % err )

        if problem.is_linear():
            mtx = prepare_matrix(problem, state)

        else:
            mtx = None

        # Initialize solvers (and possibly presolve the matrix).
        problem.init_solvers(nls_status=nls_status, mtx=mtx)

        if ts.is_quasistatic:
            # Ordinary solve.
            state = problem.solve(state0=state0)
            state.init_history()

        else:
            # Initialize variables with history.
            state0.init_history()

    else:
        if (ts.step == 1) and ts.is_quasistatic and problem.is_linear():
            mtx = prepare_matrix(problem, state0)
            problem.init_solvers(nls_status=nls_status, mtx=mtx)

        state = problem.solve(state0=state0)

    return state

def solve_evolutionary_op(problem,
                          save_results=True, return_history=False,
                          step_hook=None, post_process_hook=None,
                          nls_status=None):
    """TODO  return_history"""
    
    step_args = problem, nls_status
    time_solver = problem.get_time_solver(step_fun=time_step_function,
                                          step_args=step_args)

    suffix, is_save = prepare_save_data( time_solver.ts,
                                         problem.conf )

    state0 = problem.create_state()
    problem.setup_ic()
    state0.apply_ic()

    ii = 0
    for ts, state in time_solver( state0 ):

        if step_hook is not None:
            step_hook( problem, ts, state )

        if save_results and (is_save[ii] == ts.step):
            filename = problem.get_output_name( suffix = suffix % ts.step )
            problem.save_state(filename, state,
                               post_process_hook=post_process_hook,
                               file_per_var=None,
                               ts=ts)
            ii += 1

        problem.advance( ts )
    return state

def solve_stationary_op(problem, save_results=True, ts=None,
                        post_process_hook=None,
                        nls_status=None):
    if ts is None:
        try:
            ts = problem.get_time_solver().ts
        except ValueError:
            pass

    problem.time_update( ts )
    state = problem.solve(nls_status=nls_status)

    if save_results:
        problem.save_state(problem.get_output_name(), state,
                           post_process_hook=post_process_hook,
                           file_per_var=None)

    return state
    
def solve_direct(conf, options, problem=None, step_hook=None,
                 post_process_hook=None, post_process_hook_final=None,
                 nls_status=None):
    """Generic (simple) problem solver."""
    if problem is None:
        is_eqs = not options.solve_not
        problem = ProblemDefinition.from_conf(conf, init_equations=is_eqs)

        problem.setup_default_output(conf, options)

    ofn_trunk = problem.ofn_trunk

    save_names = Struct( ebc = None, regions = None,
                         regions_as_groups = None, field_meshes = None,
                         region_field_meshes = None )
    if options.save_ebc:
        save_names.ebc = ofn_trunk + '_ebc.vtk'
    if options.save_regions:
        save_names.regions = ofn_trunk + '_region'
    if options.save_regions_as_groups:
        save_names.regions_as_groups = ofn_trunk + '_regions'
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
        save_only( conf, save_names, problem=problem )

    if options.solve_not:
        return None, None, None
            
    if hasattr( conf.options, 'ts' ):
        ##
        # Time-dependent problem.
        state = solve_evolutionary_op(problem, options,
                                      step_hook=step_hook,
                                      post_process_hook=post_process_hook,
                                      nls_status=nls_status)
    else:
        ##
        # Stationary problem.
        state = solve_stationary_op(problem, options,
                                    post_process_hook=post_process_hook,
                                    nls_status=nls_status)

    if post_process_hook_final is not None: # User postprocessing.
       post_process_hook_final(problem, state)

    return problem, state
