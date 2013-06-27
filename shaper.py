#!/usr/bin/env python
# 06.04.2005, c
# 16.06.2005
from optparse import OptionParser

import numpy as nm

import sfepy
from sfepy.base.base import output, remap_dict, Struct, IndexedStruct
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.fem.evaluate import BasicEvaluator
import sfepy.base.ioutils as io
import sfepy.optimize.shapeOptim as so
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.solvers import Solver

def solve_stokes(dpb, equations_stokes, nls_conf):
    dpb.set_equations(equations_stokes)
    dpb.time_update(None)

    output('solving Stokes problem...')
    state = dpb.solve(nls_conf=nls_conf)
    output('...done')

    return state

def solve_navier_stokes(conf, options):
    opts = conf.options

    dpb = ProblemDefinition.from_conf(conf, init_equations=False)
    equations = getattr(conf, '_'.join(('equations_direct', opts.problem)))
    dpb.set_equations(equations)

    ls_conf = dpb.get_solver_conf( opts.ls )
    nls_conf = dpb.get_solver_conf(opts.nls_direct)

    method = opts.direct_method
    if method == 'stationary':
        data = {}
        dpb.time_update(None)
        state_dp = dpb.solve(nls_conf=nls_conf)

    elif method == 'transient':
        ls = Solver.any_from_conf( ls_conf )
        ts_conf = dpb.get_solver_conf( opts.ts_direct )

        data = {'ts' : Struct( dt = ts_conf.dt )}

        # Plug in mass term.
        mequations = {}
        for key, eq in equations.iteritems():
            if 'dw_div_grad' in eq:
                eq = '+'.join( (ts_conf.mass_term, eq) ).replace( '++', '+')
            mequations[key] = eq

        if ts_conf.stokes_init:
            state_dp0 = solve_stokes( dpb, conf.equations_direct_stokes, nls_conf )
            dpb.set_equations( mequations )
        else:
            dpb.set_equations( mequations )
            state_dp0 = dpb.create_state()
            dpb.time_update( None )
            state_dp0.apply_ebc()

        from sfepy.base.log import Log

        log = Log.from_conf( Struct( is_plot = True ),
                            ([r'$||u||$'], [r'$||p||$']) )

        output( 'Navier-Stokes...' )
        ev = BasicEvaluator( dpb, ts = Struct( dt = ts_conf.dt ) )
        nls = Solver.any_from_conf( nls_conf, evaluator = ev, lin_solver = ls )

        n_step = ts_conf.n_step
        step = 0
        while 1:
            for ii in xrange( n_step ):
                output( step )

                vec_u = state_dp0('w')
                vec_p = state_dp0('r')
                log( nm.linalg.norm( vec_u ), nm.linalg.norm( vec_p ) )

                dpb.variables.set_data_from_state('w_0', state_dp0(), 'w')
                vec_dp = nls( state_dp0() )

                step += 1
                state_dp = state_dp0.copy()
                state_dp.set_reduced(vec_dp)

                state_dp0 = state_dp

            if ts_conf.interactive:
                try:
                    n_step = int( raw_input( 'continue: ' ) )
                    if n_step <= 0: break
                except:
                    break

        vec_u = state_dp('w')
        vec_p = state_dp('r')
        log( nm.linalg.norm( vec_u ), nm.linalg.norm( vec_p ), finished = True )

    else:
        raise 'unknown Navier-Stokes solution method (%s)!'  % method

    return dpb, state_dp, data

def solve_generic_direct(conf, options):
    opts = conf.options

    dpb = ProblemDefinition.from_conf(conf, init_equations=False)
    equations = getattr(conf, '_'.join(('equations_direct', opts.problem)))
    dpb.set_equations(equations)

    dpb.time_update(None)

    nls_conf = dpb.get_solver_conf(opts.nls_direct)
    state_dp = dpb.solve(nls_conf=nls_conf)

    return dpb, state_dp, {}

##
# c: 22.11.2006, r: 15.04.2008
def solve_direct( conf, options ):
    """
    Solve the direct (nonlinear) problem.
    """
    opts = conf.options
    if hasattr( opts, 'problem' ):
        if opts.problem == 'navier_stokes':
            dpb, state_dp, data = solve_navier_stokes( conf, options )
        else:
            output( 'unknown problem type (%s), using generic solver.'\
                    % opts.problem )
            dpb, state_dp, data = solve_generic_direct( conf, options )
    else: # Generic direct problem.
        dpb, state_dp, data = solve_generic_direct( conf, options )

    trunk = io.get_trunk( conf.filename_mesh )
    dpb.save_state( trunk + '_direct.vtk', state_dp )

    if options.dump_filename is not None:
        import tables as pt
        import numarray as nar

        fd = pt.openFile( options.dump_filename, mode = 'w',
                          title = "Dump file" )
        out = state_dp.create_output_dict()
        for key, val in out.iteritems():
            fd.createArray( fd.root, key, nar.asarray( val.data ),
                            '%s data' % val.mode )
        fd.close()

    if options.pert_mesh_filename is not None:
        coors0 = dpb.get_mesh_coors()
        # !!!
        # 'u' is here for displacements of le.py!
        vec_u = state_dp('u').copy()
        vec_u = vec_u.reshape( coors0.shape )
        coors = coors0 + vec_u
        dpb.set_mesh_coors( coors )
        dpb.domain.mesh.write( options.pert_mesh_filename, io = 'auto' )

    return dpb, state_dp, data

def solve_adjoint(conf, options, dpb, state_dp, data):
    """
    Solve the adjoint (linear) problem.
    """
    opts = conf.options

    if dpb:
        apb = dpb.copy('adjoint')

    else:
        apb = ProblemDefinition.from_conf(conf, init_equations=False)

    equations = getattr(conf, '_'.join(('equations_adjoint',
                                        opts.problem,
                                        opts.objective_function)))
    apb.set_equations(equations)
    apb.time_update(None)
    apb.ebcs.zero_dofs()
    apb.update_equations(None, ebcs=apb.ebcs)

    var_data = state_dp.get_parts()
    var_data = remap_dict(var_data, opts.var_map)

    nls_conf = apb.get_solver_conf(opts.nls_adjoint)
    state_ap = apb.solve(nls_conf=nls_conf, var_data=var_data)

    trunk = io.get_trunk(conf.filename_mesh)
    apb.save_state(trunk + '_adjoint.vtk', state_ap)

    shape_opt = so.ShapeOptimFlowCase.from_conf(conf, dpb, apb)

    if options.test is not None:
        ##
        # Test shape sensitivity.
        if shape_opt.test_terms_if_test:
            so.test_terms([options.test], opts.term_delta, shape_opt,
                          var_data, state_ap)

        shape_opt.check_sensitivity([options.test], opts.delta,
                                    var_data, state_ap)
    ##
    # Compute objective function.
    val = shape_opt.obj_fun(state_dp)
    print 'actual obj_fun:', val

    ##
    # Compute shape sensitivity.
    vec_sa = shape_opt.sensitivity(var_data, state_ap)
    print 'actual sensitivity:', vec_sa

##
# c: 22.11.2006, r: 15.04.2008
def solve_optimize( conf, options ):
    opts = conf.options
    trunk = io.get_trunk( conf.filename_mesh )
    data = {}

    dpb = ProblemDefinition.from_conf( conf, init_equations = False )
    equations = getattr( conf, '_'.join( ('equations_direct', opts.problem) ) )

    dpb.set_equations( equations )

    dpb.name = 'direct'
    dpb.time_update(None)

    apb = dpb.copy('adjoint')
    equations = getattr( conf, '_'.join( ('equations_adjoint',
                                          opts.problem,
                                          opts.objective_function) ) )

    apb.set_equations( equations )
    apb.time_update(None)
    apb.ebcs.zero_dofs()
    apb.update_equations(None, ebcs=apb.ebcs)

    ls_conf = dpb.get_solver_conf(opts.ls)
    dnls_conf = dpb.get_solver_conf(opts.nls_direct)
    anls_conf = dpb.get_solver_conf(opts.nls_adjoint)
    opt_conf = dpb.get_solver_conf(opts.optimizer)

    dpb.init_solvers(ls_conf=ls_conf, nls_conf=dnls_conf)

    apb.init_solvers(ls_conf=ls_conf, nls_conf=anls_conf)

    shape_opt = so.ShapeOptimFlowCase.from_conf(conf, dpb, apb)
    design0 = shape_opt.dsg_vars.val
    shape_opt.cache = Struct(design=design0 + 100,
                             state=None,
                             i_mesh=-1)

    opt_status = IndexedStruct()
    optimizer = Solver.any_from_conf(opt_conf,
                                     obj_fun=so.obj_fun,
                                     obj_fun_grad=so.obj_fun_grad,
                                     status=opt_status,
                                     obj_args=(shape_opt, opts))

    ##
    # State problem solution for the initial design.
    vec_dp0 = so.solve_problem_for_design(dpb, design0, shape_opt, opts)

    dpb.save_state( trunk + '_direct_initial.vtk', vec_dp0 )

    ##
    # Optimize.
    des = optimizer( design0 )
    print opt_status

    ##
    # Save final state (for "optimal" design).
    dpb.domain.mesh.write( trunk + '_opt.mesh', io = 'auto' )
    dpb.save_state(trunk + '_direct_current.vtk', shape_opt.cache.state)

    print des

usage = """%prog [options] filename_in"""

help = {
    'server_mode' :
    "run in server mode [default: %default], N/A",
    'adjoint' :
    "solve adjoint problem [default: %default]",
    'direct' :
    "solve direct problem [default: %default]",
    'test' :
    "test sensitivity by finite difference,"
    " using design variable idsg; switches on -a, -d",
    'dump' :
    "dump direct problem state to filename",
    'pert':
    "save displacement-perturbed mesh to filename",
    'optimize' :
    "full shape optimization problem",
}

##
# created:       13.06.2005
# last revision: 15.04.2008
def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "-s", "--server",
                       action = "store_true", dest = "server_mode",
                       default = False, help = help['server_mode'] )
    parser.add_option( "-a", "--adjoint",
                       action = "store_true", dest = "adjoint",
                       default = False, help = help['adjoint'] )
    parser.add_option( "-d", "--direct",
                       action = "store_true", dest = "direct",
                       default = False, help = help['direct'] )
    parser.add_option( "-t", "--test", type = int, metavar = 'idsg',
                       action = "store", dest = "test",
                       default = None, help = help['test'] )
    parser.add_option( "", "--dump", metavar = 'filename',
                       action = "store", dest = "dump_filename",
                       default = None, help = help['dump'] )
    parser.add_option( "", "--pert-mesh", metavar = 'filename',
                       action = "store", dest = "pert_mesh_filename",
                       default = None, help = help['pert'] )
    parser.add_option( "-f", "--full",
                       action = "store_true", dest = "optimize",
                       default = False, help = help['optimize'] )

    options, args = parser.parse_args()

    if options.test is not None:
        options.adjoint = options.direct = True

    if options.optimize:
        options.adjoint = options.direct = False

    if ((len( args ) == 1)
        and (options.direct or options.adjoint or options.optimize)):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    required, other = get_standard_keywords()
    required.remove('equations')
    if options.adjoint:
        required += ['equations_adjoint_.*', 'filename_vp',
                     'equations_direct_.*']
        options.direct = True
    elif options.direct:
        required += ['equations_direct_.*']
    elif options.optimize:
        required += ['equations_direct_.*', 'equations_adjoint_.*',
                     'equations_sensitivity_.*',
                     'filename_vp']

    conf = ProblemConf.from_file( filename_in, required, other )

    if options.direct:
        dpb, state_dp, data = solve_direct( conf, options )
    else:
        dpb, state_dp, data = None, None, None

    if options.adjoint:
        solve_adjoint( conf, options, dpb, state_dp, data )

    if options.optimize:
        solve_optimize( conf, options )

if __name__ == '__main__':
    main()
