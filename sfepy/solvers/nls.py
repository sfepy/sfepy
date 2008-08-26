from sfepy.base.base import *
from sfepy.solvers.solvers import NonlinearSolver
import sfepy.base.plotutils as plu

##
# 13.12.2005, c
# 14.12.2005
# 02.10.2007
def check_tangent_matrix( conf, vec_x0, mtx_a0, evaluator ):
    vec_x = vec_x0.copy()
    delta = conf.delta

    vec_r, status = evaluator.eval_residual( vec_x ) # Update state.
    mtx_a0, status = evaluator.eval_tangent_matrix( vec_x, mtx_a0 )

    mtx_a = mtx_a0.tocsc()
    mtx_d = mtx_a.copy()
    mtx_d.data[:] = 0.0

    vec_dx = nm.zeros_like( vec_r )

    for ic in range( vec_dx.shape[0] ):
        vec_dx[ic] = delta
        xx = vec_x.copy()
        evaluator.update_vec( xx, vec_dx )
        vec_r1, status = evaluator.eval_residual( xx )

        vec_dx[ic] = -delta
        xx = vec_x.copy()
        evaluator.update_vec( xx, vec_dx )
        vec_r2, status = evaluator.eval_residual( xx )

        vec_dx[ic] = 0.0;

        vec = 0.5 * (vec_r2 - vec_r1) / delta

##         ir = mtx_a.indices[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]]
##         for ii in ir:
##             mtx_d[ii,ic] = vec[ii]

        ir = mtx_a.indices[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]]
        mtx_d.data[mtx_a.indptr[ic]:mtx_a.indptr[ic+1]] = vec[ir]


    vec_r, status = evaluator.eval_residual( vec_x ) # Restore.

    tt = time.clock()
    print mtx_a, '.. analytical'
    print mtx_d, '.. difference'
    plu.plot_matrix_diff( mtx_d, mtx_a, delta, ['difference', 'analytical'],
                        conf.check )

    return time.clock() - tt

##
# c: 02.12.2005, r: 02.04.2008
def conv_test( conf, it, err, err0 ):

    status = -1
    if (abs( err0 ) < conf.macheps):
        err_r = 0.0
    else:
        err_r = err / err0

    output( 'nls: iter: %d, residual: %e (rel: %e)' % (it, err, err_r) )
    if it > 0:
        if (err < conf.eps_a) and (err_r < conf.eps_r):
            status = 0
    else:
        if err < conf.eps_a:
            status = 0

    if (status == -1) and (it >= conf.i_max):
        status = 1

    return status

##
# 10.10.2007, c
class Newton( NonlinearSolver ):
    name = 'nls.newton'

    ##
    # 10.10.2007, c
    def __init__( self, conf, **kwargs ):
        NonlinearSolver.__init__( self, conf, **kwargs )

    ##
    # c: 02.12.2005, r: 04.04.2008
    # 10.10.2007, from newton()
    def __call__( self, vec_x0, conf = None, evaluator = None,
                  lin_solver = None, status = None ):
        """setting conf.problem == 'linear' means 1 iteration and no rezidual
        check!
        """
        conf = get_default( conf, self.conf )
        evaluator = get_default( evaluator, self.evaluator )
        lin_solver = get_default( lin_solver, self.lin_solver )
        status = get_default( status, self.status )

        time_stats = {}

        vec_x = vec_x0.copy()
        vec_x_last = vec_x0.copy()
        vec_dx = None

        err0 = -1.0
        err_last = -1.0
        it = 0
        while 1:

            ls = 1.0
            vec_dx0 = vec_dx;
            while 1:
                tt = time.clock()
                vec_r, ret = evaluator.eval_residual( vec_x )
                time_stats['rezidual'] = time.clock() - tt
                if ret == 0: # OK.
                    try:
                        err = nla.norm( vec_r )
                    except:
                        output( 'infs or nans in the residual:', vec_r )
                        output( nm.isfinite( vec_r ).all() )
                        debug()
                    if it == 0:
                        err0 = err;
                        break
                    if err < (err_last * conf.ls_on): break
                    red = conf.ls_red;
                    output( 'linesearch: iter %d, (%.5e < %.5e) (new ls: %e)'\
                            % (it, err, err_last * conf.ls_on, red * ls) )
                else: # Failure.
                    red = conf.ls_red_warp;
                    output(  'rezidual computation failed for iter %d'
                             ' (new ls: %e)!' % (it, red * ls) )
                    if (it == 0):
                        raise RuntimeError, 'giving up...'

                if ls < conf.ls_min:
                    if ret != 0:
                        raise RuntimeError, 'giving up...'
                    output( 'linesearch failed, continuing anyway' )
                    break

                ls *= red;

                vec_dx = ls * vec_dx0;
                vec_x = vec_x_last.copy()
                evaluator.update_vec( vec_x, vec_dx )
            # End residual loop.

            err_last = err;
            vec_x_last = vec_x.copy()

            condition = conv_test( conf, it, err, err0 )
            if condition >= 0:
                break

            tt = time.clock()
            if conf.problem == 'nonlinear':
                mtx_a, ret = evaluator.eval_tangent_matrix( vec_x )
            else:
                mtx_a, ret = evaluator.mtx, 0
            time_stats['matrix'] = time.clock() - tt
            if ret != 0:
                raise RuntimeError, 'giving up...'

            if conf.check:
                tt = time.clock()
                wt = check_tangent_matrix( conf, vec_x, mtx_a, evaluator )
                time_stats['check'] = time.clock() - tt - wt
    ##            if conf.check == 2: pause()

            tt = time.clock() 
            vec_dx = lin_solver( vec_r, mtx = mtx_a )
            time_stats['solve'] = time.clock() - tt

            for kv in time_stats.iteritems():
                output( '%10s: %7.2f [s]' % kv )

            vec_e = mtx_a * vec_dx - vec_r
            lerr = nla.norm( vec_e )
            if lerr > (conf.eps_a * conf.lin_red):
                output( 'linear system not solved! (err = %e)' % lerr )
    #            raise RuntimeError, 'linear system not solved! (err = %e)' % lerr

            evaluator.update_vec( vec_x, vec_dx )

            if conf.is_plot:
                plu.pylab.ion()
                plu.pylab.gcf().clear()
                plu.pylab.subplot( 2, 2, 1 )
                plu.pylab.plot( vec_x_last )
                plu.pylab.ylabel( r'$x_{i-1}$' )
                plu.pylab.subplot( 2, 2, 2 )
                plu.pylab.plot( vec_r )
                plu.pylab.ylabel( r'$r$' )
                plu.pylab.subplot( 2, 2, 4 )
                plu.pylab.plot( vec_dx )
                plu.pylab.ylabel( r'$\_delta x$' )
                plu.pylab.subplot( 2, 2, 3 )
                plu.pylab.plot( vec_x )
                plu.pylab.ylabel( r'$x_i$' )
                plu.pylab.draw()
                plu.pylab.ioff()
                pause()

            it += 1

##         import pylab as p
##         problem = evaluator.problem
##         r0 = problem.variables.make_full_vec( vec_r, force_value = 0.0 )
##         dx = nm.zeros_like( vec_dx )
##         ii = problem.variables.get_indx( 'r', stripped = True )
##         dx[ii] = 1.0
##         r1 = problem.variables.make_full_vec( mtx_a * dx, force_value = 0.0 )
##         p.plot( r0 )
##         p.plot( r1 )

##         vv = nm.where( nm.abs( r1 ) > 1e-12, 1.0, 0.0 )
##         problem.save_state_to_vtk( 'sd.vtk', vv )
##         nodes = problem.variables.get_nodes_of_global_dofs( nm.where( vv > 0.5 )[0] )
##         print nodes
## #        problem.save_regions( 'asdsd' )
##         p.show()

        if status is not None:
            status['time_stats'] = time_stats
            status['err0'] = err0
            status['err'] = err
            status['condition'] = condition

        return vec_x
