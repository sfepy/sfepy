from sfepy.base.base import *
from sfepy.solvers.solvers import OptimizationSolver
from sfepy.base.log import Log

import scipy.optimize as sopt
import scipy.optimize.linesearch as linesearch

##
# 19.04.2006, c
# 26.04.2006
# 28.04.2006
def conv_test( conf, it, of, of0, ofg_norm = None ):
    """
    -1 ... continue
     0 ... small OF -> stop
     1 ... i_max reached -> stop
     2 ... small OFG -> stop
     3 ... small relative decrase of OF
     """

    status = -1
    output( 'opt: iter: %d, of: %e (||ofg||: %e)' % (it, of, ofg_norm) )
#    print (of0 - of), (conf.eps_rd * of0)

    if (abs( of ) < conf.eps_of):
        status = 0
    elif ofg_norm and (ofg_norm < conf.eps_ofg):
        status = 2
    elif (it > 0) and (abs(of0 - of) < (conf.eps_rd * abs( of0 ))):
        status = 3
        
    if (status == -1) and (it >= conf.i_max):
        status = 1

    return status

##
# 19.04.2006, from scipy.optimize
# 21.04.2006
# 27.03.2007
def wrap_function( function, args ):
    ncalls = [0]
    times = []
    def function_wrapper( x ):
        ncalls[0] += 1
        tt = time.time()
        out = function( x, *args )
        tt2 = time.time()
        if tt2 < tt:
            raise RuntimeError, '%f >= %f' % (tt, tt2)
        times.append( tt2 - tt )
        return out
    return ncalls, times, function_wrapper

##
# 20.04.2006, c
def check_gradient( xit, aofg, fn_of, delta, check ):

    dofg = nm.zeros_like( aofg )
    xd = xit.copy()
    for ii in xrange( xit.shape[0] ):
        xd[ii] = xit[ii] + delta
        ofp = fn_of( xd )

        xd[ii] = xit[ii] - delta
        ofm = fn_of( xd )

        xd[ii] = xit[ii]

        dofg[ii] = 0.5 * (ofp - ofm) / delta

        output( '**********', ii, aofg[ii], dofg[ii] )

    diff = abs( aofg - dofg )
    aux = nm.concatenate( (aofg[:,nm.newaxis], dofg[:,nm.newaxis],
                           diff[:,nm.newaxis]), 1 )
    output( aux )
    output( nla.norm( diff, nm.Inf ) )
    aofg.tofile( 'aofg.txt', ' ' )
    dofg.tofile( 'dofg.txt', ' ' )
    diff.tofile( 'diff.txt', ' ' )
    if check == 2:
        import pylab
        pylab.plot( aofg )
        pylab.plot( dofg )
        pylab.legend( ('analytical', 'finite difference') )
        pylab.show()
    pause( 'gradient checking done' )

##
# 17.10.2007, c
class FMinSteepestDescent( OptimizationSolver ):
    name = 'opt.fmin_sd'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_0 = {
            'name'      : 'fmin_sd',
            'kind'      : 'opt.fmin_sd',

            'i_max'      : 10,
            'eps_rd'     : 1e-5, # Relative delta of objective function
            'eps_of'     : 1e-4,
            'eps_ofg'    : 1e-8,
            'norm'      : nm.Inf,
            'ls'        : True, # Linesearch.
            'ls_method'  : 'backtracking', # 'backtracking' or 'full'
            'ls0'       : 0.25,
            'ls_red'     : 0.5,
            'ls_red_warp' : 0.1,
            'ls_on'      : 0.99999,
            'ls_min'     : 1e-5,
            'check'     : 0,
            'delta'     : 1e-6,
            'output'    : None, # 'itc'
            'log'       : True,
            'yscales'   : ['linear', 'log', 'log'],
        }
        """
        get = conf.get_default_attr

        i_max = get( 'i_max', 10 )
        eps_rd = get( 'eps_rd', 1e-5 )
        eps_of = get( 'eps_of', 1e-4 )
        eps_ofg = get( 'eps_ofg', 1e-8 )
        norm = get( 'norm', nm.Inf )
        ls = get( 'ls', True )
        ls_method = get( 'ls_method', 'backtracking' )
        ls0 = get( 'ls0', 0.25 )
        ls_red = get( 'ls_red', 0.5 )
        ls_red_warp = get( 'ls_red_warp', 0.1 )
        ls_on = get( 'ls_on', 0.99999 )
        ls_min = get( 'ls_min', 1e-5 )
        check = get( 'check', 0 )
        delta = get( 'delta', 1e-6)
        output = get( 'output', None )
        log = get( 'log', True )
        yscales = get( 'yscales', ['linear', 'log', 'log'] )

        common = OptimizationSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

    ##
    # 17.10.2007, c
    def __init__( self, conf, **kwargs ):
        OptimizationSolver.__init__( self, conf, **kwargs )

    ##
    # 19.04.2006, c
    # 20.04.2006
    # 21.04.2006
    # 26.04.2006
    # 06.06.2006
    # 07.06.2006
    # 04.09.2006
    # 21.03.2007
    # 17.10.2007, from fmin_sd()
    def __call__( self, x0, conf = None, obj_fun = None, obj_fun_grad = None,
                  status = None, obj_args = None ):
#    def fmin_sd( conf, x0, fn_of, fn_ofg, args = () ):

        conf = get_default( conf, self.conf )
        obj_fun = get_default( obj_fun, self.obj_fun )
        obj_fun_grad = get_default( obj_fun_grad, self.obj_fun_grad )
        status = get_default( status, self.status )
        obj_args = get_default( obj_args, self.obj_args )

        if conf.output:
            globals()['output'] = conf.output

        output( 'entering optimization loop...' )

        nc_of, tt_of, fn_of = wrap_function( obj_fun, obj_args )
        nc_ofg, tt_ofg, fn_ofg = wrap_function( obj_fun_grad, obj_args )

        time_stats = {'of' : tt_of, 'ofg': tt_ofg, 'check' : []}

        if conf.log:
            log = Log.from_conf( conf, ([r'of'], [r'$||$ofg$||$'], [r'alpha']) )
        else:
            log = None

        ofg = None

        it = 0
        xit = x0.copy()
        while 1:

            of = fn_of( xit )

            if it == 0:
                of0 = ofit0 = of_prev = of
                of_prev_prev = of + 5000.0

            if ofg is None:
    #            ofg = 1
                ofg = fn_ofg( xit )

            if conf.check:
                tt = time.clock()
                check_gradient( xit, ofg, fn_of, conf.delta, conf.check )
                time_stats['check'].append( time.clock() - tt )

            ofg_norm = nla.norm( ofg, conf.norm )

            ret = conv_test( conf, it, of, ofit0, ofg_norm )
            if ret >= 0:
                break
            ofit0 = of

            ##
            # Backtrack (on errors).
            alpha = conf.ls0
            can_ls = True
            while 1:
                xit2 = xit - alpha * ofg
                aux = fn_of( xit2 )
                if aux is None:
                    alpha *= conf.ls_red_warp
                    can_ls = False
                    output( 'warp: reducing step (%f)' % alpha )
                elif conf.ls and conf.ls_method == 'backtracking':
                    if aux < of * conf.ls_on: break
                    alpha *= conf.ls_red
                    output( 'backtracking: reducing step (%f)' % alpha )
                else:
                    of_prev_prev = of_prev
                    of_prev = aux
                    break

                if alpha < conf.ls_min:
                    if aux is None:
                        raise RuntimeError, 'giving up...'
                    output( 'linesearch failed, continuing anyway' )
                    break

            # These values are modified by the line search, even if it fails
            of_prev_bak = of_prev
            of_prev_prev_bak = of_prev_prev

            if conf.ls and can_ls and conf.ls_method == 'full':
                output( 'full linesearch...' )
                alpha, fc, gc, of_prev, of_prev_prev, ofg1 = \
                       linesearch.line_search(fn_of,fn_ofg,xit,
                                              -ofg,ofg,of_prev,of_prev_prev,
                                              c2=0.4)
                if alpha is None:  # line search failed -- use different one.
                    alpha, fc, gc, of_prev, of_prev_prev, ofg1 = \
                           sopt.line_search(fn_of,fn_ofg,xit,
                                            -ofg,ofg,of_prev_bak,
                                            of_prev_prev_bak)
                    if alpha is None or alpha == 0:
                        # This line search also failed to find a better solution.
                        ret = 3
                        break
                output( ' -> alpha: %.8e' % alpha )
            else:
                if conf.ls_method == 'full':
                    output( 'full linesearch off (%s and %s)' % (conf.ls,
                                                                 can_ls) )
                ofg1 = None

            if conf.log:
                log( of, ofg_norm, alpha )

            xit = xit - alpha * ofg
            if ofg1 is None:
                ofg = None
            else:
                ofg = ofg1.copy()

            for key, val in time_stats.iteritems():
                if len( val ):
                    output( '%10s: %7.2f [s]' % (key, val[-1]) )

            it = it + 1

        output( 'status:               %d' % ret )
        output( 'initial value:        %.8e' % of0 )
        output( 'current value:        %.8e' % of )
        output( 'iterations:           %d' % it )
        output( 'function evaluations: %d in %.2f [s]' \
              % (nc_of[0], nm.sum( time_stats['of'] ) ) )
        output( 'gradient evaluations: %d in %.2f [s]' \
              % (nc_ofg[0], nm.sum( time_stats['ofg'] ) ) )

        if conf.log:
            log( of, ofg_norm, alpha, finished = True )

        if status is not None:
            status['log'] = log
            status['status'] = status
            status['of0'] = of0
            status['of'] = of
            status['it'] = it
            status['nc_of'] = nc_of[0]
            status['nc_ofg'] = nc_ofg[0]
            status['time_stats'] = time_stats

        return xit
