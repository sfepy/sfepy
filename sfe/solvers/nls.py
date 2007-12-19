from sfe.base.base import *
from sfe.solvers.solvers import NonlinearSolver
import sfe.base.plotutils as plu

##
# 13.12.2005, c
# 14.12.2005
# 02.10.2007
def checkTangentMatrix( conf, vecX0, mtxA0, evaluator ):
    vecX = vecX0.copy()
    delta = conf.delta

    vecR, status = evaluator.evalResidual( vecX ) # Update state.
    mtxA0, status = evaluator.evalTangentMatrix( vecX, mtxA0 )

    mtxA = mtxA0.tocsc()
    mtxD = mtxA.copy()
    mtxD.data[:] = 0.0

    vecDX = nm.zeros_like( vecR )

    for ic in range( vecDX.shape[0] ):
        vecDX[ic] = delta
        xx = vecX.copy()
        evaluator.updateVec( xx, vecDX )
        vecR1, status = evaluator.evalResidual( xx )

        vecDX[ic] = -delta
        xx = vecX.copy()
        evaluator.updateVec( xx, vecDX )
        vecR2, status = evaluator.evalResidual( xx )

        vecDX[ic] = 0.0;

        vec = 0.5 * (vecR2 - vecR1) / delta

##         ir = mtxA.indices[mtxA.indptr[ic]:mtxA.indptr[ic+1]]
##         for ii in ir:
##             mtxD[ii,ic] = vec[ii]
            
        ir = mtxA.indices[mtxA.indptr[ic]:mtxA.indptr[ic+1]]
        mtxD.data[mtxA.indptr[ic]:mtxA.indptr[ic+1]] = vec[ir]


    vecR, status = evaluator.evalResidual( vecX ) # Restore.

    tt = time.clock()
    print mtxA, mtxD
    plu.plotMatrixDiff( mtxD, mtxA, delta, ['difference', 'analytical'],
                        conf.check )

    return time.clock() - tt

##
# 02.12.2005, c
def convTest( conf, it, err, err0 ):

    status = -1
    if (abs( err0 ) < conf.macheps):
        errR = 0.0
    else:
        errR = err / err0

    print 'nls: iter: %d, out-of-balance: %e (rel: %e)' % (it, err, errR)
    if it > 0:
        if (err < conf.epsA) and (errR < conf.epsR):
            status = 0
    else:
        if err < conf.epsA:
            status = 0

    if (status == -1) and (it >= conf.iMax):
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
    # 02.12.2005, c
    # 13.12.2005
    # 14.12.2005
    # 14.04.2006
    # 18.09.2006
    # 05.10.2006
    # 18.07.2007
    # 02.10.2007
    # 10.10.2007, from newton()
    def __call__( self, vecX0, conf = None, evaluator = None,
                  linSolver = None, status = None ):
        """setting conf.problem == 'linear' means 1 iteration and no rezidual
        check!
        """
        conf = getDefault( conf, self.conf )
        evaluator = getDefault( evaluator, self.evaluator )
        linSolver = getDefault( linSolver, self.linSolver )
        status = getDefault( status, self.status )

        timeStats = {}

        vecX = vecX0.copy()
        vecXLast = vecX0.copy()
        vecDX = None

        err0 = -1.0
        errLast = -1.0
        it = 0
        while 1:

            ls = 1.0
            vecDX0 = vecDX;
            while 1:
                tt = time.clock()
                vecR, ret = evaluator.evalResidual( vecX )
                timeStats['rezidual'] = time.clock() - tt
                if ret == 0: # OK.
                    err = nla.norm( vecR )
                    if it == 0:
                        err0 = err;
                        break
                    if err < (errLast * conf.lsOn): break
                    red = conf.lsRed;
                    print 'linesearch: iter %d, (%.5e < %.5e) (new ls: %e)'\
                          % (it, err, errLast * conf.lsOn, red * ls)
                else: # Failure.
                    red = conf.lsRedWarp;
                    print 'rezidual computation failed for iter %d (new ls: %e)!'\
                          % (it, red * ls) 
                    if (it == 0):
                        raise RuntimeError, 'giving up...'

                if ls < conf.lsMin:
                    if ret != 0:
                        raise RuntimeError, 'giving up...'
                    print 'linesearch failed, continuing anyway'
                    break

                ls *= red;

                vecDX = ls * vecDX0;
                vecX = vecXLast.copy()
                evaluator.updateVec( vecX, vecDX )
            # End residual loop.

            errLast = err;
            vecXLast = vecX.copy()

            condition = convTest( conf, it, err, err0 )
            if condition >= 0:
                break

            tt = time.clock()
            if conf.matrix == 'internal':
                mtxA, ret = evaluator.evalTangentMatrix( vecX )
            else:
                mtxA, ret = evaluator.mtxA, 0
            timeStats['matrix'] = time.clock() - tt
            if ret != 0:
                raise RuntimeError, 'giving up...'

            if conf.check:
                tt = time.clock()
                wt = checkTangentMatrix( conf, vecX, mtxA, evaluator )
                timeStats['check'] = time.clock() - tt - wt
    ##            if conf.check == 2: pause()

    ##         print mtxA.__repr__()
    ##         pause()
    ##        print nla.eig( mtxA.A )
    ##         mtxA.save( 'mtxNSe-5', format = '%d %d %e\n' )
    ##         pause()
    ##         plu.spy( mtxA, eps = 1e-12 )
    ##         plu.pylab.show()

            tt = time.clock() 
            vecDX = linSolver( vecR, mtx = mtxA )
            timeStats['solve'] = time.clock() - tt

            for kv in timeStats.iteritems():
                print '%10s: %7.2f [s]' % kv

            vecE = mtxA * vecDX - vecR
            lerr = nla.norm( vecE )
            if lerr > (conf.epsA * conf.linRed):
                print 'linear system not solved! (err = %e)' % lerr
    #            raise RuntimeError, 'linear system not solved! (err = %e)' % lerr

            evaluator.updateVec( vecX, vecDX )

            if conf.isPlot:
                plu.pylab.ion()
                plu.pylab.gcf().clear()
                plu.pylab.subplot( 2, 2, 1 )
                plu.pylab.plot( vecXLast )
                plu.pylab.ylabel( r'$x_{i-1}$' )
                plu.pylab.subplot( 2, 2, 2 )
                plu.pylab.plot( vecR )
                plu.pylab.ylabel( r'$r$' )
                plu.pylab.subplot( 2, 2, 4 )
                plu.pylab.plot( vecDX )
                plu.pylab.ylabel( r'$\Delta x$' )
                plu.pylab.subplot( 2, 2, 3 )
                plu.pylab.plot( vecX )
                plu.pylab.ylabel( r'$x_i$' )
                plu.pylab.draw()
                plu.pylab.ioff()
                pause()

            if conf.problem == 'linear':
                break

            it += 1

        if status is not None:
            status['timeStats'] = timeStats
            status['err0'] = err0
            status['err'] = err
            status['condition'] = condition

        return vecX
