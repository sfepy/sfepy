from sfe.base.base import *
from sfe.solvers.solvers import NonlinearSolver
from nls import convTest
import sfe.base.plotutils as plu

##
# 26.07.2007, c
def areClose( a, b, rtol = 0.2, atol = 1e-8 ):
    return False
#    return abs( a - b ) <= max( atol, rtol * abs( b ) )

##
# 26.07.2007, c
def scaleMatrix( mtx, indx, factor ):
    ptr0 = mtx.indptr[indx.start]
    ptr1 = mtx.indptr[indx.stop]
    mtx.data[ptr0:ptr1] *= factor

_dimaterModes = {'edge' : 0, 'volume' : 1, 'max' : 2}

##
# c: 01.08.2007, r: 15.01.2008
def createStabilData( problem, fluidName, stabilName, eqName1, eqName2 ):

    ns = {}
    term = problem.equations[eqName1].terms['dw_lin_convect']

    ns['fluid'] = 'fluid'
    ns['v'] = term.getVirtualName()
    ns['b'] = term.getParameterNames()[0]
    ns['u'] = term.getStateNames()[0]
    ns['omega'] = term.region.name
    
    term = problem.equations[eqName1].terms['dw_grad']
    ns['p'] = term.getStateNames()[0]

    term = problem.equations[eqName2].terms['dw_div']
    ns['q'] = term.getVirtualName()

    ii = {}
    ii['u'] = problem.variables.getIndx( ns['u'] )
    ii['us'] = problem.variables.getIndx( ns['u'], stripped = True )
    ii['ps'] = problem.variables.getIndx( ns['p'], stripped = True )

    stabil = problem.materials[stabilName]
    mat = problem.materials[ns['fluid']]

    viscosity = mat.viscosity

    cFriedrichs = problem.domain.getDiameter()
    sigma = 1e-12 # 1 / dt.

#    print cFriedrichs

    def matFun( ts, coor, region, ig, bNorm = 1.0, fixedData = None ):
        if fixedData is not None:
            return fixedData[ig]

        print '|b|_max (matFun):', bNorm
        gamma = viscosity + bNorm * cFriedrichs

        data = {}
        if stabil.gamma is None:
            data['gamma'] = stabil.gammaMul * gamma
        else:
            data['gamma'] = nm.asarray( stabil.gammaMul * stabil.gamma,
                                        dtype = nm.float64 )
        

        if stabil.delta is None:
            term = problem.equations[eqName1].terms['dw_lin_convect']
            for ig in term.iterGroups():
                # This sets term.ig - for 1 group only!!!
                break
            var = problem.variables[ns['u']]
            ap, vg = var.getApproximation( term.getCurrentGroup(), 'Volume' )
            delta = 1.0
            mode = _dimaterModes[stabil.diameterMode]
            cells = stabil.region.getCells( ig )
            diameters2 = problem.domain.getElementDiameters( ig, cells, vg,
                                                             mode )
            val1 = min( 1.0, 1.0 / sigma )
            val2 = sigma * cFriedrichs**2
            val3 = (bNorm**2) * min( (cFriedrichs**2) / viscosity, 1.0 / sigma )
#            print val1, gamma, val2, val3
            delta = stabil.deltaMul * val1 * diameters2 / (gamma + val2 + val3)
            data['diameters2'] = diameters2
            data['delta'] = delta
        else:
            val = stabil.deltaMul * stabil.delta
            for ii in range( len( stabil.igs ) ):
                data['delta'] = nm.asarray( val, dtype = nm.float64 )
        
        if stabil.tau is None:
            data['tau'] = stabil.tauRed * data['delta']
        else:
            data['tau'] = nm.asarray( stabil.tauMul * stabil.tau,
                                      dtype = nm.float64 )

        return data

    stabil.setFunction( matFun )

    return stabil, ns, ii

##
# 11.10.2007, c
class Oseen( NonlinearSolver ):
    name = 'nls.oseen'

    ##
    # 10.10.2007, c
    def __init__( self, conf, **kwargs ):
        NonlinearSolver.__init__( self, conf, **kwargs )

    ##
    # 26.07.2007, c
    # 31.07.2007
    # 01.08.2007
    # 11.10.2007, from oseen()
    # 29.10.2007
    # 30.10.2007
    # 31.10.2007
    def __call__( self, vecX0, conf = None, evaluator = None,
                  linSolver = None, status = None ):
        """"""
        conf = getDefault( conf, self.conf )
        evaluator = getDefault( evaluator, self.evaluator )
        linSolver = getDefault( linSolver, self.linSolver )
        status = getDefault( status, self.status )

        if hasattr( conf, 'fixedData' ):
            fixedData = conf.fixedData
        else:
            fixedData = None

        timeStats = {}

        problem = evaluator.problem

        stabil, ns, ii = createStabilData( problem, conf.fluidMatName,
                                           conf.stabilMatName,
                                           conf.linConvectEqName,
                                           conf.divEqName )
        updateVar = problem.variables.nonStateDataFromState

        print 'problem size:'
        print '    velocity: %s' % ii['us']
        print '    pressure: %s' % ii['ps']

        vecX = vecX0.copy()
        vecXPrev = vecX0.copy()
        vecDX = None

        err0 = -1.0
        it = 0
        while 1:
            updateVar( ns['b'], vecXPrev, ns['u'] )
            vecB = vecXPrev[ii['u']]
            bNorm = nla.norm( vecB, nm.inf )
            print '|b|_max: %.12e' % bNorm

            vecU = vecX[ii['u']]
            uNorm = nla.norm( vecU, nm.inf )
            print '|u|_max: %.2e' % uNorm

            stabil.timeUpdate( None, None, problem.domain,
                               bNorm = bNorm, fixedData = fixedData )
            maxPars = stabil.reduceOnDatas( lambda a, b: max( a, b.max() ) )
            print 'stabilization parameters:'
            print '                   gamma: %.12e' % maxPars['gamma']
            print '            max( delta ): %.12e' % maxPars['delta']
            print '              max( tau ): %.12e' % maxPars['tau']
            try:
                print '              max( h^2 ): %.12e' % maxPars['diameters2']
            except:
                pass

            if (not areClose( bNorm, 1.0 )) and conf.adimensionalize:
                adimensionalize = True
            else:
                adimensionalize = False

            tt = time.clock()
            vecR, ret = evaluator.evalResidual( vecX )
            timeStats['rezidual'] = time.clock() - tt
            if ret == 0: # OK.
                err = nla.norm( vecR )
                if it == 0:
                    err0 = err;
                else:
                    err += nla.norm( vecDX )
            else: # Failure.
                print 'rezidual computation failed for iter %d!' % it
                raise RuntimeError, 'giving up...'

            condition = convTest( conf, it, err, err0 )
            if condition >= 0:
                break

            if adimensionalize:
                print 'adimensionalizing'
                mat.viscosity = viscosity / bNorm
                vecR[indxUS] /= bNorm

            tt = time.clock()
            mtxA, ret = evaluator.evalTangentMatrix( vecX )
            timeStats['matrix'] = time.clock() - tt
            if ret != 0:
                raise RuntimeError, 'giving up...'

            tt = time.clock() 
            vecDX = linSolver( vecR, mtx = mtxA )
            timeStats['solve'] = time.clock() - tt

            vecE = mtxA * vecDX - vecR
            lerr = nla.norm( vecE )
            if lerr > (conf.epsA * conf.linRed):
                print 'linear system not solved! (err = %e)' % lerr
    #            raise RuntimeError, 'linear system not solved! (err = %e)' % lerr

            if adimensionalize:
                print 'restoring pressure...'
                vecDX[indxPS] *= bNorm

            dxNorm = nla.norm( vecDX )
            print '||dx||: %.2e' % dxNorm

            for kv in timeStats.iteritems():
                print '%10s: %7.2f [s]' % kv


            vecXPrev = vecX.copy()
            evaluator.updateVec( vecX, vecDX )

            if conf.isPlot:
                plu.pylab.ion()
                plu.pylab.gcf().clear()
                plu.pylab.subplot( 2, 2, 1 )
                plu.pylab.plot( vecXPrev )
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

            it += 1

        if conf.checkNavierStokesRezidual:
    ##         updateVar( bName, vecXPrev, uName )
    ## #        updateVar( bName, vecX, uName )
    ##         vecRNS1, ret = residual( vecX, context )
    ##         errNS = nla.norm( vecRNS1 )
    ##         print '"Oseen" rezidual: %.8e' % errNS


            t1 = '+ dw_div_grad.%s( %s, %s, %s )' % (ns['omega'],
                                                     ns['fluid'],
                                                     ns['v'], ns['u'])
##             t2 = '+ dw_lin_convect.%s( %s, %s, %s )' % (ns['omega'],
##                                                         ns['v'], bName, ns['u'])
            t2 = '+ dw_convect.%s( %s, %s )' % (ns['omega'], ns['v'], ns['u'])
            t3 = '- dw_grad.%s( %s, %s )' % (ns['omega'], ns['v'], ns['p'])
            t4 = 'dw_div.%s( %s, %s )' % (ns['omega'], ns['q'], ns['u'])
            equations = {
                'balance' : ' '.join( (t1, t2, t3) ),
                'incompressibility' : t4,
            }
            problem.setEquations( equations )
            vecRNS0, ret = evaluator.evalResidual( vecX0 )
            vecRNS, ret = evaluator.evalResidual( vecX )
            if ret:
                print 'Navier-Stokes rezidual computation failed!'
                errNS = errNS0 = None
            else:
                errNS0 = nla.norm( vecRNS0 )
                errNS = nla.norm( vecRNS )
            print 'Navier-Stokes rezidual0: %.8e' % errNS0
            print 'Navier-Stokes rezidual : %.8e' % errNS
            print 'b - u: %.8e' % nla.norm( vecB - vecU )
            print condition
    ##         print vecRNS - vecRNS1
            plu.pylab.ion()
            plu.pylab.gcf().clear()
            plu.pylab.plot( vecRNS )
    ##         plu.pylab.gcf().clear()
    ##         plu.pylab.plot( vecRNS1 )
            plu.pylab.draw()
            plu.pylab.ioff()
            pause()
        else:
            errNS = None

        if status is not None:
            status['timeStats'] = timeStats
            status['err0'] = err0
            status['err'] = err
            status['errNS'] = errNS
            status['condition'] = condition

        return vecX
