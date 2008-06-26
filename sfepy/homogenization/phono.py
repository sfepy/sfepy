from sfepy.base.plotutils import pylab

from sfepy.base.base import *
from sfepy.base.la import eig
from sfepy.fem.evaluate import evalTermOP
from sfepy.base.progressbar import MyBar

##
# 26.09.2007, c
# 27.09.2007
# 01.10.2007
# 02.10.2007
def processOptions( options, nEigs ):
    try:
        save = options.saveEigVectors
    except:
        save = (nEigs, nEigs)

    try:
        eigRange = options.eigRange
        if eigRange[-1] < 0:
            eigRange[-1] += nEigs + 1
    except:
        eigRange = (0, nEigs)
    assert eigRange[0] < (eigRange[1] - 1)
    assert eigRange[1] <= nEigs
    
    try:
        freqMargins = 0.01 * nm.array( options.freqMargins, dtype = nm.float64 )
    except:
        freqMargins = nm.array( (0.05, 0.05), dtype = nm.float64 )

    try:
        freqStep = 0.01 * options.freqStep
    except:
        freqStep = 0.05

    try:
        feps = options.feps
    except:
        feps = 1e-8

    try:
        zeps = options.zeps
    except:
        zeps = 1e-8

    try:
        teps = options.teps
    except:
        teps = 1e-4

    try:
        eigVectorTransform = options.eigVectorTransform
    except:
        eigVectorTransform = None

    try:
        plotTranform = options.plotTranform
    except:
        plotTranform = None

    try:
        squared = options.squared
    except:
        squared = True

    return Struct( **locals() )

##
# c: 08.04.2008, r: 08.04.2008
def getMethod( options ):
    if hasattr( options, 'method' ):
        method = options.method
    else:
        method = 'eig.sgscipy'
    return method

##
# created:       27.09.2007
# last revision: 08.04.2008
def computeAverageDensity( pb, matName1, matName2 ):

    mat1 = pb.materials[matName1]
    regionName = mat1.region.name
    vol1 = evalTermOP( None, 'd_volume.i1.%s( u1 )' % regionName, pb )
    mat2 = pb.materials[matName2]
    regionName = mat2.region.name
    vol2 = evalTermOP( None, 'd_volume.i1.%s( u )' % regionName, pb )
    output( 'volumes:', vol1, vol2, vol1 + vol2 )
    output( 'densities:', mat1.density, mat2.density )

    return mat1.density * vol1 + mat2.density * vol2

##
# c: 26.09.2007, r: 08.04.2008
def computeMassComponents( pb, mtxPhi, threshold,
                           transform = None, pbar = None ):
    """Eigenmomenta..."""
    dim = pb.domain.mesh.dim
    nDof, nEigs = mtxPhi.shape
    nNod = nDof / dim
    
    term = pb.equations[0].terms[0]
    matName = term.getMaterialNames()[0]
    mat = pb.materials[matName]
    
    ucName = 'uc'
    dName = 'd'
    dotTerm = 'd_volume_dot.i1.%s( %s, %s )' % (mat.region.name, dName, ucName)

    density = nm.empty( (nDof,), dtype = nm.float64 )
    density.fill( mat.density )
    pb.variables[dName].dataFromData( density, slice( 0, nDof ) )
    
    masses = nm.empty( (nEigs, dim), dtype = nm.float64 )

    if pbar is not None:
        pbar.init( nEigs - 1 )
        
    for ii in xrange( nEigs ):
        if pbar is not None:
            pbar.update( ii )
        else:
            if (ii % 100) == 0:
                output( '%d of %d (%f%%)' % (ii, nEigs,
                                             100. * ii / (nEigs - 1)) )
            
        if transform is None:
            vecPhi, isZero = mtxPhi[:,ii], False
        else:
            vecPhi, isZero = transform( mtxPhi[:,ii], (nNod, dim) )
           
        if isZero:
            masses[ii,:] = 0.0
        else:
            for ir in range( dim ):
                vec = vecPhi[ir::dim].copy()
                pb.variables[ucName].dataFromData( vec, slice( 0, nNod ) )
                val = evalTermOP( None, dotTerm, pb )
                if abs( val ) >= threshold:
                    masses[ii,ir] = val
                else:
                    masses[ii,ir] = 0.0
    #            print ii, ir, val
        
    return masses

##
# 26.09.2007, c
# 27.09.2007
def computeGeneralizedMass( freq, masses, eigs, averageDensity, squared ):
    """Assumes squared freq!"""
    dim = masses.shape[1]
    mtxMass = nm.eye( dim, dim, dtype = nm.float64 ) * averageDensity
    for ir in range( dim ):
        for ic in range( dim ):
            if ir <= ic:
                if squared:
                    val = nm.sum( masses[:,ir] * masses[:,ic]\
                                  / (freq - eigs) )
                    mtxMass[ir,ic] += - freq * val
                else:
                    val = nm.sum( masses[:,ir] * masses[:,ic]\
                                  / ((freq**2) - (eigs**2)) )
                    mtxMass[ir,ic] += - (freq**2) * val
            else:
                mtxMass[ir,ic] = mtxMass[ic,ir]
    return mtxMass


##
# c: 27.09.2007, r: 08.04.2008
def findZero( f0, f1, masses, eigs, averageDensity, opts, mode ):
    feps, zeps = opts.feps, opts.zeps
    method = getMethod( opts )

    fm, fp = f0, f1
    ieig = {0 : 0, 1 : -1}[mode]
    while 1:
        f = 0.5 * (fm + fp)
        mtxMass = computeGeneralizedMass( f, masses, eigs,
                                          averageDensity, opts.squared )
        meigs = eig( mtxMass, eigenvectors = False, method = method )
#        print meigs

        val = meigs[ieig]
#        print f, val, fp - fm

        if (abs( val ) < zeps)\
               or ((fp - fm) < (100.0 * nm.finfo( float ).eps))\
               or ((fp - fm) < feps):
            return 0, f, val

        if mode == 0:
            if (f - f0) < feps:
                return 2, f0, val
            elif (f1 - f) < feps:
                return 1, f, val
        elif mode == 1:
            if (f1 - f) < feps:
                return 1, f, val
            elif (f - f0) < feps:
                return 2, f0, val
            
        if val > 0.0:
            fp = f
        else:
            fm = f

##
# c: 27.09.2007, r: 08.04.2008
def describeGaps( gaps ):
    kinds = []
    for ii, (gmin, gmax) in enumerate( gaps ):

        if (gmin[0] == 2) and (gmax[0] == 2):
            kind = ('p', 'propagation zone')
        elif (gmin[0] == 1) and (gmax[0] == 2):
            kind = ('w', 'full weak band gap')
        elif (gmin[0] == 0) and (gmax[0] == 2):
            kind = ('wp', 'weak band gap + propagation zone')
        elif (gmin[0] == 1) and (gmax[0] == 1):
            kind = ('s', 'full strong band gap (due to end of freq. range or'
                    ' too large thresholds)')
        elif (gmin[0] == 1) and (gmax[0] == 0):
            kind = ('sw', 'strong band gap + weak band gap')
        elif (gmin[0] == 0) and (gmax[0] == 0):
            kind = ('swp', 'strong band gap + weak band gap + propagation zone')
        else:
            output( 'impossible band gap combination:' )
            output( gmin, gmax )
            raise ValueError
        kinds.append( kind )
    return kinds

##
# created:       01.10.2007
# last revision: 13.12.2007
def transformPlotData( datas, plotTranform, funmod ):
    if plotTranform is not None:
        fun = getattr( funmod, plotTranform[0] )

    dmin, dmax = 1e+10, -1e+10
    tdatas = []
    for data in datas:
        tdata = data.copy()
        if plotTranform is not None:
            tdata[:,1:] = fun( tdata[:,1:], *plotTranform[1:] )
        dmin = min( dmin, tdata[:,1:].min() )
        dmax = max( dmax, tdata[:,1:].max() )
        tdatas.append( tdata )
    dmin, dmax = min( dmax - 1e-8, dmin ), max( dmin + 1e-8, dmax )
    return (dmin, dmax), tdatas

##
# c: 27.09.2007, r: 12.06.2008
def plotLogs( figNum, logs, freqRange, plotRange, squared, show = False ):
    if pylab is None: return

    fig = pylab.figure( figNum )
    ax = fig.add_subplot( 111 )

    for f in freqRange:
        l0 = ax.plot( [f, f], plotRange, 'r' )

    for log in logs:
        l1 = ax.plot( log[:,0], log[:,1], 'b--' )
        l2 = ax.plot( log[:,0], log[:,2], 'b-' )

    fmin, fmax = logs[0][0,0], logs[-1][-1,0]
    ax.plot( [fmin, fmax], [0, 0], 'k--' )
    ax.legend( (l0, l1, l2),
               ('eigenfrequencies', 'min eig($A^*$)', 'max eig($A^*$)') )
    if squared:
        ax.set_xlabel( r'$\lambda$, $\omega^2$' )
    else:
        ax.set_xlabel( r'$\sqrt{\lambda}$, $\omega$' )
    ax.set_ylabel( r'eigenvalues of mass matrix $A^*$' )

    if show:
        ax.set_xlim( [fmin, fmax] )
        ax.set_ylim( plotRange )
        pylab.show()
    
##
# c: 27.09.2007, r: 12.06.2008
def plotGaps( figNum, gaps, kinds, freqRange, plotRange, show = False ):
    if pylab is None: return

    def drawRect( ax, x, y, color ):
        ax.fill( nm.asarray( x )[[0,1,1,0]],
                 nm.asarray( y )[[0,0,1,1]],
                 fc = color, linewidth = 0 )

    fig = pylab.figure( figNum )
    ax = fig.add_subplot( 111 )

    # Colors.
    strong = (1, 1, 0.5)
    weak = (1, 1, 1)
    propagation = (0.5, 1, 0.5)

    for ii in xrange( len( freqRange ) - 1 ):
        f0, f1 = freqRange[[ii, ii+1]]
        gmin, gmax = gaps[ii]
        kind, kindDesc = kinds[ii]

        if kind == 'p':
            drawRect( ax, (f0, f1), plotRange, propagation )
            info = [(f0, f1)]
        elif kind == 'w':
            drawRect( ax, (f0, f1), plotRange, weak )
            info = [(f0, f1)]
        elif kind == 'wp':
            drawRect( ax, (f0, gmin[1]), plotRange, weak )
            drawRect( ax, (gmin[1], f1), plotRange, propagation )
            info = [(f0, gmin[1]), (gmin[1], f1)]
        elif kind == 's':
            drawRect( ax, (f0, f1), plotRange, strong )
            info = [(f0, f1)]
        elif kind == 'sw':
            drawRect( ax, (f0, gmax[1]), plotRange, strong )
            drawRect( ax, (gmax[1], f1), plotRange, weak )
            info = [(f0, gmax[1]), (gmax[1], f1)]
        elif kind == 'swp':
            drawRect( ax, (f0, gmax[1]), plotRange, strong )
            drawRect( ax, (gmax[1], gmin[1]), plotRange, weak )
            drawRect( ax, (gmin[1], f1), plotRange, propagation )
            info = [(f0, gmax[1]), (gmax[1], gmin[1]), (gmin[1], f1)]
        else:
            output( 'impossible band gap combination:' )
            output( gmin, gmax )
            raise ValueError

        output( ii, gmin[0], gmax[0], '%.8f' % f0, '%.8f' % f1 )
        output( ' -> %s\n    %s' %(kindDesc, info) )

    if show:
        ax.set_xlim( [freqRange[0], freqRange[-1]] )
        ax.set_ylim( plotRange )
        pylab.show()
    
##
# c: 27.09.2007, r: 08.04.2008
def detectBandGaps( pb, eigs, mtxPhi, conf, options ):
    
    averageDensity = computeAverageDensity( pb, 'matrix', 'inclusion' )
    output( 'average density:', averageDensity )

    nEigs = eigs.shape[0]
    opts = processOptions( conf.options, nEigs )
    method = getMethod( conf.options )
    output( 'method:', method )
    
    if not opts.squared:
        eigs = nm.sqrt( eigs )

    freqRange = eigs[slice( *opts.eigRange )]
    nFreq = freqRange.shape[0]
    minFreq, maxFreq = freqRange[0], freqRange[-1]
    margins = opts.freqMargins * (maxFreq - minFreq)

    prevEig = minFreq - margins[0]
    nextEig = maxFreq + margins[1]
    if opts.eigRange[0] > 0:
        prevEig = max( eigs[opts.eigRange[0]-1] + opts.feps, prevEig )
    if opts.eigRange[1] < nEigs:
        nextEig = min( eigs[opts.eigRange[1]] - opts.feps, nextEig )
    prevEig = max( opts.feps, prevEig )
    nextEig = max( opts.feps, nextEig, prevEig + opts.feps )
    freqRangeMargins = nm.r_[prevEig, freqRange, nextEig]

    output( 'freq. range             : [%8.3f, %8.3f]' % (minFreq, maxFreq) )
    output( 'freq. range with margins: [%8.3f, %8.3f]'\
          % tuple( freqRangeMargins[[0,-1]] ) )

##     print freqRange
##     print freqRangeMargins
##     pause()

    if opts.eigVectorTransform is not None:
        fun = getattr( conf.funmod, opts.eigVectorTransform[0] )
        def _wrapTransform( vec, shape ):
            return fun( vec, shape, *opts.eigVectorTransform[1:] )
    else:
        _wrapTransform = None
    output( 'mass matrix components...')
    pbar = MyBar( 'computing:' )
    tt = time.clock()
    masses = computeMassComponents( pb, mtxPhi, opts.teps, _wrapTransform, pbar )
    output( '...done in %.2f s' % (time.clock() - tt) )

    logs = []
    gaps = []
    df = opts.freqStep * (maxFreq - minFreq)
    for ii in xrange( nFreq + 1 ):
        f0, f1 = freqRangeMargins[[ii, ii+1]]
        output( 'interval: ]%.8f, %.8f[...' % (f0, f1) )

        log = []
        num = max( 5, (f1 - f0) / df )
        logFreqs = nm.linspace( f0 + opts.feps, f1 - opts.feps, num )
        for f in logFreqs:
            mtxMass = computeGeneralizedMass( f, masses, eigs,
                                              averageDensity, opts.squared )
            meigs = eig( mtxMass, eigenvectors = False, method = method )
            log.append( [f, meigs[0], meigs[-1]] )
        log0, log1 = log[0], log[-1]
        if log0[1] > 0.0: # No gaps.
            gap = ([2, f0, log0[1]], [2, f0, log0[2]])
        elif log1[2] < 0.0: # Full interval strog gap.
            gap = ([1, f1, log1[1]], [1, f1, log1[2]])
        else:
            output( 'finding zero of the largest eig...' )
            smax, fmax, vmax = findZero( f0, f1, masses, eigs,
                                         averageDensity, opts, 1 )
            output( '...done' )
            if smax in [0, 2]:
                output( 'finding zero of the smallest eig...' )
                smin, fmin, vmin = findZero( fmax, f1, masses, eigs,
                                             averageDensity, opts, 0 )
                output( '...done' )
            elif smax == 1:
                smin = 1 # both are negative everywhere.
                fmin, vmin = fmax, vmax

            gap = ([smin, fmin, vmin], [smax, fmax, vmax])

        output( gap[0] )
        output( gap[1] )
#        pause()
        gaps.append( gap )
        logs.append( nm.array( log, dtype = nm.float64 ) )
        output( '...done' )

    kinds = describeGaps( gaps )

    return Struct( logs = logs, gaps = gaps, kinds = kinds,
                   freqRange = freqRange, freqRangeMargins = freqRangeMargins,
                   opts = opts )
