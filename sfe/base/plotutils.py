from base import *

try:
    import pylab
    import matplotlib as mpl
except ImportError:
    print 'matplotlib import failed!'

##
# 13.12.2005, c
# 16.12.2005
# 04.08.2006
# 22.08.2006
def spy( mtx, eps = None, color = 'b', **kwargs ):
    aux = mtx.tocoo()
    ij, val = nm.concatenate( (aux.row[:,nm.newaxis],
                               aux.col[:,nm.newaxis]), 1 ), aux.data
    nItem = aux.getnnz()
    nRow, nCol = aux.shape

    if eps is not None:
        print 'using', eps
        ij = nm.compress( nm.absolute( val ) > eps, ij, 0 )
        nItem = ij.shape[0]
    else:
        print 'showing all'

    print nItem
    if nItem:
        args = {'marker' : '.', 'markersize' : 0.5, 'markeredgewidth' : 0.1}
        args.update( kwargs )
        pylab.plot( ij[:,1] + 0.5, ij[:,0] + 0.5, color, linestyle = 'None',
                    **args )
    pylab.axis( [-0.5, nRow+0.5, -0.5, nCol+0.5] )
    pylab.axis( 'image' )
    pylab.xlabel( '%d x %d: %d nnz, %.2f\%% fill'
                  % (nRow, nCol, nItem, 100. * nItem /
                     (float( nRow ) * float( nCol )) ) )
    ax = pylab.gca()
    ax.set_ylim( ax.get_ylim()[::-1] )

##
# 13.12.2005, c
def printMatrixDiff( title, legend, mtx1, mtx2, mtxDA, mtxDR, iis ):
    print '%s: ir, ic, %s, %s, adiff, rdiff' % ((title,) + tuple( legend ))
    for ii in iis:
        [ir, ic] = mtxDA.rowcol( ii )
        print '%5d %5d %11.4e %11.4e %9.2e %9.2e'\
              % (ir, ic, mtx1[ir,ic], mtx2[ir,ic], mtxDA[ir,ic], mtxDR[ir,ic] )

    print 'total: %d' % len( iis )

##
# 13.12.2005, c
# 14.12.2005
# 15.12.2005
# 18.07.2007
def plotMatrixDiff( mtx1, mtx2, delta, legend, mode ):

    eps = 1e-16

    print nm.amin( mtx1.data ), nm.amin( mtx2.data )
    print nm.amax( mtx1.data ), nm.amax( mtx2.data )

    mtxDA = mtx1.copy() # To preserve structure of mtx1.
    mtxDA.data[:] = nm.abs( mtx1.data - mtx2.data )

    mtxDR = mtxDA.copy()
    mtxDR.data[:] = -1
    iin = nm.where( nm.abs( mtx1.data ) > eps )[0]
    mtxDR.data[iin] = mtxDA.data[iin] / nm.abs( mtx1.data[iin] )

    print nm.amin( mtxDA.data ), nm.amax( mtxDA.data )
    print nm.amin( mtxDR.data ), nm.amax( mtxDR.data )

    epsilon = max( 1e-5, 10 * delta )

    print 'epsilon:', epsilon
    pause()

    ija = nm.where( mtxDA.data > epsilon )[0]
    printMatrixDiff( '--- absolute diff', legend,
                     mtx1, mtx2, mtxDA, mtxDR, ija )
    pause()

    iin = nm.where( nm.abs( mtx1.data ) > epsilon )[0]
    ij = nm.where( nm.abs( mtxDR.data[iin] ) > epsilon )[0]
    ij = iin[ij]
    printMatrixDiff( '--- relative diff', legend,
                     mtx1, mtx2, mtxDA, mtxDR, ij )
    pause()

    ijb = nm.intersect1d( ija, ij )
    printMatrixDiff( '--- a-r', legend,
                     mtx1, mtx2, mtxDA, mtxDR, ijb )
    pause()

    ii = nm.argsort( mtxDR.data[ijb] )
    nS = min( 20, len( ii ) )
    ijbs = ijb[ii[-1:-nS-1:-1]]
    printMatrixDiff( '--- a-r 20 biggest (by r)', legend,
                     mtx1, mtx2, mtxDA, mtxDR, ijbs )
    pause()

    if mode < 2: return
    
    h = 100
    pylab.figure( h ); pylab.clf

    pylab.axes( [0.04, 0.6, 0.3, 0.3], frameon = True )
    spy( mtxDA, epsilon )
    pylab.title( 'absolute diff' )

    pylab.axes( [0.68, 0.6, 0.3, 0.3], frameon = True )
    iia = nm.where( mtxDR.data )[0]
    mtxDR.data[nm.setdiff1d( iia, iin )] = 0.0
    spy( mtxDR, epsilon )
    pylab.title( 'relative diff' )

    pylab.axes( [0.36, 0.6, 0.3, 0.3], frameon = True )
    mtx = mtxDR.copy()
    mtx.data[:] = 0.0
    ii = nm.intersect1d( nm.where( mtxDR.data > epsilon )[0],
                           nm.where( mtxDA.data > epsilon )[0] )
    mtx.data[ii] = 1.0
    spy( mtx, epsilon )
    pylab.title( 'a-r intersection' )

    pylab.axes( [0.04, 0.08, 0.42, 0.42], frameon = True )
    spy( mtx1, epsilon )
    pylab.title( legend[0] )

    pylab.axes( [0.54, 0.08, 0.42, 0.42], frameon = True )
    spy( mtx2, epsilon )
    pylab.title( legend[1] )

    pylab.show()

##
# 02.05.2006, c
def setAxesFontSize( ax, size ):
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size( size )

##
# 27.09.2006, c
def fontSize( size ):
    return mpl.font_manager.FontProperties( size = size )

##
# 28.08.2007, c
def iplot( *args, **kwargs ):
    pylab.ion()
    pylab.plot( *args, **kwargs )
    pylab.draw()
    pylab.ioff()
    pause()
