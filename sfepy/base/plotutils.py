from base import *

try:
    import matplotlib as mpl
    mpl.use('GTKAgg')
    import pylab
except (ImportError, RuntimeError):
    pylab = None
    #print 'matplotlib import failed!'

##
# 13.12.2005, c
# 16.12.2005
# 04.08.2006
# 22.08.2006
def spy( mtx, eps = None, color = 'b', **kwargs ):
    aux = mtx.tocoo()
    ij, val = nm.concatenate( (aux.row[:,nm.newaxis],
                               aux.col[:,nm.newaxis]), 1 ), aux.data
    n_item = aux.getnnz()
    n_row, n_col = aux.shape

    if eps is not None:
        print 'using', eps
        ij = nm.compress( nm.absolute( val ) > eps, ij, 0 )
        n_item = ij.shape[0]
    else:
        print 'showing all'

    print n_item
    if n_item:
        args = {'marker' : '.', 'markersize' : 0.5, 'markeredgewidth' : 0.1}
        args.update( kwargs )
        pylab.plot( ij[:,1] + 0.5, ij[:,0] + 0.5, color, linestyle = 'None',
                    **args )
    pylab.axis( [-0.5, n_row+0.5, -0.5, n_col+0.5] )
    pylab.axis( 'image' )
    pylab.xlabel( '%d x %d: %d nnz, %.2f\%% fill'
                  % (n_row, n_col, n_item, 100. * n_item /
                     (float( n_row ) * float( n_col )) ) )
    ax = pylab.gca()
    ax.set_ylim( ax.get_ylim()[::-1] )

##
# 13.12.2005, c
def print_matrix_diff( title, legend, mtx1, mtx2, mtx_da, mtx_dr, iis ):
    print '%s: ir, ic, %s, %s, adiff, rdiff' % ((title,) + tuple( legend ))
    for ii in iis:
        [ir, ic] = mtx_da.rowcol( ii )
        print '%5d %5d %11.4e %11.4e %9.2e %9.2e'\
              % (ir, ic, mtx1[ir,ic], mtx2[ir,ic], mtx_da[ir,ic], mtx_dr[ir,ic] )

    print 'total: %d' % len( iis )

##
# 13.12.2005, c
# 14.12.2005
# 15.12.2005
# 18.07.2007
def plot_matrix_diff( mtx1, mtx2, delta, legend, mode ):

    eps = 1e-16

    print nm.amin( mtx1.data ), nm.amin( mtx2.data )
    print nm.amax( mtx1.data ), nm.amax( mtx2.data )

    mtx_da = mtx1.copy() # To preserve structure of mtx1.
    mtx_da.data[:] = nm.abs( mtx1.data - mtx2.data )

    mtx_dr = mtx_da.copy()
    mtx_dr.data[:] = -1
    iin = nm.where( nm.abs( mtx1.data ) > eps )[0]
    mtx_dr.data[iin] = mtx_da.data[iin] / nm.abs( mtx1.data[iin] )

    print nm.amin( mtx_da.data ), nm.amax( mtx_da.data )
    print nm.amin( mtx_dr.data ), nm.amax( mtx_dr.data )

    epsilon = max( 1e-5, 10 * delta )

    print 'epsilon:', epsilon
    pause()

    ija = nm.where( mtx_da.data > epsilon )[0]
    print_matrix_diff( '--- absolute diff', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ija )
    pause()

    iin = nm.where( nm.abs( mtx1.data ) > epsilon )[0]
    ij = nm.where( nm.abs( mtx_dr.data[iin] ) > epsilon )[0]
    ij = iin[ij]
    print_matrix_diff( '--- relative diff', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ij )
    pause()

    ijb = nm.intersect1d( ija, ij )
    print_matrix_diff( '--- a-r', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ijb )
    pause()

    ii = nm.argsort( mtx_dr.data[ijb] )
    n_s = min( 20, len( ii ) )
    ijbs = ijb[ii[-1:-n_s-1:-1]]
    print_matrix_diff( '--- a-r 20 biggest (by r)', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ijbs )
    pause()

    if mode < 2: return
    
    h = 100
    pylab.figure( h ); pylab.clf

    pylab.axes( [0.04, 0.6, 0.3, 0.3], frameon = True )
    spy( mtx_da, epsilon )
    pylab.title( 'absolute diff' )

    pylab.axes( [0.68, 0.6, 0.3, 0.3], frameon = True )
    iia = nm.where( mtx_dr.data )[0]
    mtx_dr.data[nm.setdiff1d( iia, iin )] = 0.0
    spy( mtx_dr, epsilon )
    pylab.title( 'relative diff' )

    pylab.axes( [0.36, 0.6, 0.3, 0.3], frameon = True )
    mtx = mtx_dr.copy()
    mtx.data[:] = 0.0
    ii = nm.intersect1d( nm.where( mtx_dr.data > epsilon )[0],
                           nm.where( mtx_da.data > epsilon )[0] )
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
def set_axes_font_size( ax, size ):
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size( size )

##
# 27.09.2006, c
def font_size( size ):
    return mpl.font_manager.FontProperties( size = size )

##
# 28.08.2007, c
def iplot( *args, **kwargs ):
    pylab.ion()
    pylab.plot( *args, **kwargs )
    pylab.draw()
    pylab.ioff()
    pause()
