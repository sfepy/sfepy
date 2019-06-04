from __future__ import print_function
import numpy as nm

try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except (ImportError, RuntimeError):
    plt = mpl = None
    #print 'matplotlib import failed!'

from sfepy.base.base import output, pause

def spy(mtx, eps=None, color='b', **kwargs):
    """
    Show sparsity structure of a `scipy.sparse` matrix.
    """
    aux = mtx.tocoo()
    ij, val = nm.concatenate((aux.row[:,nm.newaxis],
                              aux.col[:,nm.newaxis]), 1), aux.data
    n_item = aux.getnnz()
    n_row, n_col = aux.shape

    if eps is not None:
        output('using eps =', eps)
        ij = nm.compress(nm.absolute(val) > eps, ij, 0)
        n_item = ij.shape[0]
    else:
        output('showing all')

    output('n_item:', n_item)
    if n_item:
        args = {'marker' : '.'}
        args.update(kwargs)
        plt.plot(ij[:,1], ij[:,0], color, linestyle='None', **args)

    plt.axis('image')
    plt.axis([-0.5, n_row+0.5, -0.5, n_col+0.5])
    plt.xlabel(r'%d x %d: %d nnz, %.2f%% fill'
               % (n_row, n_col, n_item, 100. * n_item /
                  (float(n_row) * float(n_col))))
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])

def spy_and_show(mtx, **kwargs):
    spy(mtx, **kwargs)
    plt.show()

##
# 13.12.2005, c
def print_matrix_diff(title, legend, mtx1, mtx2, mtx_da, mtx_dr, iis):
    print('%s: ir, ic, %s, %s, adiff, rdiff' % ((title,) + tuple(legend)))

    aux = mtx_da.copy().tocsc() # mtx_da should be CSC, cast for safety anyway.
    aux.data = nm.ones(mtx_da.data.shape[0])
    ics, irs = aux.nonzero()

    for ii in iis:
        ir, ic = irs[ii], ics[ii]
        print('%5d %5d %11.4e %11.4e %9.2e %9.2e'
              % (ir, ic, mtx1[ir,ic], mtx2[ir,ic],
                 mtx_da[ir,ic], mtx_dr[ir,ic]))

    print('total: %d' % len(iis))

##
# 13.12.2005, c
# 14.12.2005
# 15.12.2005
# 18.07.2007
def plot_matrix_diff(mtx1, mtx2, delta, legend, mode):

    eps = 1e-16

    print("min", legend[0] , legend[1], ":", nm.amin(mtx1.data), nm.amin(mtx2.data))
    print("max", legend[0] , legend[1], ":", nm.amax(mtx1.data), nm.amax(mtx2.data))

    mtx_da = mtx1.copy() # To preserve structure of mtx1.
    mtx_da.data[:] = nm.abs(mtx1.data - mtx2.data)

    mtx_dr = mtx_da.copy()
    mtx_dr.data[:] = -1
    iin = nm.where(nm.abs(mtx1.data) > eps)[0]
    mtx_dr.data[iin] = mtx_da.data[iin] / nm.abs(mtx1.data[iin])

    print("err abs min max:", nm.amin(mtx_da.data), nm.amax(mtx_da.data))
    print("err rel min max:", nm.amin(mtx_dr.data), nm.amax(mtx_dr.data))

    epsilon = max(1e-5, 10 * delta)

    # FIXME - hot fix to see results - findout why does not pause work
    print('epsilon:', epsilon)
    # pause()

    ija = nm.where(mtx_da.data > epsilon)[0]
    print_matrix_diff('--- absolute diff', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ija)
    # pause()

    iin = nm.where(nm.abs(mtx1.data) > epsilon)[0]
    ij = nm.where(nm.abs(mtx_dr.data[iin]) > epsilon)[0]
    ij = iin[ij]
    print_matrix_diff('--- relative diff', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ij)
    # pause()

    ijb = nm.intersect1d(ija, ij)
    print_matrix_diff('--- a-r', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ijb)
    # pause()

    ii = nm.argsort(mtx_dr.data[ijb])
    n_s = min(20, len(ii))
    ijbs = ijb[ii[-1:-n_s-1:-1]]
    print_matrix_diff('--- a-r 20 biggest (by r)', legend,
                     mtx1, mtx2, mtx_da, mtx_dr, ijbs)
    # pause()

    if mode < 2: return

    h = 100
    plt.figure(h); plt.clf()

    plt.axes([0.04, 0.6, 0.3, 0.3], frameon=True)
    spy(mtx_da, epsilon)
    plt.title('absolute diff')

    plt.axes([0.68, 0.6, 0.3, 0.3], frameon=True)
    iia = nm.where(mtx_dr.data)[0]
    mtx_dr.data[nm.setdiff1d(iia, iin)] = 0.0
    spy(mtx_dr, epsilon)
    plt.title('relative diff')

    plt.axes([0.36, 0.6, 0.3, 0.3], frameon=True)
    mtx = mtx_dr.copy()
    mtx.data[:] = 0.0
    ii = nm.intersect1d(nm.where(mtx_dr.data > epsilon)[0],
                           nm.where(mtx_da.data > epsilon)[0])
    mtx.data[ii] = 1.0
    spy(mtx, epsilon)
    plt.title('a-r intersection')

    plt.axes([0.04, 0.08, 0.42, 0.42], frameon=True)
    spy(mtx1, epsilon)
    plt.title(legend[0])

    plt.axes([0.54, 0.08, 0.42, 0.42], frameon=True)
    spy(mtx2, epsilon)
    plt.title(legend[1])

    plt.show()

##
# 02.05.2006, c
def set_axes_font_size(ax, size):
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(size)

##
# 27.09.2006, c
def font_size(size):
    return mpl.font_manager.FontProperties(size=size)

##
# 28.08.2007, c
def iplot(*args, **kwargs):
    plt.ion()
    plt.plot(*args, **kwargs)
    plt.draw()
    plt.ioff()
    pause()
