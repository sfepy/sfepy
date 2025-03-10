from itertools import product

import numpy as nm
from numpy.lib.stride_tricks import as_strided
import numpy.linalg as nla
import scipy as sc

from sfepy.base.base import assert_, output

def norm_l2_along_axis(ar, axis=1, n_item=None, squared=False):
    """Compute l2 norm of rows (axis=1) or columns (axis=0) of a 2D array.

    n_item ... use only the first n_item columns/rows
    squared ... if True, return the norm squared
    """
    assert_(axis in [0, 1])
    assert_(ar.ndim == 2)

    other = 1 - axis
    vec = nm.zeros((ar.shape[other],), dtype=nm.float64)

    if n_item is None:
        n_item = ar.shape[axis]
    else:
        n_item = min( n_item, ar.shape[axis] )

    if axis == 1:
        for ii in range( n_item ):
            vec += ar[:,ii]**2
    else:
        for ii in range( n_item ):
            vec += ar[ii,:]**2

    if not squared:
        vec = nm.sqrt( vec )

    return vec

def normalize_vectors(vecs, eps=1e-8):
    """
    Normalize an array of vectors in place.

    Parameters
    ----------
    vecs : array
        The 2D array of vectors in rows.
    eps : float
        The tolerance for considering a vector to have zero norm. Such
        vectors are left unchanged.
    """
    norms = norm_l2_along_axis(vecs, axis=1)
    ii = norms > eps
    vecs[ii] = vecs[ii] / norms[ii][:, None]

def dets_fast(a):
    """
    Fast determinant calculation of 3-dimensional array.

    Parameters
    ----------
    a : array
        The input array with shape (m, n, n).

    Returns
    -------
    out : array
        The output array with shape (m,): out[i] = det(a[i, :, :]).
    """
    from packaging import version

    if version.parse(nm.__version__) >= version.parse('1.8'):
        return nm.linalg.det(a)

    else:
        from numpy.linalg import lapack_lite
        from numpy.core import intc

        a = a.copy()
        m = a.shape[0]
        n = a.shape[1]
        lapack_routine = lapack_lite.dgetrf
        pivots = nm.zeros((m, n), intc)
        flags = nm.arange(1, n + 1).reshape(1, -1)
        for i in range(m):
            tmp = a[i]
            lapack_routine(n, n, tmp, n, pivots[i], 0)
        sign = 1. - 2. * (nm.add.reduce(pivots != flags, axis=1) % 2)
        idx = nm.arange(n)
        d = a[:, idx, idx]
        absd = nm.absolute(d)
        sign *= nm.multiply.reduce(d / absd, axis=1)
        nm.log(absd, absd)
        logdet = nm.add.reduce(absd, axis=-1)

        return sign * nm.exp(logdet)


def invs_fast(a, det=None):
    """
    Fast inversion calculation of 4-dimensional array.

    Parameters
    ----------
    a : array
        The input array with shape (c, q, n, n).
    det: array
        To speed up the calculation, enter the already calculated determinant.

    Returns
    -------
    out : array
        The output array with shape (c, q, n, n):
        out[c, q] = inv(a[c, q, :, :]).
    """
    ax = nm.einsum("ij...->...ij", a)
    dim = a.shape[-1]
    inv_ax = nm.empty_like(ax)
    if det is None:
        det_a = dets_fast(a)[..., None, None]
    else:
        det_a = det.reshape(a.shape[:2] + (1, 1))

    if dim == 3:
        inv_ax[0, 0] = -ax[1, 2] * ax[2, 1] + ax[1, 1] * ax[2, 2]
        inv_ax[1, 0] = ax[1, 2] * ax[2, 0] - ax[1, 0] * ax[2, 2]
        inv_ax[2, 0] = -ax[1, 1] * ax[2, 0] + ax[1, 0] * ax[2, 1]

        inv_ax[0, 1] = ax[0, 2] * ax[2, 1] - ax[0, 1] * ax[2, 2]
        inv_ax[1, 1] = -ax[0, 2] * ax[2, 0] + ax[0, 0] * ax[2, 2]
        inv_ax[2, 1] = ax[0, 1] * ax[2, 0] - ax[0, 0] * ax[2, 1]

        inv_ax[0, 2] = -ax[0, 2] * ax[1, 1] + ax[0, 1] * ax[1, 2]
        inv_ax[1, 2] = ax[0, 2] * ax[1, 0] - ax[0, 0] * ax[1, 2]
        inv_ax[2, 2] = -ax[0, 1] * ax[1, 0] + ax[0, 0] * ax[1, 1]
    elif dim == 2:
        inv_ax[0, 0] = ax[1, 1]
        inv_ax[1, 0] = -ax[1, 0]
        inv_ax[0, 1] = -ax[0, 1]
        inv_ax[1, 1] = ax[0, 0]
    elif dim == 1:
        inv_ax[0, 0] = 1.
    else:
        raise NotImplementedError(f'matrix dimension {dim}x{dim}')

    return nm.einsum("...ij->ij...", inv_ax) / det_a

def get_blocks_stats(blocks, *args):
    """
    Return statistics of array/matrix `blocks` defined by indices in `args`.

    Returns
    -------
    stats: structured array
        The array with 'shape', 'min', 'mean' and 'max' fields at positions of
        each matrix block.

    Examples
    --------
    >>> import numpy as nm
    >>> from sfepy.linalg.utils import get_blocks_stats
    >>>
    >>> A = nm.eye(3)
    >>> B = nm.full((3,2), 2)
    >>> C = nm.full((1,3), 3)
    >>> D = nm.full((1,2), 4)
    >>> M = nm.block([[A, B], [C, D]])
    >>>
    >>> sr = [slice(0, 3), slice(3, 5)]
    >>> sc = [slice(0, 3), slice(3, 4)]
    >>> stats = get_blocks_stats(M, sr, sc)
    >>>
    >>> print(stats['shape'])
    [[(3, 3) (3, 1)]
     [(1, 3) (1, 1)]]
    >>>
    >>> print(stats['min'])
    [[0. 2.]
     [3. 4.]]
    """
    bindices = args
    idim = len(args)
    bdim = blocks.ndim
    bshape = [len(indices) for indices in bindices]
    if idim == 1:
        bshape = bshape * blocks.ndim
        bindices *= blocks.ndim

    elif idim != bdim:
        raise ValueError('wrong number of dimensions of block indices!'
                         f' (can be 1 or {bdim}, is {idim})')

    dt = blocks.dtype
    sizes = nm.empty(bshape,
                     dtype=[('shape', tuple), ('min', dt), ('mean', dt),
                            ('max', dt), ('maxabs', dt)])
    for iflat, ii in enumerate(product(*bindices)):
        key = nm.unravel_index(iflat, bshape)
        block = blocks[ii]
        sizes[key] = ((block.shape, block.min(), block.mean(), block.max(),
                       nm.abs(block).max()))
    return sizes

def print_array_info(ar):
    """
    Print array shape and other basic information.
    """
    ar = nm.asanyarray(ar)

    print(ar.shape, 'c_contiguous:', ar.flags.c_contiguous, \
          'f_contiguous:', ar.flags.f_contiguous)
    print('min:', ar.min(), 'mean:', ar.mean(), 'max:', ar.max())

def output_array_stats(ar, name, verbose=True):
    ar = nm.asarray(ar)
    if not len(ar):
        output('%s: empty' % name)

    elif nm.isrealobj(ar):
        output('%s\nmin: % .6e mean: % .6e median: % .6e max: % .6e'
               % (name, ar.min(), ar.mean(), nm.median(ar), ar.max()),
               verbose=verbose)

    else:
        output_array_stats(ar.real, 'Re(%s)' % name , verbose=verbose)
        output_array_stats(ar.imag, 'Im(%s)' % name, verbose=verbose)

def max_diff_csr(mtx1, mtx2):
    aux = nm.abs((mtx1 - mtx2).data)
    return aux.max() if len(aux) else 0.0

##
# 14.12.2005, c
def cycle(bounds):
    """
    Cycles through all combinations of bounds, returns a generator.

    More specifically, let bounds=[a, b, c, ...], so cycle returns all
    combinations of lists [0<=i<a, 0<=j<b, 0<=k<c, ...] for all i,j,k,...

    Examples:
    In [9]: list(cycle([3, 2]))
    Out[9]: [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]

    In [14]: list(cycle([3, 4]))
    [[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [2, 0],
    [2, 1], [2, 2], [2, 3]]

    """
    yield from product(*map(range, bounds))

def assemble1d(ar_out, indx, ar_in):
    """
    Perform `ar_out[indx] += ar_in`, where items of `ar_in`
    corresponding to duplicate indices in `indx` are summed together.
    """
    if len(indx) > 0:
        zz = nm.zeros_like(indx)
        aux = sc.sparse.coo_matrix((ar_in, (indx, zz)), dtype=ar_in.dtype)
        aux = aux.tocsr().tocoo() # This sums the duplicates.

        ar_out[aux.row] += aux.data

def unique_rows(ar, return_index=False, return_inverse=False):
    """
    Return unique rows of a two-dimensional array `ar`. The arguments follow
    `numpy.unique()`.
    """
    ar = nm.ascontiguousarray(ar)

    # View the rows as a 1D structured array.
    arv = ar.view(ar.shape[1] * [('', ar.dtype)])
    out = nm.unique(arv, return_index=return_index,
                    return_inverse=return_inverse)
    if isinstance(out, tuple):
        uarv = out[0]

    else:
        uarv = out

    # Restore the original dimensions.
    uar = uarv.view(ar.dtype).reshape((-1, ar.shape[1]))

    if isinstance(out, tuple):
        out = (uar,) + out[1:]

    else:
        out = uar

    return out

def argsort_rows(seq):
    """
    Returns an index array that sorts the sequence `seq`. Works along
    rows if `seq` is two-dimensional.
    """
    seq = nm.asanyarray(seq)
    if seq.ndim == 1:
        ii = nm.argsort(seq)

    else:
        ii = nm.lexsort(seq.T[::-1])

    return ii

def map_permutations(seq1, seq2, check_same_items=False):
    """
    Returns an index array `imap` such that `seq1[imap] == seq2`, if
    both sequences have the same items - this is not checked by default!

    In other words, finds the indices of items of `seq2` in `seq1`.
    """
    assert_(len(seq1) == len(seq2))

    seq1 = nm.asanyarray(seq1)
    seq2 = nm.asanyarray(seq2)

    i1 = argsort_rows(seq1)
    i2 = argsort_rows(seq2)

    if check_same_items:
        assert_(seq1.shape == seq2.shape)
        assert_((seq1[i1] == seq2[i2]).all())

    ii = nm.argsort(i2)

    imap = i1[ii]

    return imap

def mini_newton( fun, x0, dfun, i_max = 100, eps = 1e-8 ):
    x = x0
    ii = 0
    while ii < i_max:
        r = fun( x )
        err = nla.norm( r )
##         print ii, x, r, err
        if err < eps: break

        mtx = dfun( x )
        try:
            dx = nm.dot( nla.inv( mtx.T ), r )
        except:
            break
        x = x - dx
        ii += 1
    return x

def chunk_arrays(arrs, chunk_size):
    """
    Yield consecutive views into arrays `arrs` of the same lengths and most
    `chunk_size` long.
    """
    ii = 0
    while ii < len(arrs[0]):
        yield [arr[ii:ii+chunk_size] for arr in arrs]

        ii += chunk_size

def insert_strided_axis(ar, axis, length):
    """
    Insert a new axis of given length into an array using numpy stride
    tricks, i.e. no copy is made.

    Parameters
    ----------
    ar : array
        The input array.
    axis : int
        The axis before which the new axis will be inserted.
    length : int
        The length of the inserted axis.

    Returns
    -------
    out : array
        The output array sharing data with `ar`.

    Examples
    --------
    >>> import numpy as nm
    >>> from sfepy.linalg import insert_strided_axis
    >>> ar = nm.random.rand(2, 1, 2)
    >>> ar
    array([[[ 0.18905119,  0.44552425]],

           [[ 0.78593989,  0.71852473]]])
    >>> ar.shape
    (2, 1, 2)
    >>> ar2 = insert_strided_axis(ar, 1, 3)
    >>> ar2
    array([[[[ 0.18905119,  0.44552425]],

            [[ 0.18905119,  0.44552425]],

            [[ 0.18905119,  0.44552425]]],


           [[[ 0.78593989,  0.71852473]],

            [[ 0.78593989,  0.71852473]],

            [[ 0.78593989,  0.71852473]]]])
    >>> ar2.shape
    (2, 3, 1, 2)
    """
    shape = list(ar.shape)
    shape.insert(axis, length)

    strides = list(ar.strides)
    strides.insert(axis, 0)

    out = as_strided(ar, shape=shape, strides=strides)
    return out

def dot_sequences(mtx, vec, mode='AB'):
    """
    Computes dot product for each pair of items in the two sequences.

    Equivalent to

    >>> out = nm.empty((vec.shape[0], mtx.shape[1], vec.shape[2]),
    >>>                dtype=vec.dtype)
    >>> for ir in range(mtx.shape[0]):
    >>>     out[ir] = nm.dot(mtx[ir], vec[ir])

    Parameters
    ----------
    mtx : array
        The array of matrices with shape `(n_item, m, n)`.
    vec : array
        The array of vectors with shape `(n_item, a)` or matrices with shape
        `(n_item, a, b)`.
    mode : one of 'AB', 'ATB', 'ABT', 'ATBT'

        The mode of the dot product - the corresponding axes are dotted
        together:

        'AB'   : `a = n`
        'ATB'  : `a = m`
        'ABT'  : `b = n` (*)
        'ATBT' : `b = m` (*)

        (*) The 'BT' part is ignored for the vector second argument.

    Returns
    -------
    out : array
       The resulting array.

    Notes
    -----
    Uses `numpy.matmul()` via the `@` operator.
    """
    if vec.ndim == mtx.ndim:
        squeeze = False

    else:
        squeeze = True
        vec = vec[..., None]

    if 'BT' in mode:
        ax = list(range(vec.ndim))
        vec = vec.transpose((ax[:-2]) + [ax[-1], ax[-2]])

    if 'AT' in mode:
        ax = list(range(mtx.ndim))
        mtx = mtx.transpose((ax[:-2]) + [ax[-1], ax[-2]])

    out = mtx @ vec

    if squeeze:
        out = out[..., 0]

    return out

def apply_to_sequence(seq, fun, ndim, out_item_shape):
    """
    Applies function `fun()` to each item of the sequence `seq`. An item
    corresponds to the last `ndim` dimensions of `seq`.

    Parameters
    ----------
    seq : array
        The sequence array with shape `(n_1, ..., n_r, m_1, ..., m_{ndim})`.
    fun : function
        The function taking an array argument of shape of length `ndim`.
    ndim : int
        The number of dimensions of an item in `seq`.
    out_item_shape : tuple
        The shape an output item.

    Returns
    -------
    out : array
       The resulting array of shape `(n_1, ..., n_r) + out_item_shape`. The
       `out_item_shape` must be compatible with the `fun`.
    """
    n_seq = nm.prod(seq.shape[0:-ndim], dtype=int)
    aux = nm.reshape(seq, (n_seq,) + seq.shape[-ndim:])

    out = nm.empty((n_seq,) + out_item_shape, dtype=seq.dtype)
    for ii, item in enumerate(aux):
        out[ii,:] = fun(item)

    out = nm.reshape(out, seq.shape[0:-ndim] + out_item_shape)

    return out
