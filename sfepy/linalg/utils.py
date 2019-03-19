from __future__ import print_function
from __future__ import absolute_import
import numpy as nm
from numpy.lib.stride_tricks import as_strided
import numpy.linalg as nla
import scipy as sc
from six.moves import range

try:
    from numpy.core.umath_tests import matrix_multiply

except:
    matrix_multiply = None

from sfepy.base.base import assert_, insert_method, output, Struct

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
    from distutils.version import LooseVersion

    if LooseVersion(nm.__version__) >= '1.8':
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
    if nm.isrealobj(ar):
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
# 21.11.2005, c
def split_range( n_item, step ):
    num = n_item / step
    out = [step] * num
    aux = sum( out )
    if aux < n_item:
        out.append( n_item - aux )

    return out

##
# Inspired on net (ASPN Recipec).
# 14.12.2005, c
def permutations( seq ):

    ls = len( seq )

    if ls <= 1:
        yield seq
    else:
        for ii in range( ls ):
            for perm in permutations( seq[:ii] + seq[ii+1:] ):
                yield [seq[ii]] + perm

##
# 14.12.2005, c
def cycle( bounds ):
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

    nb  = len( bounds )
    if nb == 1:
        for ii in range( bounds[0] ):
            yield [ii]
    else:
        for ii in range( bounds[0] ):
            for perm in cycle( bounds[1:] ):
                yield [ii] + perm

def combine( seqs ):
    """Same as cycle, but with general sequences.

    Example:

    In [19]: c = combine( [['a', 'x'], ['b', 'c'], ['dd']] )

    In [20]: list(c)
    Out[20]: [['a', 'b', 'dd'], ['a', 'c', 'dd'], ['x', 'b', 'dd'],
    ['x', 'c', 'dd']]
    """
    nb  = len( seqs )
    if nb == 1:
        for ii in seqs[0]:
            yield [ii]
    else:
        for ii in seqs[0]:
            for perm in combine( seqs[1:] ):
                yield [ii] + perm

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
    Uses `numpy.core.umath_tests.matrix_multiply()` if available, which is much
    faster than the default implementation.

    The default implementation uses `numpy.sum()` and element-wise
    multiplication. For r-D arrays `(n_1, ..., n_r, ?, ?)` the arrays
    are first reshaped to `(n_1 * ... * n_r, ?, ?)`, then the dot is
    performed, and finally the shape is restored to `(n_1, ..., n_r, ?, ?)`.
    """
    if matrix_multiply is not None:
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

        out = matrix_multiply(mtx, vec)
        if squeeze:
            out = out[..., 0]

    else:
        if (vec.ndim == 2) and (mtx.ndim == 3):
            if mode in ('AB', 'ABT'):
                out = nm.sum(mtx * vec[:, None, :], axis=2)

            else:
                out = nm.sum(mtx * vec[:, :, None], axis=1)

        elif (vec.ndim == 3) and (mtx.ndim == 3):

            if mode == 'AB':
                out = nm.empty((vec.shape[0], mtx.shape[1], vec.shape[2]),
                               dtype=vec.dtype)

                for ic in range(vec.shape[2]):
                    out[:, :, ic] = dot_sequences(mtx, vec[:, :, ic], mode=mode)

            elif mode == 'ABT':
                out = nm.empty((vec.shape[0], mtx.shape[1], vec.shape[1]),
                               dtype=vec.dtype)

                for ic in range(vec.shape[1]):
                    out[:, :, ic] = dot_sequences(mtx, vec[:, ic, :], mode=mode)


            elif mode == 'ATB':
                out = nm.empty((vec.shape[0], mtx.shape[2], vec.shape[2]),
                               dtype=vec.dtype)

                for ic in range(vec.shape[2]):
                    out[:, :, ic] = dot_sequences(mtx, vec[:, :, ic], mode=mode)

            elif mode == 'ATBT':
                out = nm.empty((vec.shape[0], mtx.shape[2], vec.shape[1]),
                               dtype=vec.dtype)

                for ic in range(vec.shape[1]):
                    out[:, :, ic] = dot_sequences(mtx, vec[:, ic, :], mode=mode)

            else:
                raise ValueError('unknown dot mode! (%s)' % mode)

        elif (vec.ndim >= 4) and (mtx.ndim >= 4) and (vec.ndim == mtx.ndim):
            mtx_seq = nm.reshape(mtx,
                                 (nm.prod(mtx.shape[0:-2], dtype=int),)
                                 + mtx.shape[-2:])

            vec_seq = nm.reshape(vec,
                                 (nm.prod(vec.shape[0:-2], dtype=int),)
                                 + vec.shape[-2:])

            out_seq = dot_sequences(mtx_seq, vec_seq, mode=mode)
            out = nm.reshape(out_seq, mtx.shape[0:-2] + out_seq.shape[-2:])

        else:
            raise ValueError('unsupported operand shape')

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

##
# 30.08.2007, c
class MatrixAction( Struct ):

    ##
    # 30.08.2007, c
    def from_function( fun, expected_shape, dtype ):
        def call( self, vec ):
            aux = fun( vec )
            assert_( aux.shape[0] == self.shape[0] )
            return nm.asanyarray( aux, dtype = self.dtype )
        obj = MatrixAction( shape = expected_shape,
                            dtype = dtype,
                            kind = 'function' )
        insert_method( obj, call )
        return obj
    from_function = staticmethod( from_function )

    ##
    # 30.08.2007, c
    def from_array( arr ):
        def call( self, vec ):
            return nm.asarray( sc.dot( self.arr, vec ) )
        obj = MatrixAction( shape = arr.shape,
                            dtype = arr.dtype,
                            arr = arr,
                            kind = 'array' )
        insert_method( obj, call )
        return obj
    from_array = staticmethod( from_array )

    ##
    # 30.08.2007, c
    def __call__( self, vec ):
        return self.call( vec )

    ##
    # 30.08.2007, c
    def to_array( self ):
        if self.kind == 'array':
            return self.arr
        else:
            print('cannot make array from MatrixAction of kind %s!' % self.kind)
            raise ValueError
