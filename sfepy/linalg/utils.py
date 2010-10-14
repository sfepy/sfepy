import numpy as nm
import numpy.linalg as nla
import scipy as sc

from sfepy.base.base import assert_, insert_method, Struct

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
        for ii in xrange( n_item ):
            vec += ar[:,ii]**2
    else:
        for ii in xrange( n_item ):
            vec += ar[ii,:]**2

    if not squared:
        vec = nm.sqrt( vec )

    return vec

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
        for ii in xrange( ls ):
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
        for ii in xrange( bounds[0] ):
            yield [ii]
    else:
        for ii in xrange( bounds[0] ):
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

    i1 = argsort_rows(seq1)
    i2 = argsort_rows(seq2)

    if check_same_items:
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

def dot_sequences(mtx, vec, use_rows=False):
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
        The array of vectors with shape `(n_item, n)` or matrices with shape
        `(n_item, n, k)`.
    use_rows : bool
        If `vec` is an array of matrices with `(n_item, k, n)` shape, 

    Returns
    -------
    out : array
       The resulting array of shape `(n_item, m)` or `(n_item, m, k)`.

    Notes
    -----
    For r-D arrays `(n_1, ..., n_r, ?, ?)` the arrays are first reshaped to
    `(n_1 * ... * n_r, ?, ?)`, then the dot is performed, and finally the shape
    is restored to `(n_1, ..., n_r, ?, ?)`.
    """
    if (vec.ndim == 2) and (mtx.ndim == 3):
        out = nm.sum(mtx * vec[:,None,:], axis=2)

    elif (vec.ndim == 3) and (mtx.ndim == 3):
        out = nm.empty((vec.shape[0], mtx.shape[1], vec.shape[2]),
                       dtype=vec.dtype)

        if use_rows:
            for ic in range(vec.shape[2]):
                out[:,:,ic] = dot_sequences(mtx, vec[:,:,ic])

        else:
            for ic in range(vec.shape[2]):
                out[:,:,ic] = dot_sequences(mtx, vec[:,ic,:])

    elif (vec.ndim >= 4) and (mtx.ndim >= 4) and (vec.ndim == mtx.ndim):
        mtx_seq = nm.reshape(mtx,
                             (nm.prod(mtx.shape[0:-2], dtype=int),)
                             + mtx.shape[-2:])

        vec_seq = nm.reshape(vec,
                             (nm.prod(vec.shape[0:-2], dtype=int),)
                             + vec.shape[-2:])

        out_seq = dot_sequences(mtx_seq, vec_seq, use_rows=use_rows)
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
            print 'cannot make array from MatrixAction of kind %s!' % self.kind
            raise ValueError
            
