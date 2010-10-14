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

def inverse_element_mapping( coors, e_coors, eval_base, ref_coors,
                             suppress_errors = False ):
    """
    Given spatial element coordinates, find the inverse mapping for
    points with coordinats X = X(xi), i.e. xi = xi(X).
        Return:
            xi : reference element coordinates
    """
    n_v, dim = e_coors.shape
    if coors.ndim == 2:
        n_c, dim2 = coors.shape
    else:
        n_c, dim2 = 1, coors.shape[0]

    assert_( dim == dim2 )

    if n_v == (dim + 1): # Simplex.
        bc = barycentric_coors( coors, e_coors )
        xi = nm.dot( bc.T, ref_coors )
        
    else: # Tensor-product and other.
        def residual( xi ):
            bf = eval_base(xi[nm.newaxis,:],
                           suppress_errors=suppress_errors).squeeze()
            res = coors - nm.dot( bf, e_coors )
            return res.squeeze()
        
        def matrix( xi ):
            bfg = eval_base(xi[nm.newaxis,:], diff=True,
                            suppress_errors=suppress_errors).squeeze()
            mtx = - nm.dot( bfg, e_coors )
            return mtx

        xi0 = nm.array([0.0, 0.0, 0.0])
        xi = mini_newton( residual, xi0, matrix )

    return xi

##
# 01.09.2007, c
def barycentric_coors( coors, s_coors ):
    """
    Simplex elements:
        Return:
            bc : barycentric (area in 2D, volume in 3D) coordinates
        Then reference element coordinates xi = dot(bc.T, ref_coors).
    """
    n_v, dim = s_coors.shape
    n_c, dim2 = coors.shape
    assert_( dim == dim2 )
    assert_( n_v == (dim + 1) )
    mtx = nm.ones( (n_v, n_v), nm.float64 )
    mtx[0:n_v-1,:] = s_coors.T
    rhs = nm.empty( (n_v,n_c), nm.float64 )
    rhs[0:n_v-1,:] = coors.T
    rhs[n_v-1,:] = 1.0
    bc = nla.solve( mtx, rhs )
##     print bc.T
    return bc

##
# 30.08.2007, c
# 01.09.2007
def points_in_simplex( coors, s_coors, eps = 1e-8 ):
    n_c, dim = coors.shape
    bc = barycentric_coors( coors, s_coors )
    flag = nm.ones( (n_c,), dtype = nm.bool )
    for idim in xrange( dim + 1 ):
        flag &= nm.where( (bc[idim,:] > -eps)
                          & (bc[idim,:] < (1.0 + eps)), True, False )
    return flag

##
# c: 18.01.2008, r: 18.01.2008
def rotation_matrix2d( angle ):
    angle *= nm.pi / 180.0
    mtx = nm.array( [[nm.cos( angle ), -nm.sin( angle )],
                     [nm.sin( angle ), nm.cos( angle )]], dtype = nm.float64 )
    return mtx

def make_axis_rotation_matrix(direction, angle):
    """
    Create a rotation matrix corresponding to the rotation around a general
    axis by a specified angle.
    
    R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)
    
    Parameters:

    direction : array
        The rotation axis direction vector "d".
    angle : float
        The rotation angle "a".
    """
    d = nm.array(direction, dtype=nm.float64)
    d /= nm.linalg.norm(d)
    
    eye = nm.eye(3, dtype=nm.float64)
    ddt = nm.outer(d, d)
    skew = nm.array([[    0,  d[2],  -d[1]],
                     [-d[2],     0,  d[0]],
                     [d[1], -d[0],    0]], dtype=nm.float64)

    mtx = ddt + nm.cos(angle) * (eye - ddt) + nm.sin(angle) * skew
    return mtx

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
            
