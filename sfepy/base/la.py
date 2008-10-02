from base import *
from sfepy.solvers import Solver

##
# Intersection of 1D arrays with unique elements.
# 01.11.2005, c
def intersect1d( array1, array2 ):
    aux = nm.sort( nm.concatenate( (array1, array2 ) ) )
    return nm.compress( (aux[1:] - aux[:-1]) == 0, aux )

##
# Intersection of 1D arrays with any elements.
# 01.11.2005, c
def intersect1d_nu( array1, array2 ):
    aux = nm.sort( nm.concatenate( (unique1d( array1 ), unique1d( array2  )) ) )
    return nm.compress( (aux[1:] - aux[:-1]) == 0, aux )

##
# Unique elements of 1D array.
# 01.11.2005, c
# 02.11.2005
# 16.12.2005
def unique1d( array1, ret_indx = False ):
    if len( array1 ) == 0:
        if ret_indx:
            return array1.copy(), array1.copy()
        else:
            return array1.copy()
            
    ar = nm.array( array1 ).flat
    if ret_indx:
        perm = nm.argsort( ar )
        aux = nm.take( ar, perm )
        ic = nm.empty( aux.shape, dtype = aux.dtype )
        ic[-1] = 1
        ic[:-1] = (aux[1:] - aux[:-1])
        flag = ic != 0
        return nm.compress( flag, perm ), nm.compress( flag, aux )
    else:
        aux = nm.sort( ar )
        ic = nm.empty( aux.shape, dtype = aux.dtype )
        ic[-1] = 1
        ic[:-1] = (aux[1:] - aux[:-1])
        return nm.compress( ic != 0, aux )
        
##
# 18.02.2005, c
# 21.02.2005
# 22.02.2005
def unique( array_in, mode = 'flat' ):
    from sfepy.fem.extmods.meshutils import sort_rows

    if mode == 'flat':
        aux = nm.sort( array_in.flat )
        ic = nm.zeros( aux.shape, aux.dtype )
        ic[-1] = 1
        ic[:-1] = nm.where( diff( aux ), 1, 0 )
    elif mode == 'rows':
        aux = array_in.copy()
        sort_rows( aux, nm.arange( aux.shape[1], dtype = nm.int32 ) )
        ic = nm.zeros( aux.shape[0], aux.dtype )
        ic[-1] = 1
        ic[:-1] = nm.where( nm.sum( nm.abs( \
                aux[1:,:] - aux[:-1,:] ), 1 ), 1, 0 )
    else:
        print 'unknown unique mode: %s' % mode
        raise ValueError

    array_out = aux[nm.where( ic )[0]]

    return( array_out )

##
# 14.01.2005, c
def as_unique_set( obj ):
    obj_s = nm.zeros( (0,), nm.int32 )
    for ii in obj:
        obj_s = nm.concatenate( (obj_s, ii.flat) )
    obj_s = nm.sort( obj_s )
#    obj_s = nm.where( obj_s == -1, , obj_s )
    flag = nm.zeros( obj_s.shape, nm.int32 )
    flag[obj_s] = 1;
    set = flag.nonzero()[0]
    return( set )

##
# 01.04.2005, c
# 04.09.2006
def rect( array, ir, ic ):
    ind2, ind1 = nm.meshgrid( ic, ir )
    return( array[ind1,ind2] )

def diff( obj ):
    return( obj[1:] - obj[:-1] )

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
# c: 25.09.2007, r: 08.04.2008
def eig( mtx_a, mtx_b = None, num = None, eigenvectors = True,
         return_time = None, method = 'eig.scipy', **ckwargs ):

    kwargs = {'name' : 'aux', 'kind' : method}
    kwargs.update( ckwargs )
    conf = Struct( **kwargs )
    solver = Solver.any_from_conf( conf )

    status = {}
    out = solver( mtx_a, mtx_b, num, eigenvectors, status )
    if return_time is not None:
        return_time[0] = status['time']
        
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
        out = [0] * nb
        for ii in xrange( bounds[0] ):
            for perm in cycle( bounds[1:] ):
                yield [ii] + perm


##
# 01.09.2007, c
def barycentric_coors( coors, s_coors ):
    n_v, dim = s_coors.shape
    n_c, dim2 = coors.shape
    assert_( dim == dim2 )
    assert_( ((dim + 1) * dim / 2) == n_v )

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
            
