from base import *
from sfepy.fem.extmods.meshutils import sortRows
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
def unique1d( array1, retIndx = False ):
    if len( array1 ) == 0:
        if retIndx:
            return array1.copy(), array1.copy()
        else:
            return array1.copy()
            
    ar = nm.array( array1 ).flat
    if retIndx:
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
def unique( arrayIn, mode = 'flat' ):

    if mode == 'flat':
        aux = nm.sort( arrayIn.flat )
        ic = nm.zeros( aux.shape, aux.dtype )
        ic[-1] = 1
        ic[:-1] = nm.where( diff( aux ), 1, 0 )
    elif mode == 'rows':
        aux = arrayIn.copy()
        sortRows( aux, nm.arange( aux.shape[1], dtype = nm.int32 ) )
        ic = nm.zeros( aux.shape[0], aux.dtype )
        ic[-1] = 1
        ic[:-1] = nm.where( nm.sum( nm.abs( \
                aux[1:,:] - aux[:-1,:] ), 1 ), 1, 0 )
    else:
        print 'unknown unique mode: %s' % mode
        raise ValueError

    arrayOut = aux[nm.where( ic )[0]]

    return( arrayOut )

##
# 14.01.2005, c
def asUniqueSet( obj ):
    objS = nm.zeros( (0,), nm.int32 )
    for ii in obj:
        objS = nm.concatenate( (objS, ii.flat) )
    objS = nm.sort( objS )
#    objS = nm.where( objS == -1, , objS )
    flag = nm.zeros( objS.shape, nm.int32 )
    flag[objS] = 1;
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
def splitRange( nItem, step ):
    num = nItem / step
    out = [step] * num
    aux = sum( out )
    if aux < nItem:
        out.append( nItem - aux )

    return out

##
# c: 25.09.2007, r: 08.04.2008
def eig( mtxA, mtxB = None, num = None, eigenvectors = True, returnTime = None,
         method = 'eig.scipy', **ckwargs ):

    kwargs = {'name' : 'aux', 'kind' : method}
    kwargs.update( ckwargs )
    conf = Struct( **kwargs )
    solver = Solver.anyFromConf( conf )

    status = {}
    out = solver( mtxA, mtxB, num, eigenvectors, status )
    if returnTime is not None:
        returnTime[0] = status['time']
        
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
def barycentricCoors( coors, sCoors ):
    nV, dim = sCoors.shape
    nC, dim2 = coors.shape
    assert dim == dim2
    assert ((dim + 1) * dim / 2) == nV

    mtx = nm.ones( (nV, nV), nm.float64 )
    mtx[0:nV-1,:] = sCoors.T
    rhs = nm.empty( (nV,nC), nm.float64 )
    rhs[0:nV-1,:] = coors.T
    rhs[nV-1,:] = 1.0
    bc = nla.solve( mtx, rhs )
##     print bc.T
    return bc

##
# 30.08.2007, c
# 01.09.2007
def pointsInSimplex( coors, sCoors, eps = 1e-8 ):
    nC, dim = coors.shape
    bc = barycentricCoors( coors, sCoors )
    flag = nm.ones( (nC,), dtype = nm.bool )
    for idim in xrange( dim + 1 ):
        flag &= nm.where( (bc[idim,:] > -eps)
                          & (bc[idim,:] < (1.0 + eps)), True, False )
    return flag

##
# c: 18.01.2008, r: 18.01.2008
def rotationMatrix2D( angle ):
    angle *= nm.pi / 180.0
    mtx = nm.array( [[nm.cos( angle ), -nm.sin( angle )],
                     [nm.sin( angle ), nm.cos( angle )]], dtype = nm.float64 )
    return mtx

##
# 30.08.2007, c
class MatrixAction( Struct ):

    ##
    # 30.08.2007, c
    def fromFunction( fun, expectedShape, dtype ):
        def call( self, vec ):
            aux = fun( vec )
            assert aux.shape[0] == self.shape[0] 
            return nm.asanyarray( aux, dtype = self.dtype )
        obj = MatrixAction( shape = expectedShape,
                            dtype = dtype,
                            kind = 'function' )
        insertMethod( obj, call )
        return obj
    fromFunction = staticmethod( fromFunction )

    ##
    # 30.08.2007, c
    def fromArray( arr ):
        def call( self, vec ):
            return nm.asarray( sc.dot( self.arr, vec ) )
        obj = MatrixAction( shape = arr.shape,
                            dtype = arr.dtype,
                            arr = arr,
                            kind = 'array' )
        insertMethod( obj, call )
        return obj
    fromArray = staticmethod( fromArray )
    
    ##
    # 30.08.2007, c
    def __call__( self, vec ):
        return self.call( vec )

    ##
    # 30.08.2007, c
    def toArray( self ):
        if self.kind == 'array':
            return self.arr
        else:
            print 'cannot make array from MatrixAction of kind %s!' % self.kind
            raise ValueError
            
