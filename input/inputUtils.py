import numpy as nm

##
# 18.10.2006, c
# last revision: 04.06.2007
def matchGridLine( coor1, coor2, which ):
    if coor1.shape != coor2.shape:
        raise ValueError, 'incompatible shapes: %s == %s'\
              % ( coor1.shape, coor2.shape)

    c1 = coor1[:,which]
    c2 = coor2[:,which]
    i1 = nm.argsort( c1 )
    i2 = nm.argsort( c2 )

    if not nm.all( nm.abs(c1[i1] - c2[i2]) < 1e-12 ):
        print 'cannot match nodes\n(%s,\n %s), %e' % \
              (c1[i1], c2[i2], nm.abs(c1[i1] - c2[i2]).max())
        raise ValueError

    return i1, i2

##
# 18.10.2006, c
# last revision: 18.10.2006
def matchXLine( coor1, coor2 ):
    return matchGridLine( coor1, coor2, 0 )
def matchYLine( coor1, coor2 ):
    return matchGridLine( coor1, coor2, 1 )
def matchZLine( coor1, coor2 ):
    return matchGridLine( coor1, coor2, 2 )

##
# 01.06.2007, c
# last revision: 01.06.2007
def matchGridPlane( coor1, coor2, which ):
    from sfe.fem.mesh import findMap
    
    if coor1.shape != coor2.shape:
        raise ValueError, 'incompatible shapes: %s == %s'\
              % ( coor1.shape, coor2.shape)

    offset = coor1[0,which] - coor2[0,which]
    aux = coor2.copy()
    aux[:,which] += offset
    i1, i2 = findMap( coor1, aux, join = False )

    if i1.shape[0] != coor1.shape[0]:
        print 'cannot match nodes\n(%s,\n %s)' % (coor1[i1], coor2[i2])
        print i1
        print coor1
        print i2
        print coor2
        raise ValueError

    return i1, i2

##
# 01.06.2007, c
# last revision: 01.06.2007
def matchXPlane( coor1, coor2 ):
    return matchGridPlane( coor1, coor2, 0 )
def matchYPlane( coor1, coor2 ):
    return matchGridPlane( coor1, coor2, 1 )
def matchZPlane( coor1, coor2 ):
    return matchGridPlane( coor1, coor2, 2 )
