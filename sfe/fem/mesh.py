from sfe.base.base import *
import sfe.base.la as la
from meshio import MeshIO

import os.path as op

##
# 28.05.2007, c
def makePointCells( indx, dim ):
    conn = nm.zeros( (indx.shape[0], dim + 1), dtype = nm.int32 )
    for ii in range( 0, dim + 1 ):
        conn[:,ii] = indx
    return conn

##
# 23.05.2007, updated from matlab version, r: 05.05.2008
def findMap( x1, x2, eps = 1e-8, allowDouble = False, join = True ):
    """
    Find a mapping between common coordinates in x1 and x2, such that
    x1[cmap[:,0]] == x2[cmap[:,1]]
    """
    off, dim = x1.shape
    ir = nm.zeros( (off + x2.shape[0],), dtype = nm.int32 )
    ir[off:] = off

    x1 = nm.round( x1.T / eps ) * eps
    x2 = nm.round( x2.T / eps ) * eps
    xx = nm.c_[x1, x2]

    keys = [xx[ii] for ii in range( dim )]
    iis = nm.lexsort( keys = keys )

    xs = xx.T[iis]
##     import scipy.io as io
##     io.write_array( 'sss1', x1.T )
##     io.write_array( 'sss2', x2.T )
##     io.write_array( 'sss', xs, precision = 16 )
##     pause()
    xd = nm.sqrt( nm.sum( nm.diff( xs, axis = 0 )**2.0, axis = 1 ) )

    ii = nm.where( xd < eps )[0]
    off1, off2 = ir[iis][ii], ir[iis][ii+1]
    i1, i2 = iis[ii] - off1, iis[ii+1] - off2
    dns = nm.where( off1 == off2 )[0]
    if dns.size:
        print 'double node(s) in:'
        for dn in dns:
            if off1[dn] == 0:
                print 'x1: %d %d -> %s %s' % (i1[dn], i2[dn],
                                              x1[:,i1[dn]], x1[:,i2[dn]])
            else:
                print 'x2: %d %d -> %s %s' % (i1[dn], i2[dn],
                                              x2[:,i1[dn]], x2[:,i2[dn]])
        if not allowDouble:
            raise ValueError

    if join:
        cmap = nm.c_[i1, i2]
        return cmap
    else:
        return i1, i2

##
# 23.05.2007, updated from matlab version
# 24.05.2007
# 25.05.2007
# 28.05.2007
def mergeMesh( x1, conns1, x2, conns2, cmap, eps = 1e-8 ):
    nc = cmap.shape[0]
    n1 = x1.shape[0]
    n2 = x2.shape[0]

    err = nm.sum( nm.sum( nm.abs( x1[cmap[:,0],:-1] - x2[cmap[:,1],:-1] ) ) )
    if abs( err ) > (10.0 * eps):
        print 'nonmatching meshes!', err
        raise ValueError

    mask = nm.ones( (n2,), dtype = nm.int32 )
    mask[cmap[:,1]] = 0
#    print mask, nm.cumsum( mask )
    remap = nm.cumsum( mask ) + n1 - 1
    remap[cmap[:,1]] = cmap[:,0]
#    print remap

    i2 = nm.setdiff1d( nm.arange(  n2, dtype = nm.int32 ),
                       cmap[:,1] )
    xx = nm.r_[x1, x2[i2]]

    conns = []
    for ii in xrange( len( conns1 ) ):
        conn = nm.vstack( (conns1[ii], remap[conns2[ii]]) )
        conns.append( conn )
    
    return xx, conns

##
# 24.05.2007, c
def makeMesh( coor, conns, meshIn ):
    nod0 = nm.c_[coor, nm.zeros( (coor.shape[0],), dtype = nm.float64 )]

    matIds = []
    for ii, conn in enumerate( conns ):
        matId = nm.empty( (conn.shape[0],), dtype = nm.int32 )
        matId.fill( meshIn.matIds[ii][0] )
        matIds.append( matId )
        
    meshOut = Mesh.fromData( 'merged mesh', nod0, conns,
                             matIds, meshIn.descs )
    return meshOut

##
# 30.08.2007, c
# 31.08.2007
# 05.09.2007
def findRefinement( coor, conns, fcoor, fconns, eps, checkRefined = True ):
    ##
    # Find mesh vertices in the refined mesh.
    vmap = findMap( coor, fcoor )
#    print vmap
    if vmap.shape[0] != coor.shape[0]:
        print 'nonconforming meshes!'
        print vmap.shape, coor.shape, fcoor.shape
        raise ValueError
    ii = nm.argsort( vmap[:,0] )
    vmap = vmap[ii,1]

    ##
    # Inverted connectivity for the refined mesh (= eonlist in mesh_graph()!).
    # Elements are numbered locally per group.
    print fcoor.shape[0]
    iconns = []
    for ig, conn in enumerate( fconns ):
        iconn = [[] for ii in xrange( fcoor.shape[0] )]
        for iel, row in enumerate( conn ):
            for node in row:
                iconn[node].append( iel )
        iconns.append( iconn )
##    print iconn

    ##
    # Create element -> fine elements map.
    dim = coor.shape[1]
    emaps = []
    for ig, conn in enumerate( conns ):
        emap = []
        iconn = iconns[ig]
        fconn = fconns[ig]

        nEP = fconn.shape[1]
        for iel, row in enumerate( conn ):
            fels = []
#            print '*** iel:', iel, row

            seeds = vmap[row]
            while seeds.size:
#                print 'seeds:', seeds

                newSeeds = []
                for fnode in seeds:
#                    print fnode
                    pfels = [fel for fel in iconn[fnode] if not fel in fels]
                    if not pfels:
                        continue
                    pfels = nm.array( pfels )
#                    print pfels
                    fnodes = nm.array( fconn[pfels,:] )
#                    print fnodes
##                 print fcoor[fnodes,:]

                    # Use average element coordinate (simplex -> always convex).
                    eflag = la.pointsInSimplex( fcoor[fnodes,:].sum( 1 ) / nEP,
                                                coor[row,:], eps )
#                    print eflag
                    if eflag.any():
                        fels.append( pfels[eflag] )
                        newSeeds.append( fnodes[eflag] )
#                    print pfels[eflag], '->', fels

#                print 'new:', newSeeds
                seeds = nm.setdiff1d( nm.asarray( newSeeds ).ravel(), seeds )
#                pause()
            emap.append( nm.array( fels ).ravel() )
##             print '->>', iel, fels
##             pause()
##         print '->>>', emap
        if checkRefined:
            nRef = 2**dim
            for row in emap:
                assert len( row ) == nRef
        emaps.append( emap )
##         pause()

    iemaps = []
    nFEls = []
    for ig, emap in enumerate( emaps ):
        iemap = nm.empty( (fconns[ig].shape[0],), dtype = nm.int32 )
        iemap.fill( -1 )
        nFEl = nm.empty( (len( emap ),), dtype = nm.int32 )
        for iel, row in enumerate( emap ):
            for fiel in row:
                iemap[fiel] = iel
            nFEl[iel] = len( row )
        iemaps.append( iemap )
        nFEls.append( nFEl )

    return emaps, iemaps, nFEls

##
# Mesh.
# 13.12.2004, c
# 02.01.2005
class Mesh( Struct ):
    """
    Contains the FEM mesh together with all utilities related to it.

    Input and output is handled by the MeshIO class and subclasses.
    The Mesh class only contains the real mesh - nodes, connectivity,
    regions, plus methods for doing operations on this mesh.

    Example of creating and working with a mesh:

    >>> from sfe.fem.mesh import Mesh
    >>> m = Mesh.fromFile("database/simple.vtk")
    sfe: reading mesh (database/simple.vtk)...
    sfe: ...done in 0.06 s
    >>> m.nod0
    array([[  1.00000000e-01,   2.00000000e-02,  -1.22460635e-18,
              0.00000000e+00],
           [  1.00000000e-01,   1.80193774e-02,   8.67767478e-03,
              0.00000000e+00],
           [  1.00000000e-01,   1.24697960e-02,   1.56366296e-02,
              0.00000000e+00],
           ..., 
           [  8.00298527e-02,   5.21598617e-03,  -9.77772215e-05,
              0.00000000e+00],
           [  7.02544004e-02,   3.61610291e-04,  -1.16903153e-04,
              0.00000000e+00],
           [  3.19633596e-02,  -1.00335972e-02,   9.60460305e-03,
              0.00000000e+00]])
    >>> m.conns
    [array([[ 28,  60,  45,  29],
           [ 28,  60,  57,  45],
           [ 28,  57,  27,  45],
           ..., 
           [353, 343, 260, 296],
           [353, 139, 181, 140],
           [353, 295, 139, 140]])]
    >>> m.matIds
    [array([6, 6, 6, ..., 6, 6, 6])]
    >>> m.descs
    ['3_4']
    >>> m
    Mesh:database/simple
    >>> print m
    Mesh:database/simple
      nEPs:
        [4]
      dim:
        3
      nEl:
        1348
      name:
        database/simple
      descs:
        ['3_4']
      nNod:
        354
      matIds:
        [array([6, 6, 6, ..., 6, 6, 6])]
      nEls:
        [1348]
      nod0:
        [[  1.00000000e-01   2.00000000e-02  -1.22460635e-18   0.00000000e+00]
         [  1.00000000e-01   1.80193774e-02   8.67767478e-03   0.00000000e+00]
         [  1.00000000e-01   1.24697960e-02   1.56366296e-02   0.00000000e+00]
         ..., 
         [  8.00298527e-02   5.21598617e-03  -9.77772215e-05   0.00000000e+00]
         [  7.02544004e-02   3.61610291e-04  -1.16903153e-04   0.00000000e+00]
         [  3.19633596e-02  -1.00335972e-02   9.60460305e-03   0.00000000e+00]]
      io:
        None
      conns:
        [array([[ 28,  60,  45,  29],
               [ 28,  60,  57,  45],
               [ 28,  57,  27,  45],
               ..., 
               [353, 343, 260, 296],
               [353, 139, 181, 140],
               [353, 295, 139, 140]])]
      setupDone:
        0
      elOffsets:
        [   0 1348]


    The Mesh().nod0 is an array of nodes and Mesh().conns is the list of
    elements of each type (see Mesh().desc), so for example if you want to know
    the coordinates of the nodes of the fifth finite element of the type 3_4 do:

    In [1]: a.descs
    Out[1]: ['3_4']

    So now you know that the finite elements of the type 3_4 are in a.conns[0]:

    In [2]: [a.nod0[n] for n in a.conns[0][4]]
    Out[2]:
    [array([ 3.0877856 , -4.40913864, -1.58148163,  0.        ]),
     array([ 3.28954489, -4.53265378, -2.07926241,  0.        ]),
     array([ 3.00343981, -4.09445003, -2.14632505,  0.        ]),
     array([ 3.51217117, -4.05946689, -1.68843294,  0.        ])]



    """
    ##
    # 19.01.2005, c
    # 03.03.2005
    # 08.03.2005
    # 05.10.2005
    # 04.08.2006
    # 29.08.2007
    def fromSurface( surfFaces, meshIn ):

        inod = la.asUniqueSet( surfFaces )
        nNod = len( inod )
        nNodM, nCol = meshIn.nod0.shape

        aux = nm.arange( nNod )
        remap = nm.zeros( (nNodM,), nm.int32 )
        remap[inod] = aux

        mesh = Mesh( meshIn.name + "_surf" )

        mesh.nod0 = meshIn.nod0[inod,:nCol]

        sfm = {3 : "2_3", 4 : "2_4"}
        mesh.conns = []
        mesh.descs = []
        mesh.matIds = []
        for ii, sf in enumerate( surfFaces ):
            nEl, nFP = sf.shape

            conn = remap[sf]
            matId = nm.empty( (conn.shape[0],), dtype = nm.int32 )
            matId.fill( ii )

            mesh.descs.append( sfm[nFP] )
            mesh.conns.append( conn )
            mesh.matIds.append( matId )

        mesh._setShapeInfo()
        
        return mesh
    fromSurface = staticmethod( fromSurface )

    ##
    # c: 25.01.2006, r: 23.06.2008
    def fromFile( fileName = None, io = 'auto' ):
        """passing *MeshIO instance has precedence over fileName"""
        if io == 'auto':
            if fileName is None:
                output( 'fileName or io must be specified!' )
                raise ValueError
            else:
                io = MeshIO.anyFromFileName( fileName )
                if isinstance( fileName, file ):
                    trunk = 'fromDescriptor'
                else:
                    trunk = op.splitext( fileName )[0]
        else:
            trunk = io.fileName

        output( 'reading mesh (%s)...' % (fileName) )
        tt = time.clock()
        mesh = Mesh( trunk )
        mesh = io.read( mesh )
        output( '...done in %.2f s' % (time.clock() - tt) )
        mesh._setShapeInfo()
        return mesh
    fromFile = staticmethod( fromFile )

    ##
    # c: 17.02.2006, r: 28.04.2008
    def fromRegion( region, meshIn, ed = None, fa = None, localize = None ):
        mesh = Mesh( meshIn.name + "_reg" )
        mesh.nod0 = meshIn.nod0.copy()
        
        mesh.conns = []
        mesh.descs = []
        mesh.matIds = []
        if region.hasCells():
            for ig in region.igs:
                mesh.descs.append( meshIn.descs[ig] )
                els = region.getCells( ig )
                mesh.matIds.append( meshIn.matIds[ig][els,:].copy() )
                mesh.conns.append( meshIn.conns[ig][els,:].copy() )

        if ed is not None:
            for ig in region.igs:
                edges = region.getEdges( ig )
                mesh.descs.append( '1_2' )
                mesh.matIds.append( ed.data[edges,0] + 1 )
                mesh.conns.append( ed.data[edges,-2:].copy() )

        if fa is not None:
            for ig in region.igs:
                faces = region.getFaces( ig )
                fdata = fa.data[faces]
                i3 = nm.where( fdata[:,-1] == -1 )[0]
                i4 = nm.where( fdata[:,-1] != -1 )[0]
                if i3.size:
                    mesh.descs.append( '2_3' )
                    mesh.matIds.append( fdata[i3,0] + 1 )
                    mesh.conns.append( fdata[i3,-4:-1].copy() )
                if i4.size:
                    mesh.descs.append( '2_4' )
                    mesh.matIds.append( fdata[i4,0] + 1 )
                    mesh.conns.append( fdata[i4,-4:].copy() )

        if (ed is not None) or (fa is not None):
            mesh.descs.append( {2 : '2_3', 3 : '3_4'}[meshIn.dim] )
            mesh.matIds.append( -nm.ones_like( region.allVertices ) )
            mesh.conns.append( makePointCells( region.allVertices, meshIn.dim ) )

        if localize:
            mesh.localize( region.allVertices )
        mesh._setShapeInfo()
        return mesh
    fromRegion = staticmethod( fromRegion )

    ##
    # c: 02.01.2008, r: 02.01.2008
    def fromRegionAndField( region, field ):
        mesh, ed, fa = field.domain.mesh, field.domain.ed, field.domain.fa
        mesh = Mesh.fromRegion( region, mesh, ed, fa )
        mesh.name = mesh.name + '_field'

        nodes = region.getFieldNodes( field, merge = True )

        aux = field.getExtraNodesAsSimplices( nodes )
        mesh.nod0 = field.aps.coors
        mesh.descs.append( aux[0] )
        mesh.matIds.append( aux[1] )
        mesh.conns.append( aux[2] )

        mesh.localize( nodes )
        mesh._setShapeInfo()
        return mesh
    fromRegionAndField = staticmethod( fromRegionAndField )

    ##
    # c: 21.02.2007, r: 08.02.2008
    def fromData( name, coors, conns, matIds, descs, igs = None ):
        if igs is None:
            igs = range( len( conns ) )
        mesh = Mesh( name = name,
                     nod0 = coors,
                     conns = [conns[ig] for ig in igs],
                     matIds = [matIds[ig] for ig in igs],
                     descs = [descs[ig] for ig in igs] )
        mesh._setShapeInfo()
        return mesh
    fromData = staticmethod( fromData )
        

    ##
    # 22.02.2005
    # 16.06.2005
    # 26.09.2006
    def __init__( self, name = 'mesh', **kwargs ):
        Struct.__init__( self, **kwargs )
        self.name = name
        self.io = None
        self.setupDone = 0

    ##
    # 04.08.2006, c
    # 29.09.2006
    def _setShapeInfo( self ):
        self.nNod = self.nod0.shape[0]
        self.dim = self.nod0.shape[1] - 1
        self.nEls = nm.array( [conn.shape[0] for conn in self.conns] )
        self.nEPs = nm.array( [conn.shape[1] for conn in self.conns] )
        self.elOffsets = nm.cumsum( nm.r_[0, self.nEls] )
        self.nEl = nm.sum( self.nEls )

    ##
    # c: 15.02.2008, r: 15.02.2008
    def _setData( self, coors, conns, matIds, descs ):
        self.nod0 = coors
        self.conns = conns
        self.matIds = matIds
        self.descs = descs
        
    ##
    # c: 23.01.2006, r: 23.06.2008
    def write( self, fileName = None, io = None,
               coors = None, igs = None, out = None, **kwargs ):
        """Write mesh + optional results in 'out'.

        'io' == 'auto' respects the extension of 'fileName'
        'coors' can be used instead of mesh coordinates,
        providing 'igs' filters some groups only"""
        if fileName is None:
            fileName = self.name + '.mesh'

        if io is None:
            io = self.io
        else:
            if io == 'auto':
                io = MeshIO.anyFromFileName( fileName )

        if coors is None:
            coors = self.nod0

        if igs is None:
            igs = range( len( self.conns ) )

        auxMesh = Mesh.fromData( self.name, coors,
                                 self.conns, self.matIds, self.descs, igs )
        io.write( fileName, auxMesh, out, **kwargs )

    ##
    # 23.05.2007, c
    def getBoundingBox( self ):
        return nm.array( [nm.amin( self.nod0[:,:-1], 0 ),
                          nm.amax( self.nod0[:,:-1], 0 )] )


    ##
    # c: 02.01.2008, r: 02.01.2008
    def localize( self, inod ):
        """Strips nodes not in inod and remaps connectivities.
        TODO: fix the case when remap[conn] contains -1..."""
        remap = nm.empty( (self.nod0.shape[0],), dtype = nm.int32 )
        remap.fill( -1 )
        remap[inod] = nm.arange( inod.shape[0], dtype = nm.int32 )

        self.nod0 = self.nod0[inod]
        conns = []
        for conn in self.conns:
            conns.append( remap[conn] )
        self.conns = conns


    ##
    # c: 18.01.2008, r: 18.01.2008
    def transformCoords( self, mtxT, refCoords = None ):
        """x = T * x."""
        if refCoords is None:
            refCoords = self.nod0[:,:-1]

        self.nod0[:,:-1] = nm.dot( refCoords, mtxT.T )
