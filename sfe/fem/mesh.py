from sfe.base.base import *
import sfe.base.la as la
from sfe.base.ioutils import writeMesh, readMeshHDF5
from meshio import MeshIO

import os.path as op

##
# connsIn must be sorted by matId within a group!
# 16.06.2005, c
# 04.08.2005
def mesh_splitByMatId( connsIn, matIdsIn, descsIn ):

    conns = []
    matIds = []
    descs = []

    for ig, conn in enumerate( connsIn ):
        one = nm.array( [-1], nm.int32 )
        ii = la.diff( nm.concatenate( (one, matIdsIn[ig], one) ) ).nonzero()[0]
        nGr = len( ii ) - 1;
#        print ii, nGr
        for igr in range( 0, nGr ):
            conns.append( conn[ii[igr]:ii[igr+1],:].copy() )
            matIds.append( matIdsIn[ig][ii[igr]:ii[igr+1]] )
            descs.append( descsIn[ig] )
            
    return (conns, matIds, descs)

##
# Join groups of the same type.
# 03.10.2005, c
# 04.10.2005
# 22.02.2007
def mesh_joinGroups( conns, descs ):

    el = dictFromKeysInit( descs, list )
    for ig, desc in enumerate( descs ):
        el[desc].append( ig )
    groups = [ii for ii in el.values() if ii]
##     print el, groups

    descsOut, connsOut = [], []
    for group in groups:
        nEP = conns[group[0]].shape[1]
        conn = nm.zeros( (0, nEP), nm.int32 )
        for ig in group:
            conn = nm.concatenate( (conn, conns[ig]) )
        connsOut.append( conn )
        descsOut.append( descs[group[0]] )

    return connsOut, descsOut

##
# 28.05.2007, c
def makePointCells( indx, dim ):
    conn = nm.zeros( (indx.shape[0], dim + 1), dtype = nm.int32 )
    for ii in range( 0, dim + 1 ):
        conn[:,ii] = indx
    return conn

##
# 23.05.2007, updated from matlab version
# 24.05.2007
# 01.06.2007
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
    xd = nm.sum( nm.diff( xs, axis = 0 )**2.0, axis = 1 )

    ii = nm.where( xd < eps )[0]
    off1, off2 = ir[iis][ii], ir[iis][ii+1]
    i1, i2 = iis[ii] - off1, iis[ii+1] - off2

    dns = nm.where( off1 == off2 )[0]
    if dns.size:
        print 'double node(s) in:'
        for dn in dns:
            if off1[dn] == 0:
                print 'x1: %s %s' % (i1[dn], i2[dn])
            else:
                print 'x2: %s %s' % (i1[dn], i2[dn])
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
    # 25.01.2006, c
    # 17.02.2006
    # 04.08.2006
    # 19.02.2007
    # 29.08.2007
    def fromFile( fileName ):
        if isinstance( fileName, file ):
            trunk = 'fromDescriptor'
        else:
            trunk = op.splitext( fileName )[0]

        mesh = Mesh( trunk )

        tt = time.clock() 
        mesh.read( fileName )
        print "mesh.read = ", time.clock() - tt
        mesh._setShapeInfo()

        tt = time.clock() 
        mesh.split( 'matId' )
        print "mesh.split = ", time.clock() - tt
        mesh._setShapeInfo()
        return mesh
    fromFile = staticmethod( fromFile )

    ##
    # 17.02.2006, c
    # 04.08.2006
    # 03.10.2006
    # 30.03.2007
    # 28.05.2007
    # 05.06.2007
    def fromRegion( region, meshIn, ed = None, fa = None ):
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

        mesh.descs.append( {2 : '2_3', 3 : '3_4'}[meshIn.dim] )
        mesh.matIds.append( -nm.ones_like( region.allVertices ) )
        mesh.conns.append( makePointCells( region.allVertices, meshIn.dim ) )

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
    # 26.09.2006, c
    # 29.09.2006
    # 21.02.2007
    def fromFileHDF5( fileName ):
        name, nod0, conns, matIds, descs = readMeshHDF5( fileName )
        return Mesh.fromData( name, nod0, conns, matIds, descs )
    fromFileHDF5 = staticmethod( fromFileHDF5 )

    ##
    # 21.02.2007, c
    def fromData( name, coors, conns, matIds, descs ):
        mesh = Mesh( name = name,
                     nod0 = coors,
                     conns = conns,
                     matIds = matIds,
                     descs = descs )
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
    # c: 17.01.2005, r: 05.02.2008
    def read( self, fileName = None, io = None ):
        """passing *MeshIO instance has precedence over fileName"""
        if io is None:
            if fileName is None:
                output( 'fileName or io must be specified!' )
                raise ValueError
            else:
                io = MeshIO.anyFromFileName( fileName )
        self.io = io
        self.nod0, connsIn, descs = io.read()

        # Detect wedges and pyramides -> separate groups.
        if ('3_8' in descs):
            ic = descs.index( '3_8' )

            connIn = connsIn.pop( ic )
            flag = nm.zeros( (connIn.shape[0],), nm.int32 )
            for ii, el in enumerate( connIn ):
                if (el[4] == el[5]):
                    if (el[5] == el[6]):
                        flag[ii] = 2
                    else:
                        flag[ii] = 1

            conn = []
            desc = []

            ib = nm.where( flag == 0 )[0]
            if (len( ib ) > 0):
                conn.append( connIn[ib] )
                desc.append( '3_8' )

            iw = nm.where( flag == 1 )[0]
            if (len( iw ) > 0):
                ar = nm.array( [0,1,2,3,4,6,8], nm.int32 )
                conn.append( la.rect( connIn, iw, ar ) )
                desc.append( '3_6' )

            ip = nm.where( flag == 2 )[0]
            if (len( ip ) > 0):
                ar = nm.array( [0,1,2,3,4,8], nm.int32 )
                conn.append( la.rect( connIn, ip, ar ) )
                desc.append( '3_5' )

##             print "brick split:", ic, ":", ib, iw, ip, desc

            connsIn[ic:ic] = conn
            del( descs[ic] )
            descs[ic:ic] = desc

        # Sort by matId within a group, preserve order.
        conns = []
        matIds = []
        for ig, conn in enumerate( connsIn ):
            ii = nm.argsort( conn[:,-1], kind = 'mergesort' )
            conn = conn[ii]
            
            conns.append( conn[:,:-1].copy() )
            matIds.append( conn[:,-1].copy() )


        self.conns = conns
        self.matIds = matIds
        self.descs = descs

    ##
    # 16.06.2005, c
    # 04.08.2005
    def split( self, on, inPlace = 1 ):

        if on == 'matId':
            conns, matIds, descs = \
                   mesh_splitByMatId( self.conns, self.matIds, self.descs )
        else:
            print 'unknown splitting mode: %s' % on
            raise ValueError
        
        if inPlace:
            self.conns = conns
            self.matIds = matIds
            self.descs = descs
            return None
        else:
            obj = deepcopy( self )
            obj.conns = conns
            obj.matIds = matIds
            obj.descs = descs
            return obj

    ##
    # 23.01.2006, c
    # 17.02.2006
    def write( self, fileName = None, fd = None, nod = None, igs = None ):

        if fileName is None:
            fileName = self.name + '.mesh'
            
        if fd is None:
            try:
                fd_ = open( fileName, "w" );
            except:
                print "Cannot open " + fileName + " for writing!";
                raise "ERR_FileOpen"

        if nod is None:
            nod = self.nod0

        if igs is None:
            igs = range( len( self.conns ) )

        conns = [nm.concatenate( (self.conns[ig],
                                  self.matIds[ig][:,nm.newaxis]), 1 )
                 for ig in igs]
        descs = [self.descs[ig] for ig in igs]
        aux1, aux2 = mesh_joinGroups( conns, descs )
        writeMesh( fd_, nod, aux1, aux2 )

        if fd is None:
            fd_.close()

    ##
    # 23.05.2007, c
    def getBoundingBox( self ):
        return nm.array( [nm.amin( self.nod0[:,:-1], 0 ),
                          nm.amax( self.nod0[:,:-1], 0 )] )


    ##
    # c: 02.01.2008, r: 02.01.2008
    def localize( self, inod ):
        """Strips nodes not in inod and remaps connectivities."""
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
