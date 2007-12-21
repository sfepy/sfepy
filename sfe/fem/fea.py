from sfe.base.base import *
from feGenerators import ap_nodeGenerators, ap_bfGenerators
from baseFunction import BaseFunction
from quadratures import collectQuadratures
import extmods.meshutils as mu
import extmods.geometry as gm

##
# created:       26.01.2006
# last revision: 21.12.2007
def setMeshCoors( domain, fields, geometries, coors, updateState = False ):
    domain.mesh.nod0[:,:-1] = coors
    if updateState:
        fields.setupCoors()
        for field in fields:
            field.aps.updateGeometry( field, domain.regions, geometries )


##
# 04.08.2005, c
def _interpToFaces( vertexVals, bfs, faces ):
    dim = vertexVals.shape[1]
    nFace = faces.shape[0]
    nQP = bfs.shape[0]
    
    facesVals = nm.zeros( (nFace, nQP, dim), nm.float64 )
    for ii, face in enumerate( faces ):
        vals = vertexVals[face,:dim]
        facesVals[ii,:,:] = nm.dot( bfs[:,0,:], vals )

    return( facesVals )

def _getIName( iname, integral, key ):
    if iname is None:
        if integral is None:
            print 'no integral name given for key "%s"' % key
            raise ValueError
        else:
            iname = integral.name
    return iname

##
# 02.08.2005, c
# 30.09.2005
# 03.10.2005
# 07.10.2005
# 19.12.2005
# 02.08.2006
# 11.04.2007
# 03.05.2007
# 04.05.2007
# 24.05.2007
def evalBF( qp, baseFun, nodes, derivative ):
    """bfB(G) indexing is (ifa,iqp,:,nEP) -> can use FMF_SetCell"""
    fun, varSet = baseFun.fun, baseFun.varSet
    if (qp.vals.ndim == 2):
        if derivative == 0:
            val = fun.value( qp.vals, nodes )
        else:
            val = fun.value( qp.vals, nodes, varSet )

    else: # Boundary QP.
        sh = qp.vals.shape[0:2]
        if derivative == 0:
            val = nm.zeros( sh + (1, len( nodes )), nm.float64 )
            for ifa, bqp in enumerate( qp.vals ):
                val[ifa,:,:,:] = fun.value( bqp, nodes )
        else:
            val = nm.zeros( sh + (len( varSet ), len( nodes )), nm.float64 )
            for ifa, bqp in enumerate( qp.vals ):
                val[ifa,:,:,:] = fun.value( bqp, nodes, varSet )


    if derivative == 0:
        aux = nm.sum( val, -1 ).ravel()
        if not nm.alltrue( nm.absolute( aux - 1.0 ) < 1e-14 ):
            print aux
            raise AssertionError

    return val

##
# Field approximation.
# 29.11.2004, c
# 30.11.2004
# 01.12.2004
# 04.02.2005
class Interpolant( Struct ):

    ##
    # 31.03.2005
    # 16.06.2005
    def __init__( self ):
        self.setupDone = 0

    ##
    # 09.03.2005, c
    # 31.03.2005
    # 22.06.2005
    # 19.07.2005
    # 20.07.2005
    # 02.08.2005
    # 18.07.2006
    # 23.08.2006
    # 11.04.2007
    def setup( self, gel = None ):
        if (self.setupDone): return

        if gel is not None:
            self.gel = gel

        # Input transformation and defaults.
        for key, nod in self.nodes.iteritems():
            if not nod.has_key( 'share' ):
                nod['share'] = 1

        self.desc = dictToStruct( self.desc, flag = (1,) )
        self.nodes, self.baseFuns\
                    = dictToStruct( self.nodes, self.baseFuns,
                                    flag = (0, 1) )
        name = self.desc.family
        self.iKeyMap = iKeyMap = invertDict( self.keyMap )
        # Nodes.
        for key, nod in self.nodes.iteritems():
            if nod.mode == 'generate':
                try:
                    gen = ap_nodeGenerators[name]
                except:
                    print self.desc.family
                    raise NotImplementedError, name
                if not key == 'v':
                    gen = gen[1]
                else:
                    gen = gen[0]
                aux = gen( self.desc.approxOrder, iKeyMap[key],
                           nod, self.gel, wantBarCoors = 1 )
                # ntx := node types, a row = (type, entity (vefb)) index
                nod.nts, nod.vals, nod.barCoors = aux
            else:
                print 'unknown mode', nod.mode
                raise NotImplementedError

        # Base functions.
        for key, bf in self.baseFuns.items():
            gd = self.gel.data[iKeyMap[key]]
            bf.varSet = range( gd.dim )
            
            if bf.mode == 'generate':
                try:
                    gen = ap_bfGenerators[name]
                except:
                    print self.desc.family
                    raise NotImplementedError, name
                if not key == 'v':
                    gen = gen[1]
                else:
                    gen = gen[0]
                bfgen = gen( gd.coors, bf, self.nodes[key].vals, bf.varSet )
                bf.fun = BaseFunction( bfgen, self.nodes[key].vals, bf.varSet )

            else:
                print 'unknown mode', bf.mode
                raise NotImplementedError


##         for key, nod in self.nodes.iteritems():
##             print key, nod
##         pause()
            
        self.setupDone = 1



    ##
    # 02.08.2005, c
    # 22.09.2005
    # 23.09.2005
    # 30.09.2005
    # 03.10.2005
    def listExtraNodeTypes( self, et, ft ):
        gd = self.gel.data['v']
        nodes = self.nodes['v']
        maxAO = nm.amax( nm.sum( nodes.vals, 1 ) )

##         print nodes
##         print 'sdsd', maxAO
##         pause()

        for ii, nt in enumerate( nodes.nts ):
            if (nt[0] == 1): # Edge node.
                edge = gd.edges[nt[1]]
                tpar = float( nodes.vals[ii,edge[1]] )
                key = int( round( tpar / maxAO, 5 ) * 1e5 )
                if not et.has_key( key ): et[key] = len( et )
##                 print ii, nt
##                 print tpar
##                 print nodes.vals[ii], edge

            elif (nt[0] == 2): # Face node.
                face = gd.faces[nt[1]]
                upar = float( nodes.vals[ii,face[1]] )
                vpar = float( nodes.vals[ii,face[-1]] )
                key = [int( round( ii / maxAO, 5 ) * 1e5 )
                       for ii in [upar, vpar]]
                key = tuple( key )
                if not ft.has_key( key ): ft[key] = len( ft )
##                 print upar, vpar
##                 print nodes.vals[ii], face

        return (et, ft)


    ##
    # 1. col ... iep, 2. col ... isExtra (simplified for now)
    # 03.10.2005, c
    def describeNodes( self ):

        def _describeNodes_aux( inds ):
            objs = []
            ii = 0
            while (ii < len( inds )):
                nt = nts[inds[ii]]
                first = nt[1]
                obj = []
                while (first == nt[1]):
                    obj.append( [inds[ii], 1] )
                    ii += 1
                    if (ii >= len( inds )):
                        break
                    nt = nts[inds[ii]]
                objs.append( nm.array( obj ) )

            return objs

        nts = self.nodes['v'].nts

        nodeDesc = Struct()

        # Vertex node.
        ii = (nts[:,0] == 0).nonzero()[0]
        zz = nm.zeros( (ii.shape[0], 1), dtype = nm.int32 )
        nodeDesc.vertex = nm.concatenate( (ii[:,nm.newaxis], zz), 1 )

        # Edge nodes.
        inds = (nts[:,0] == 1).nonzero()[0]
        nodeDesc.edge = _describeNodes_aux( inds )
        
        # Face nodes.
        inds = (nts[:,0] == 2).nonzero()[0]
        nodeDesc.face = _describeNodes_aux( inds )

        # Bubble nodes.
        inds = (nts[:,0] == 3).nonzero()[0]
        nodeDesc.bubble = _describeNodes_aux( inds )
        
        return nodeDesc

    ##
    # 16.11.2007, c
    def getNNodes( self ):
        nn = {}
        for key, nodes in self.nodes.iteritems():
            nn[key] = nodes.vals.shape[0]
        return nn

##
# 18.07.2006, c
class Approximation( Struct ):
    ##
    # 18.07.2006, c
    # 10.10.2006
    # 11.07.2007
    # 17.07.2007
    def __init__( self, name, interp, region, ig ):
        """interp, region are borrowed."""

        self.name = name
        self.interp = interp
        self.region = region
        self.ig = ig
        self.surfaceData = {}
        self.edgeData = {}
        self.pointData = {}
        self.qpCoors = {}
        self.bf = {}
        self.integrals = {}
        self.nEP = self.interp.getNNodes()

    ##
    # 19.07.2006, c
    def describeNodes( self, nodeDesc ):
        self.nodeDesc = nodeDesc

        self.hasExtraEdgeNodes = self.hasExtraFaceNodes = False
        if self.nodeDesc.edge:
            self.hasExtraEdgeNodes = True
        if self.nodeDesc.face:
            self.hasExtraFaceNodes = True

##         print self.nodeDesc
##        pause()

    ##
    # 03.10.2005, c
    # 04.10.2005
    # 10.10.2005
    # 26.10.2005
    # 19.07.2006
    # 02.08.2006
    # 04.08.2006
    # 13.02.2007
    # 20.02.2007
    def evalExtraCoor( self, coors, mesh ):
        """The last coors column (id) is not evaluated."""
        nodeOffsets = self.nodeOffsets
        nNod = nm.sum( nodeOffsets[1:,1] - nodeOffsets[1:,0] )
#        print self.name
#        print nNod
##         pause()

        if nNod == 0: return

        ##
        # Evaluate geometry interpolation base functions in extra nodes.
        ginterp = self.interp.gel.interp
        nodes = self.interp.nodes['v']

        iex = (nodes.nts[:,0] > 0).nonzero()[0]
        qpCoors = nodes.barCoors[iex,:]

        qp = Struct( vals = qpCoors )
        bf = evalBF( qp, ginterp.baseFuns['v'], ginterp.nodes['v'].vals, 0 )
        bf = bf[:,0,:].copy()

        ##
        # Evaulate extra coordinates with 'bf'.
#        print self.econn.shape, self.sub.nEl

        econn = nm.zeros_like( self.econn[:,iex] ).copy()
        for ii in range( 1, 4 ):
            ix = nm.where( nodes.nts[:,0] == ii )[0]
            if not len( ix ): continue
            
##             print ii, ix, iex[0], ix-iex[0]
##             pause()
            econn[:,ix-iex[0]] = self.econn[:,ix]

##         print self.econn
##         print econn
##         print coors.shape, nm.amax( econn )
##         pause()

        region = self.region
        group = region.domain.groups[self.ig]
        cells = region.getCells( self.ig )
        mu.interpVertexData( coors, econn, mesh.nod0, group.conn[cells], bf, 1 )

    ##
    # 24.07.2006, c
    def getVDataShape( self, iname = None ):
        """returns (nEl, nQP, dim, nEP)"""
        if iname:
            return (self.region.shape[self.ig].nCell,)\
                   + self.bf[(iname, 'v', 1)].shape
        else:
            return (self.region.shape[self.ig].nCell,
                    self.interp.gel.dim, self.nEP['v'])

    ##
    # 05.09.2006, c
    # 20.02.2007
    # 17.07.2007
    def getSDataShape( self, iname, key ):
        """returns (nFa, nQP, dim, nFP)"""
        if not self.surfaceData:
            return 0, 0, 0, 0
        sd = self.surfaceData[key]
        bfSG = self.bf[(iname, sd.faceType, 1)]
        nQP, dim, nFP = bfSG.shape
        assert nFP == sd.nFP
        
        return sd.nFa, nQP, dim + 1, nFP

    ##
    # 05.09.2006, c
    # 06.09.2006
    # 10.10.2006
    # 12.10.2006
    # 20.02.2007
    # 16.03.2007
    # 26.03.2007
    def setupSurfaceData( self, region ):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        faceIndices = region.fis[self.ig]
        
        faces = self.efaces[faceIndices[:,1]]
        if faces.size == 0:
            raise RuntimeError, 'region with group with no faces! (%s)'\
                  % region.name
#        print self.efaces
        try:
            ee = self.econn[faceIndices[:,0]]
        except:
            pdb.set_trace()

        econn = nm.empty( faces.shape, dtype = nm.int32 )
        for ir, face in enumerate( faces ):
#            print ir, face, ee[ir,face]
            econn[ir] = ee[ir,face]

        ef = econn.flat
        nodes = nm.unique1d( ef )

        aux = -nm.ones( (nm.max( ef ) + 1,), dtype = nm.int32 )
        aux[nodes] = nm.arange( len( nodes ) )
        leconn = aux[econn].copy()
        assert nm.alltrue( nodes[leconn] == econn )
        
        nFa, nFP = faceIndices.shape[0], faces.shape[1]
        faceType = 's%d' % nFP

        sd = Struct( name = 'surfaceData_%s' % region.name,
                     econn = econn, fis = faceIndices, nFa = nFa, nFP = nFP,
                     nodes = nodes, leconn = leconn, faceType = faceType )
        self.surfaceData[region.name] = sd
        return sd

    ##
    # 11.07.2007, c
    def setupPointData( self, field, region ):
        conn = region.getFieldNodes( field, merge = True, igs = region.igs )
##         conn = [nods]\
##                + [nm.empty( (0,), dtype = nm.int32 )]\
##                * (len( region.igs ) - 1)

        conn.shape += (1,)
        self.pointData[region.name] = conn

    ##
    # created:       29.11.2007
    # last revision: 11.12.2007
    def getQP( self, key, iname, integral = None ):
        """integrals are stored in self.integrals..."""
        qpkey = (iname, key)

        if not self.qpCoors.has_key( qpkey ):
            if integral is None:
                try:
                    integral = self.integrals[iname]
                except:
                    print self.name
                    print 'no integral given for key "%s"' % str( qpkey )
                    raise ValueError
            interp = self.interp
            if key[0] == 's':
                dim = interp.gel.dim - 1
                nFP = len( interp.gel.sEdges[interp.iKeyMap[key]][0] )
                geometry = '%d_%d' % (dim, nFP)
            else:
                geometry = interp.geometry
            vals, weights = integral.getQP( geometry )()
            self.qpCoors[qpkey] = Struct( vals = vals, weights = weights )
##             print self.name, self.qpCoors
##             pause()
        return self.qpCoors[qpkey]

    ##
    # created:       29.11.2007
    def getBase( self, key, derivative, iname = None, integral = None,
                 fromGeometry = False, baseOnly = True ):
        iname = _getIName( iname, integral, key )
        qp = self.getQP( key, iname, integral )

        if fromGeometry:
            gkey = self.interp.iKeyMap[key]
            baseFun = self.interp.gel.interp.baseFuns[gkey]
            nodes = self.interp.gel.interp.nodes[gkey].vals
            bfKey = (iname, 'g' + key, derivative)
        else:
            baseFun = self.interp.baseFuns[key]
            nodes = self.interp.nodes[key].vals
            bfKey = (iname, key, derivative)

        if not self.bf.has_key( bfKey ):
            self.bf[bfKey] = evalBF( qp, baseFun, nodes, derivative )            

        if baseOnly:
            return self.bf[bfKey]
        else:
            return self.bf[bfKey], qp.weights

    ##
    # 05.09.2006, c
    # 03.05.2007
    def describeGeometry( self, field, geomRequest, integral, coors ):
        gtype = geomRequest.gtype
        if gtype == 'Volume':
            bfVG, weights = self.getBase( 'v', 1, integral = integral,
                                          baseOnly = False )

##             print bfVG, weights
##             pause()
            vg = gm.VolumeGeometry( *self.getVDataShape( integral.name ) )
            vg.mode = gm.GM_Material
            try:
                vg.describe( coors, self.econn, bfVG, weights )
            except:
#                pdb.set_trace()
                gm.errclear()
                raise
            return vg

        elif (gtype == 'Surface') or (gtype == 'SurfaceExtra'):
##             print self.name
##             print geomRequest
            region = geomRequest.region
            fa = region.domain.getNeighbourLists( True )[2]
            region.setupFaceIndices( fa )
            region.selectCellsOfSurface()
            self.setupSurfaceData( region )
            
            sd = self.surfaceData[region.name]
##             print sd
##             print integral
            bfSG, weights = self.getBase( sd.faceType, 1, integral = integral,
                                          baseOnly = False )
            print bfSG, weights 
#            pause()

            sg = gm.SurfaceGeometry( *self.getSDataShape( integral.name,
                                                          region.name ) )
            sg.mode = gm.GM_Material

##             print sg
            
            try:
                sg.describe( coors, sd.econn, bfSG, weights )
            except:
                gm.errclear()
                raise

            if gtype == 'SurfaceExtra':
                sg.allocExtraData( self.getVDataShape()[2] )

                bfS = self.getBase( sd.faceType, 0, integral = integral,
                                    fromGeometry = True )
##                 print bfS
                bkey = self.createBQP( sd.faceType, bfS, weights, integral.name )
                bfBG = self.getBase( bkey, 1, integral = integral )
##                 print bfBG
                sg.evaluateBFBGM( bfBG, coors, sd.fis, self.econn )

##             sg.str( sys.stdout, 0 )
##             pause()
            
            return sg

        elif gtype == 'Point':
            region = geomRequest.region
##             print gtype, region.name, self.ig
##             pause()
            if self.ig == region.igs[0]:
                # Point data only in the first group to avoid multiple
                # assembling of nodes on group boundaries.
                self.setupPointData( field, region )

        else:
            print 'unknown geometry type:', gtype
            raise ValueError

    ##
    #
    def createBQP( self, skey, bfS, weights, iname ):
        interp = self.interp
        gd = interp.gel.data['v']
        bkey = 'b%s' % skey[1:]
        bqpkey = (iname, bkey)
        coors, faces = gd.coors, gd.faces

        vals = _interpToFaces( coors, bfS, faces )
        self.qpCoors[bqpkey] = Struct( name = 'BQP_%s' % bkey,
                                     vals = vals, weights = weights )
##         print self.qpCoors[bqpkey]
        interp.nodes[bkey] = interp.nodes['v']
        interp.baseFuns[bkey] = interp.baseFuns['v']
        return bkey

##
# 14.07.2006, c
class Approximations( Container ):
    """- Region can span over several groups -> different Aproximation
    instances

    - interps and hence nodeDescs are per region (must have single
    geometry!)"""

    ##
    # last revision: 21.12.2007
    def __init__( self, bases, interps, domain ):

        self.igs = []
        apsPerRegion = {}
        self.interps = {}
        for regionName, baseName in bases.iteritems():
            print regionName, baseName

            try:
                region = domain.regions[regionName]
            except:
                output( 'region %s does not exist' % regionName )
                raise

            interp = interps[baseName]
            if self.interps.has_key( region.name ):
                if self.interps[region.name] is not interp:
                    output( 'interpolation mismatch!' )
                    output( self.interps[region.name].name, interp.name )
                    raise ValueError
            else:
                interp.setup( domain.geomEls[interp.geometry] )
                self.interps[region.name] = interp
            
            apsPerGroup = {}
            objs = OneTypeList( Approximation )
            for ig in region.igs:
                if ig in self.igs:
                    output( 'regions must not intersect (%s, %d)'\
                            % (region.name, ig) )
                    raise ValueError

                self.igs.append( ig )
                
                ap = Approximation( baseName + '_%s_ig%d' % (region.name, ig),
                                    interp, region, ig )
                apsPerGroup[ig] = ap
                objs.append( ap )

            apsPerRegion[region.name] = apsPerGroup

        self.apsPerRegion = apsPerRegion
        self.update( objs )

        self.apsPerGroup = {}
        for regionName, ig, ap in self.iterAps( all = True ):
            self.apsPerGroup.setdefault( ig, {} )[regionName] = ap

        self.clearGeometries()

    ##
    # created:       21.12.2007
    # last revision: 21.12.2007
    def clearGeometries( self ):
        self.geometries = {}

    ##
    #
    def iterAps( self, all = False, igs = None ):
        if all:
            for regionName, apsPerGroup in self.apsPerRegion.iteritems():
                for ig, ap in apsPerGroup.iteritems():
                    if igs is not None:
                        if not ig in igs: continue
                    yield regionName, ig, ap
        else:
            for regionName, apsPerGroup in self.apsPerRegion.iteritems():
                yield regionName, apsPerGroup

    ##
    #
    def __len__( self ):
        nn = 0
        for apsPerGroup in self.apsPerRegion.itervalues():
            nn += len( apsPerGroup )
        return nn

    ##
    # 19.07.2006, c
    def describeNodes( self ):
        self.ent = ent = {}
        self.fnt = fnt = {}

        self.nodeDescs = {}
        for key, interp in self.interps.iteritems():
            interp.listExtraNodeTypes( ent, fnt )
            nodeDesc = interp.describeNodes()
            self.nodeDescs[key] = nodeDesc

##             print ent, fnt
##         pause()

        for key1, key2, ap in self.iterAps( all = True ):
            ap.describeNodes( self.nodeDescs[key1] )
            

    ##
    # 23.09.2005, c
    # 29.09.2005
    # 07.10.2005
    # 19.07.2006
    def setupNodes( self ):

        self.edgeOris = {}
        self.faceOris = {}
        for key1, key2, ap in self.iterAps( all = True ):
##             print ap
##             print key1, key2
            domain = ap.region.domain
            igs = ap.region.igs
            for ig in igs:
                if not self.edgeOris.has_key( ig ):
                    self.edgeOris[ig] = domain.getOrientation( ig )
                if not self.faceOris.has_key( ig ):
                    self.faceOris[ig] = None
#            pause()

        ##
        # Edge node type permutation table.
        nENT = len( self.ent )
        entt = nm.zeros( (2, nENT), nm.int32 )
        entt[0] = nm.arange( nENT )
        entt[1] = nm.arange( nENT - 1, -1, -1 )
        self.entTable = entt

##         print self.ent
##         print self.entTable
#        pause()

        ##
        # Face node type permutation table.
        
    ##
    # 30.09.2005, c
    # 03.10.2005
    # 04.10.2005
    # 07.10.2005
    # 10.10.2005
    # 01.11.2005
    # 03.12.2005
    # 16.12.2005
    # 19.07.2006
    # 23.08.2006
    # 10.10.2006
    # 13.02.2007
    # 19.02.2007
    # 20.02.2007
    def setupGlobalBase( self, domain ):
        """
        efaces: indices of faces into econn.
        """

        nodeOffsetTable = nm.zeros( (4, len( self ) + 1), dtype = nm.int32 )
##         print nodeOffsetTable.shape
        
        iVertex, iEdge, iFace, iBubble = 0, 1, 2, 3

        # Global node number.
        iseq = 0

        ##
        # Vertex nodes.
        nV = domain.shape.nNod
        cntVN = nm.empty( (nV,), dtype = nm.int32 )
        cntVN.fill( -1 )
        
        nodeOffsetTable[iVertex,0] = iseq
        ia = 0
        for regionName, ig, ap in self.iterAps( all = True ):
            region = ap.region
            nodeDesc = self.nodeDescs[regionName]
            nEP = ap.nEP['v']
            group = region.domain.groups[ig]
            
            ap.econn = nm.zeros( (region.shape[ig].nCell, nEP), nm.int32 )
##             print ap.econn.shape
#            pause()
            if nodeDesc.vertex.size:
                offset = group.shape.nEP
                vertices = region.getVertices( ig )
                nNew = (nm.where( cntVN[vertices] == -1 )[0]).shape[0]
                cntVN[vertices] = vertices
#                print nNew
                iseq += nNew
            nodeOffsetTable[iVertex,ia+1] = iseq
            ia += 1

        ##
        # Remap vertex node connectivity to field-local numbering.
        indx = nm.arange( iseq, dtype = nm.int32 )
        remap = nm.empty( (nV,), dtype = nm.int32 )
        remap.fill( -1 )
        remap[nm.where( cntVN >= 0 )[0]] = indx
##         print remap
##         pause()
##         print iseq, remap
##         pause()
        for regionName, ig, ap in self.iterAps( all = True ):
            region = ap.region
            nodeDesc = self.nodeDescs[regionName]
            group = region.domain.groups[ig]
            if nodeDesc.vertex.size:
                offset = group.shape.nEP
                cells = region.getCells( ig )
                ap.econn[:,:offset] = remap[group.conn[cells]]
##                 print group.conn, nm.amax( group.conn )
##                 print ap.econn, nm.amax( ap.econn )
##                 pause()

        ed, ned, fa, nfa = domain.getNeighbourLists()
        entt = self.entTable
        cntEN = nm.zeros( (entt.shape[1], ed.nUnique), nm.int32 ) - 1

        ##
        # Edge nodes.
        nodeOffsetTable[iEdge,0] = iseq
        ia = 0
        for regionName, ig, ap in self.iterAps( all = True ):
            region = ap.region
            nodeDesc = self.nodeDescs[regionName]
            group = region.domain.groups[ig]
            if nodeDesc.edge:
                cptr0 = ned.pel[ned.pg[ig]]
                ori = self.edgeOris[ig]
                iseq = mu.assignEdgeNodes( iseq, ap.econn, cntEN, \
                                           ori, entt, ned.uid, \
                                           nodeDesc.edge, cptr0 )[1]
##                 print ap.econn
##                 pause()
            nodeOffsetTable[iEdge,ia+1] = iseq
            ia += 1
            
        #    cntFN = nm.zeros( (fntt.shape[1], fa.nUnique), nm.int32 ) - 1
        nodeOffsetTable[iFace,0] = iseq
        ia = 0
        for regionName, ig, ap in self.iterAps( all = True ):
            nodeOffsetTable[iFace,ia+1] = iseq
            ia += 1

        #    cntBN = nm.zeros( (fntt.shape[1], fa.nUnique), nm.int32 ) - 1
        ##
        # Bubble nodes.
        nodeOffsetTable[iBubble,0] = iseq
        ia = 0
        for regionName, ig, ap in self.iterAps( all = True ):
            region = ap.region
            nodeDesc = self.nodeDescs[regionName]
            group = region.domain.groups[ig]
            if (nodeDesc.bubble):
                nBubble = nodeDesc.bubble[0].shape[0]
                nEl = region.shape[ig].nCell
                aux = nm.arange( iseq, iseq + nBubble * nEl )
                aux.shape = (nEl, nBubble)
                offset = nodeDesc.bubble[0][0,0]
                ap.econn[:,offset:] = aux[:,:]
                iseq += nBubble * nEl

            nodeOffsetTable[iBubble,ia+1] = iseq
            ia += 1
            
##         print nodeOffsetTable
        if nodeOffsetTable[-1,-1] != iseq:
            raise RuntimeError
        
        self.nodeOffsetTable = nodeOffsetTable

        ia = 0
        for regionName, ig, ap in self.iterAps( all = True ):
            ap.nodeOffsets = self.nodeOffsetTable[:,ia:ia+2]
            ia += 1
##             print ia
##             print ap.econn
##             print ap.nodeOffsets
#            pause()
            
        for regionName, ig, ap in self.iterAps( all = True ):
            nodeDesc = self.nodeDescs[regionName]

            gd = ap.interp.gel.data['v']
            ap.efaces = gd.faces.copy()

            if ap.hasExtraEdgeNodes:
                nd = nodeDesc.edge
                efs = []
                for eof in gd.edgesOfFaces:
                    ef = [nd[ie][:,0] for ie in eof]
                    efs.append( ef )
                efs = nm.array( efs ).squeeze()
                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
#                print efs
                ap.efaces = nm.hstack( (ap.efaces, efs ) )

            if ap.hasExtraFaceNodes:
                nd = nodeDesc.face
                efs = [ef[:,0] for ef in nd]
#                print nd
                efs = nm.array( efs ).squeeze()
                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
#                print efs
                ap.efaces = nm.hstack( (ap.efaces, efs ) )
                
##             print ap.efaces
##             print ap.interp.baseFuns['v'].fun
##             try:
##                 print ap.interp.baseFuns['s2'].fun
##                 print ap.interp.baseFuns['s3'].fun
##             except:
##                 pass
##             pause()

        return iseq, remap, cntVN, cntEN

    ##
    # 19.07.2006, c
    # 13.02.2007
    # 19.02.2007
    # 20.02.2007
    # 12.07.2007
    def setupCoors( self, mesh, cntVN ):
        """id column is set to 1."""
        noft = self.nodeOffsetTable

        nNod = noft[-1,-1]
        self.coors = nm.empty( (nNod, mesh.dim + 1), nm.float64 )
        self.coors[:,-1] = 1.0
#        print nNod

        # Mesh vertex nodes.
        inod = slice( noft[0,0], noft[0,-1] )
#        print inod
        if inod:
            indx = cntVN[cntVN >= 0]
            self.coors[inod,:] = mesh.nod0[indx,:]
##         print self.coors
##         print mesh.nod0
##         pause()
        for regionName, ig, ap in self.iterAps( all = True ):
            ap.evalExtraCoor( self.coors, mesh )

##         print self.coors
##         print self.coors.shape
##         pause()

    ##
    # created:       12.10.2005
    # last revision: 21.12.2007
    def describeGeometry( self, field, geometries, geomRequest, integral,
                          overWrite = False ):

        for regionName, ig, ap in self.iterAps( all = True,
                                                igs = geomRequest.region.igs ):
##             print regionName, ig, ap.name

            # Store integral for possible future base function request.
            ap.integrals[integral.name] = integral

            ##
            # Prepare common bases.
            gtype = geomRequest.gtype
            if gtype == 'Volume':
                ap.getBase( 'v', 0, integral = integral )
                ap.getBase( 'v', 1, integral = integral )
            elif (gtype == 'Surface') or (gtype == 'SurfaceExtra'):
                pass

            geomKey = (integral.name, geomRequest.gtype,
                       geomRequest.region.name, ap.name)
##            print field.name, geomKey

            if geometries.has_key( geomKey ):
                self.geometries[geomKey] = geometries[geomKey]
            else:
                if self.geometries.has_key( geomKey ):
                    geometries[geomKey] = self.geometries[geomKey]
                else:
##                print 'new geometry: %s of %s' % (geomKey, ap.name)
                    geom = ap.describeGeometry( field, geomRequest, integral,
                                                self.coors[:,:-1].copy() )
                    self.geometries[geomKey] = geometries[geomKey] = geom

        if overWrite:
            self.geometries = geometries

    ##
    # created:       21.12.2007
    # last revision: 21.12.2007
    def updateGeometry( self, field, regions, geometries ):
        for geomKey, geom in self.geometries.iteritems():
            iname, gtype, tregionName, apName = geomKey
            ap = self[apName]
            integral = ap.integrals[iname]
            geomRequest = Struct( gtype = gtype, region = regions[tregionName] )
            geom = ap.describeGeometry( field, geomRequest, integral,
                                        self.coors[:,:-1].copy() )
            self.geometries[geomKey] = geometries[geomKey] = geom
            
    ##
    # 03.05.2007, c
    # 17.07.2007
    def purgeSurfaceData( self ):
        for ap in self:
            ap.surfaceData = {}
