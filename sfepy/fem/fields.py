from sfepy.base.base import *
from sfepy.base.reader import Reader
import sfepy.base.la as la
import fea
from mesh import Mesh, makePointCells
import sfepy.terms as terms
import extmods.geometry as gm

##
# 14.07.2006, c
class Fields( Container ):

    ##
    # 14.07.2006, c
    def fromConf( conf ):

        objs = OneTypeList( Field )
        for key, val in conf.iteritems():
            objs.append( Field.fromConf( val ) )

        obj = Fields( objs )

        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 05.09.2007, c
    def fromFieldList( flist, qpCoors, names = None ):
        objs = OneTypeList( Field )
        objs.extend( flist )
        obj = Fields( objs )
        if names is not None:
            obj.names = names

        obj.qpCoors = qpCoors
        
        return obj
    fromFieldList = staticmethod( fromFieldList )
    
    ##
    # 14.07.2006, c
    # 18.07.2006
    def readInterpolants( self, componentDir ):

        read = Reader( componentDir )

        interps = {}

        for field in self:
            field.createInterpolants( interps, read )

        self.interps = interps

    ##
    # 18.07.2006, c
    # 19.02.2007
    def setupApproximations( self, domain ):

        for field in self:
            field.setupApproximations( domain )

    ##
    # 19.07.2006, c
    def setupGlobalBase( self ):
        gb = []
        for field in self:
            gb.append( field.setupGlobalBase() )

    ##
    # 19.07.2006, c
    def setupCoors( self ):
        for field in self:
            field.setupCoors()

##
# 14.07.2006, c
class Field( Struct ):

    ##
    # 14.07.2006, c
    # 26.04.2007
    def fromConf( conf ):

        obj = Field( name = conf.name,
                     dim = conf.dim,
                     flags = set( getattr( conf, 'flags', () ) ),
                     regionName = conf.domain,
                     bases = conf.bases )
        return obj
    fromConf = staticmethod( fromConf )

    ##
    # 18.07.2006, c
    def createInterpolants( self, interps, read ):
        """One interpolant for each base, shared if same."""
        self.interps = {}
        for regionName, baseName in self.bases.iteritems():

            if interps.has_key( baseName ):
                interp = interps[baseName]
            else:
                interp = read( fea.Interpolant, baseName )
                interps[baseName] = interp

            self.interps[baseName] = interp
            
    ##
    # 18.07.2006, c
    # 05.09.2006
    def setupApproximations( self, domain ):
        self.aps = fea.Approximations( self.bases, self.interps, domain )
        self.domain = domain
        self.region = domain.regions[self.regionName]

    ##
    #
    def igs( self ):
        return self.aps.igs

    ##
    # 19.07.2006, c
    def setupGlobalBase( self ):

        self.aps.describeNodes()
        self.aps.setupNodes()

        aux = self.aps.setupGlobalBase( self.domain )
        self.nNod, self.remap, self.cntVN, self.cntEN = aux

##         print self.nNod, self.cntVN, self.cntEN
#        pause()

    ##
    # 19.07.2006, c
    def setupCoors( self ):
        """Coordinates of field nodes."""
        self.aps.setupCoors( self.domain.mesh, self.cntVN )

    ##
    # c: 02.01.2008, r: 02.01.2008
    def getExtraNodesAsSimplices( self, iextra = None ):
        dim = self.domain.mesh.dim
        if iextra is None:
            noft = self.aps.nodeOffsetTable
            iextra = nm.arange( noft[1,0], noft[-1,-1], dtype = nm.int32 )
        extra = makePointCells( iextra, dim )
        return {2 : '2_3', 3 : '3_4'}[dim], -nm.ones_like( iextra ), extra

    ##
    # c: 19.07.2006, r: 27.02.2008
    def writeMesh( self, nameTemplate, fieldName = None ):
        """Extra nodes are written as zero-size simplices (= points)."""
        if fieldName is None:
            fieldName = self.name

        mesh = self.domain.mesh
        dim = mesh.dim

        conns, matIds, descs = [], [], []
        for regionName, ig, ap in self.aps.iterAps():
            region = ap.region
            group = region.domain.groups[ig]
            offset = group.shape.nEP
            conn = ap.econn[:,:offset]
            conns.append( conn )
            matIds.append( mesh.matIds[ig] )
            descs.append( mesh.descs[ig] )

        aux = self.getExtraNodesAsSimplices()
        descs.append( aux[0] )
        matIds.append( aux[1] )
        conns.append( aux[2] )

        tmp = Mesh.fromData( nameTemplate % fieldName,
                             self.aps.coors, conns, matIds, descs )
##         print tmp
##         pause()
        tmp.write( io = 'auto' )

    ##
    # c: 20.07.2006, r: 15.01.2008
    def getNodeDescs( self, region ):
        nds = {}
        for regionName, ig, ap in self.aps.iterAps():
            if ig in region.igs:
                nds[ig] = self.aps.nodeDescs[regionName]

        return nds

    ##
    # Modify me for bubble-only approximations to not generate vertex nodes.
    # 12.10.2005, c
    # 26.10.2005
    # 26.05.2006
    # 05.06.2006
    # 25.07.2006
    # 04.09.2006
    def interpCValsToNVals( self, vec ):
        """len( vec ) == domain.nEl"""
        nEls = [sub.nEl for sub in self.domain.subs]
        oel = nm.cumsum( [0] + nEls )
        if sum( nEls ) != vec.shape[0]:
            print 'incomatible shape! (%d == %d)' % (sum( nEls ), vec.shape[0])
            raise ValueError

        ##
        # Mesh vertex values. 
        nVertex = self.domain.nNod
        nodVol = nm.zeros( (nVertex,), nm.float64 )
        dim = vec.shape[1]
        nodVolVal = nm.zeros( (nVertex, dim ), nm.float64 )
        for ii, ap in enumerate( self.aps ):
            sub = ap.sub
            ig = sub.iseq
            vg = self.vgs[ii]
            volume = nm.squeeze( vg.variable( 2 ) );

            for ii in range( sub.conn.shape[1] ):
                cc = sub.conn[:,ii]
                nodVol[cc] += volume
                val = volume[:,nm.newaxis] * vec[oel[ig]:oel[ig+1],:]
                ind2, ind1 = nm.meshgrid( nm.arange( dim ), cc )
                nodVolVal[ind1,ind2] += val
        
        nodVolVal = nodVolVal / nodVol[:,nm.newaxis]

        ##
        # Field nodes values.
        enodVolVal = self.interpVValsToNVals( nodVolVal )

        return enodVolVal
    
    ##
    # 05.06.2006, c
    # 25.07.2006
    # 31.08.2006
    def interpVValsToNVals( self, vec ):
        dim = vec.shape[1]
        enodVolVal = nm.zeros( (self.nNod, dim), nm.float64 )
        for ii, ap in enumerate( self.aps ):
            sub = ap.sub

            noff = ap.nodeOffsets
            if noff[1] == noff[-1]:
                # Vertex values only...
                enodVolVal[sub.conn] = vec[sub.conn]
                continue

            econn = ap.econn

            ginterp = ap.interp.gel.interp
            coors = ap.interp.nodes['v'].barCoors

            qp = Struct( vals = coors )
            bf, aux = fea.evalBF( {'v': qp}, ginterp.baseFuns, ginterp.nodes )
            bf = bf['v'][:,0,:].copy()
            
            fea.mu.interpVertexData( enodVolVal, econn, vec, sub.conn, bf, 0 )

        return enodVolVal

    ##
    # 08.08.2006, c
    # 13.02.2007
    def getCoor( self, nods = None, igs = None ):
        """Will igs be ever needed?"""
        if nods is None:
            return self.aps.coors
        else:
            return self.aps.coors[nods]

    ##
    # Need to fix below, or remove!!!

    ##
    # 16.07.2007, c
    # 31.07.2007
    def getGeometry( self, kind, ig, silence = False ):

        if kind == 'volume':
            vgs, igs = self.vgs

            if not ig in igs:
                if silence:
                    return None
                else:
                    print 'volume geometry of field %s not defined'\
                          'in subdomain %d' % (field.name, ig)
                    raise ValueError

            ii = list( igs ).index( ig )
            return vgs[ii]

    ##
    # 31.08.2007, c
    def getBaseFunctions( self, kind = None, igs = None ):

        if igs is None:
            igs = self.igs

        bfs, bfgs = [], []
        for ig in igs:
            ii = self.igs.index( ig )
            ap = self.aps[ii]
            if kind is None:
                bfs.append( ap.bf )
                bfgs.append( ap.bfg )
            else:
                bfs.append( ap.bf[kind] )
                bfgs.append( ap.bfg[kind] )

        return bfs, bfgs

    ##
    # 31.08.2007, c
    def getQuadraturePoints( self, kind = None, igs = None ):

        if igs is None:
            igs = self.igs

        qps = []
        for ig in igs:
            ii = self.igs.index( ig )
            ap = self.aps[ii]
            if kind is None:
                qps.append( self.qpCoors )
            else:
                qps.append( self.qpCoors[kind] )

        return qps

    ##
    # 31.08.2007, c
    def getInterpolants( self, igs = None ):
        
        if igs is None:
            igs = self.igs

        interps = []
        for ig in igs:
            ii = self.igs.index( ig )
            ap = self.aps[ii]
            interps.append( ap.interp )

        return interps

    ##
    # 31.08.2007, c
    # 01.09.2007
    # 03.09.2007
    # 05.09.2007
    def getProlongBase( self, coors, coarseDomain, iemaps,
                        coarseInterps, suppressErrors = False ):

        ccoors = coarseDomain.getMeshCoors()
        cconns = coarseDomain.getConns()

        pbase = {}
        for ii, ap in enumerate( self.aps ):
            ig = self.igs[ii]

            bf = ap.bf['v'].squeeze()
            conn = ap.sub.conn

            iemap = iemaps[ig]
            cconn = cconns[ig]

            cinterp = coarseInterps[ii]
            cnodes = cinterp.nodes['v'].vals
            crefCoors = cinterp.nodes['v'].barCoors
            fun = cinterp.baseFuns['v'].fun

            cbfs = nm.empty( (ap.sub.nEl, bf.shape[0], 1, cconn.shape[1]),
                             dtype = nm.float64 )
            for iel, row in enumerate( conn ):
                ciel = iemap[iel]
                crow = cconn[ciel]
                
                # Reference QP -> spatial QP.
                ecoor = coors[row,:]
#                print ecoor
                xqp = nm.dot( bf, ecoor )
#                print xqp

                # Barycentric coordinates of spatial QP w.r.t. coarse element.
                ccoor = ccoors[crow,:]
#                print ccoor
                bcqp = la.barycentricCoors( xqp, ccoor )
#                print bcqp

                # Barycentric QP -> reference QP w.r.t. coarse element.
                cqp = nm.dot( bcqp.T, crefCoors )
#                print cqp

                # Coarse base function in fine QP.
                try:
                    cbf = fun.value( cqp, cnodes,
                                     suppressErrors = suppressErrors )
                except AssertionError:
                    print ig, iel, row
                    print bcqp
                    import pylab
                    fig = pylab.figure()
                    ax1 = fig.add_subplot( 121 )
                    ax1.scatter( ccoor[:,0], ccoor[:,1] )
                    ax1.scatter( xqp[:,0], xqp[:,1], c = 'r' )
                    ax1.axis( 'equal' )
                    ax2 = fig.add_subplot( 122 )
                    ax2.scatter( crefCoors[:,0], crefCoors[:,1] )
                    ax2.scatter( cqp[:,0], cqp[:,1], c = 'r' )
                    ax2.axis( 'equal' )
                    pylab.show()
                    
#                print cbf
#                print cbf.shape

                cbfs[iel,:] = cbf

            pbase[ig] = cbfs

        return pbase
