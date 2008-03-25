from sfe.base.base import *
import sfe.base.la as la

##
# 15.06.2006, c
# 17.07.2006
# 04.09.2006
def sortByDependency( graph ):

    out = []

    nNod = len( graph )
    idone = 0
    idone0 = -1
    while idone < nNod:

        depRemoved = 0
        for node, deps in graph.iteritems():
#            print '--', node, deps
            if (len( deps ) == 1) and not deps[0]:
                out.append( node )
                deps[0] = 1
                idone += 1
            elif not deps[0]:
#                print '--->', deps
                for ii, dep in enumerate( deps[1:] ):
                    if graph[dep][0]:
                        ir = deps.index( dep )
                        deps.pop( ir )
                        depRemoved += 1
#                print '---<', deps

##         print graph
##         print out
##         print idone, idone0, nNod, depRemoved
##         pause()

        if (idone <= idone0) and not depRemoved:
            raise ValueError, 'circular dependency'
        idone0 = idone

    return out

##
# 15.06.2006, c
def _join( def1, op, def2 ):
    return '(' + def1 + ' ' + op + ' ' + def2 + ')'

##
# 31.10.2005, c
class Region( Struct ):

    ##
    # 14.06.2006, c
    # 15.06.2006
    # 23.02.2007
    def __init__( self, name, definition, domain, parseDef ):
        """conns, vertexGroups are links to domain data"""
        Struct.__init__( self,
                         name = name, definition = definition,
                         nVMax = domain.shape.nNod, domain = domain,
                         parseDef = parseDef, allVertices = None,
                         igs = [], vertices = {}, edges = {}, faces = {},
                         cells = {}, fis = {},
                         volume = {}, surface = {}, length = {},
                         canCells = True, mustUpdate = True,
                         isComplete = False )

    ##
    # 15.06.2006, c
    def lightCopy( self, name, parseDef ):
        return Region( name, self.definition, self.domain, parseDef )

    ##
    # c: 15.06.2006, r: 04.02.2008
    def updateGroups( self, force = False ):
        """Vertices common to several groups are listed only in all of them -
        fa, ed.uniqueIndx contain no edge/face duplicates already."""
        if self.mustUpdate or force:

            self.igs = []
            self.vertices = {}
            self.cells = {}

            for group in self.domain.iterGroups():
                ig = group.ig
                vv = nm.intersect1d( group.vertices, self.allVertices )
                if len( vv ) == 0: continue

                self.igs.append( ig )
                self.vertices[ig] = vv

                if self.canCells:
                    mask = nm.zeros( self.nVMax, nm.int32 )
                    mask[vv] = 1

                    conn = group.conn
                    aux = nm.sum( mask[conn], 1, dtype = nm.int32 )
                    rcells = nm.where( aux == conn.shape[1] )[0]
                    self.cells[ig] = nm.asarray( rcells, dtype = nm.int32 )
        self.mustUpdate = False

    ##
    # 15.06.2006, c
    def updateVertices( self ):
        self.allVertices = nm.zeros( (0,), nm.int32 )
        self.vertices = {}
        for ig, group in self.domain.iterGroups( self.igs ):
            rcells = self.cells[ig]
            conn = group.conn
            nods = conn[rcells,:].ravel()
            aux = la.unique1d( nods )
            self.vertices[ig] = aux
            self.allVertices = nm.unique1d( nm.r_[self.allVertices, aux] )
        
    ##
    # 15.06.2006, c
    def setVertices( self, vertices ):

        self.allVertices = vertices
        self.updateGroups( force = True )
        self.isComplete = False

    ##
    # 15.06.2006, c
    def setCells( self, cells ):

        self.cells = cells
        self.updateVertices()
        self.isComplete = False

    ##
    # 15.06.2006, c
    def setFromGroup( self, ig, vertices, nCell ):

        self.igs = [ig]
        self.cells = {ig : nm.arange( nCell, dtype = nm.int32 )}
        self.vertices = {ig: vertices.copy()}
        self.allVertices = vertices.copy()
        self.mustUpdate = False

    ##
    # c: 23.02.2007, r: 22.01.2008
    def deleteGroups( self, digs ):
        """self.completeDescription must be called after!"""
        for ig in digs:
            try:
                del self.vertices[ig]
                del self.cells[ig]
                self.igs.remove( ig )
            except KeyError:
                pass

    ##
    # 17.07.2007, c
    def switchCells( self, canCells ):
        if self.canCells:
            self.canCells = canCells
            if not canCells:
                self.cells = {}
        else:
            self.canCells = canCells
            if canCells:
                self.updateGroups( force = True )
        
    ##
    # 15.06.2006, c
    # 02.08.2006
    # 21.02.2007
    # 23.02.2007
    # 30.03.2007
    def completeDescription( self, ed, fa ):
        """self.edges, self.faces list edge/face indices per group (pointers to
        ed.data, fa.data) - repetitions among groups are possible."""
        ##
        # Get edges, faces, etc. par subdomain.
        edges = ed.data

        if fa is not None:
            faces = fa.data.copy()
            # Treat both 3, 4 node faces.
            ii = nm.where( faces[:,-1] == -1 )[0]
            faces[ii,-1] = faces[ii,3]

        self.edges = {}
        self.faces = {}
        self.shape = {}
        for ig, group in self.domain.iterGroups( self.igs ):
            vv = self.vertices[ig]
            if self.cells.has_key( ig ):
                nCell = self.cells[ig].shape[0]
            else:
                nCell = 0
            self.shape[ig] = Struct( nVertex = vv.shape[0],
                                     nCell = nCell )
            if len( vv ) == 0: continue

            mask = nm.zeros( self.nVMax, nm.int32 )
            mask[vv] = 1

            ied = nm.arange( ed.gptr[ig], ed.gptr[ig+1], dtype = nm.int32 )
            aux = nm.sum( mask[edges[ied,3:5]], 1 )
            # Points to ed.data.
            redges = ied[nm.where( aux == 2 )[0]]
            self.edges[ig] = redges
            if fa is None: continue
            
            ifa = nm.arange( fa.gptr[ig], fa.gptr[ig+1], dtype = nm.int32 )
            aux = nm.sum( mask[faces[ifa,3:7]], 1 )
            # Points to fa.data.
            rfaces = ifa[nm.where( aux == 4 )[0]]
            self.faces[ig] = rfaces

            self.shape[ig].nEdge = redges.shape[0]
            self.shape[ig].nFace = rfaces.shape[0]
            
        self.isComplete = True

    ##
    # 24.08.2006, c
    # 16.02.2007
    # 21.02.2007
    # 23.02.2007
    def setupFaceIndices( self, fa ):
        """(iel, ifa) for each face."""
        if self.faces:
            faces = self.faces
        else:
            faces = self.edges

        self.fis = {}
        for ig in self.igs:
            rfaces = faces[ig]
            aux = fa.data[rfaces]
            assert nm.all( aux[:,0] == ig )
            fi = aux[:,1:3].copy()
            self.fis[ig] = fi

    ##
    # 05.09.2006, c
    # 22.02.2007
    # 17.07.2007
    def selectCells( self, nVerts ):
        """Select cells containing at least nVerts[ii] vertices per group ii."""
        if not self.canCells:
            print 'region %s cannot have cells!' % self.name
            raise ValueError
        
        self.cells = {}
        for ig, group in self.domain.iterGroups( self.igs ):
            vv = self.vertices[ig]
            if len( vv ) == 0: continue
            
            mask = nm.zeros( self.nVMax, nm.int32 )
            mask[vv] = 1

            aux = nm.sum( mask[group.conn], 1 )
            rcells = nm.where( aux >= nVerts[ig] )[0]
#            print rcells.shape
            self.cells[ig] = rcells

    ##
    # 22.02.2007, c
    # 17.07.2007
    def selectCellsOfSurface( self ):
        """Select cells corresponding to faces (or edges in 2D)."""
        if not self.canCells:
            print 'region %s cannot have cells!' % self.name
            raise ValueError

        if self.faces:
            faces = self.faces
        else:
            faces = self.edges

        self.cells = {}
        for ig in self.igs:
            rcells = self.fis[ig][:,0]
            self.cells[ig]= rcells

    ##
    # 02.03.2007, c
    def copy( self ):
        """Vertices-based copy."""
        tmp = self.lightCopy( 'copy', self.parseDef )
        tmp.setVertices( copy( self.allVertices ) )
        return tmp
        
    ##
    # 15.06.2006, c
    def subN( self, other ):
        tmp = self.lightCopy( 'op',
                              _join( self.parseDef, '-n', other.parseDef ) )
        tmp.setVertices( nm.setdiff1d( self.allVertices,
                                       other.allVertices ) )
        
        return tmp

    ##
    # 15.06.2006, c
    def addN( self, other ):
        tmp = self.lightCopy( 'op',
                              _join( self.parseDef, '+n', other.parseDef ) )
        tmp.setVertices( nm.union1d( self.allVertices,
                                     other.allVertices ) )
        
        return tmp

    ##
    # 15.06.2006, c
    def intersectN( self, other ):
        tmp = self.lightCopy( 'op',
                              _join( self.parseDef, '*n', other.parseDef ) )
        tmp.setVertices( nm.intersect1d( self.allVertices,
                                         other.allVertices ) )
        
        return tmp

    ##
    # c: 15.06.2006, r: 17.03.2008
    def subE( self, other ):
        tmp = self.lightCopy( 'op',
                              _join( self.parseDef, '-e', other.parseDef ) )
        for ig in self.igs:
            if ig not in other.igs:
                tmp.igs.append( ig )
                tmp.cells[ig] = self.cells[ig].copy()
                continue
            
            aux = nm.setdiff1d( self.cells[ig], other.cells[ig] )
            if not len( aux ): continue
            tmp.cells[ig] = aux

        tmp.updateVertices()
        return tmp

    ##
    # 15.06.2006, c
    def addE( self, other ):
        tmp = self.lightCopy( 'op',
                              _join( self.parseDef, '+e', other.parseDef ) )
        for ig in self.igs:
            tmp.igs.append( ig )
            if ig not in other.igs:
                tmp.cells[ig] = self.cells[ig].copy()
                continue

            tmp.cells[ig] = nm.union1d( self.cells[ig],
                                        other.cells[ig] )

        for ig in other.igs:
            if ig in tmp.igs: continue
            tmp.igs.append( ig )
            tmp.cells[ig] = other.cells[ig].copy()

        tmp.updateVertices()
        return tmp

    ##
    # 15.06.2006, c
    # 20.02.2007
    def intersectE( self, other ):
        tmp = self.lightCopy( 'op',
                              _join( self.parseDef, '*e', other.parseDef ) )
        for ig in self.igs:
            if ig not in other.igs: continue
            aux = nm.intersect1d( self.cells[ig], other.cells[ig] )
            if not len( aux ): continue
            tmp.igs.append( ig )
            tmp.cells[ig] = aux

        tmp.updateVertices()
        return tmp

    ##
    # 16.10.2006, c
    # 20.02.2007
    # 28.02.2007
    # 31.05.2007
    # 10.07.2007
    # 16.07.2007
    # 17.07.2007
    def getFieldNodes( self, field, merge = False, igs = None ):
        """For one edge node type only! (should index row of cntEN...)"""
        if igs is None:
            igs = self.igs
        cntEN = field.cntEN

        nods = []
        nodeDescs = field.getNodeDescs( self )
        for ig, nodeDesc in nodeDescs.iteritems():
            if not ig in igs:
                nods.append( None )
                continue
            
            nnew = nm.empty( (0,), dtype = nm.int32 )
            if nodeDesc.vertex.size:
                nnew = nm.concatenate( (nnew, field.remap[self.vertices[ig]]) )

            if nodeDesc.edge:
                ed = field.domain.ed
                # ed.uidI[self.edges[ii]] == ed.uid[ed.permI[self.edges[ii]]]
                enods = cntEN[:cntEN.shape[0],ed.uidI[self.edges[ig]]].ravel()
                enods = nm.compress( (enods >= 0), enods )
                nnew = nm.concatenate( (nnew, enods) )

            if nodeDesc.face:
                print self.name, field.name
                raise NotImplementedError

            if nodeDesc.bubble and self.canCells:
                noft = field.aps.nodeOffsetTable
                ia = field.aps.igs.index( ig )
                enods = self.cells[ig] + noft[3,ia]
                nnew = nm.concatenate( (nnew, enods) )

            nods.append( nnew )

        if merge:
            nods = [nn for nn in nods if nn is not None]
            nods = nm.unique1d( nm.hstack( nods ) )
            
        return nods

    ##
    # 22.02.2007, c
    def getVertices( self, ig ):
        return self.vertices[ig]

    ##
    # 05.06.2007, c
    def getEdges( self, ig ):
        return self.edges[ig]
        
    ##
    # 05.06.2007, c
    def getFaces( self, ig ):
        return self.faces[ig]
        
    ##
    # 05.06.2007, c
    def getCells( self, ig ):
        return self.cells[ig]
        
    ##
    # created:       28.05.2007
    # last revision: 11.12.2007
    def hasCells( self ):

        if self.canCells:
            for cells in self.cells.itervalues():
                if cells.size:
                    return True
            return False
        else:
            return False

    ##
    # c: 16.07.2007, r: 17.03.2008
    def updateGeometryInfo( self, field, key ):
        """
        key: iname, aregionName, ig
        TODO: surfaces, lengths
        ?call for all regions & fields in describeGeometry()?"""
        if self.hasCells():
            val = self.volume.setdefault( field.name, {} )

        iname, ig = key
        aps = field.aps
        geometries = aps.geometries
        ap = aps.apsPerGroup[ig]
        gKey = (iname, 'Volume', self.name, ap.name)

        vg = geometries[gKey]

        volume = vg.variable( 2 )
        volume.shape = (nm.prod( volume.shape ),)

        val[key] = nm.sum( volume[self.cells[ig]] )
        self.volume[field.name] = val

    ##
    # created:       16.07.2007
    # last revision: 13.12.2007
    def getVolume( self, field, key = None,
                   update = False ):
        if update:
            self.updateGeometryInfo( field, key )

        if key is None:
            return self.volume[field.name]
        else:
            return self.volume[field.name][key]

    def contains( self, other ):
        """Tests only igs for now!!!"""
        return set( other.igs ).issubset( set( self.igs ) )

    ##
    # c: 25.03.2008, r: 25.03.2008
    def getCellOffsets( self ):
        offs = {}
        off = 0
        for ig in self.igs:
            offs[ig] = off
            off += self.shape[ig].nCell
        return offs
