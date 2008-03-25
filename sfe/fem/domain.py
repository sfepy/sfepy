from sfe.base.base import *
from sfe.base.reader import Reader
from geomElement import GeomElement
from region import Region, sortByDependency
import fea
import sfe.base.la as la
import extmods.meshutils as mu

##
# 14.01.2005
# 22.02.2005
# 04.08.2005
def dm_createList( groups, sentinel, mode, isSort ):

    nGr = len( groups )

    if mode == 0:
        nObj = [group.shape.nEdgeTotal for group in groups.itervalues()]
        nn = 2
    else:
        nObj = [group.shape.nFaceTotal for group in groups.itervalues()]
        nn = 4
        
    objs = nm.zeros( (sum( nObj ) + 2, 3 + nn), nm.int32 ) - 1
    gptr = nm.zeros( (nGr + 1,), nm.int32 )
    eptr = nm.zeros( (nGr,), nm.int32 )

    ii = 0
    for ig, group in groups.iteritems():
        conn, gel = group.conn, group.gel
        gd = gel.data['v']
        if (mode == 0):
            nItem, items = gd.nEdge, gd.edges
        else:
            nItem, items = gd.nFace, gd.faces

#        print '::', ii
        ii = mu.createList( ii, objs, ig, conn, items, isSort )[1]
#        print '::', ii

        eptr[ig] = nItem
        gptr[ig+1] = ii

    aux = nm.repeat( nm.array( [sentinel], nm.int32 ), nn )
    objs[-2,3:] = aux
    objs[-1,3:] = aux + 1

##         print gptr, eptr

    if (ii != sum( nObj )):
        print 'neighbourList size mismatch! (%d == %d = sum( %s ))'\
              % (ii, sum( nObj ), nObj)
        raise AssertionError

    return( (objs, gptr, eptr) )

##
# c: 05.01.2005, r: 04.02.2008
def dm_neighbourList( objIn, groups, ic, perm, mode ):

    nGr = len( groups )
    nEl = [group.shape.nEl for group in groups.itervalues()]
    if mode == 0:
        nObj = [group.shape.nEdgeTotal for group in groups.itervalues()]
    else:
        nObj = [group.shape.nFaceTotal for group in groups.itervalues()]

    dataS = objIn.dataS

    pg = nm.zeros( (nGr + 1,), dtype = nm.int32 )
    pg[1:] = nEl
    pg = nm.cumsum( pg, dtype = nm.int32 )

    pel = nm.zeros( (sum( nEl ) + 1,), nm.int32 )
    for ig in range( nGr ):
        gd = groups[ig].gel.data['v']
        if (mode == 0):
            nItem = gd.nEdge
        else:
            nItem = gd.nFace
        for ie in range( nEl[ig] ):
            pel[pg[ig]+ie+1] = nItem
    pel = nm.cumsum( pel, dtype = nm.int32 )

    pobj = nm.zeros( (sum( nObj ) + 1,), dtype = nm.int32 ) + 1
    mu.neighbourListPtr( pobj, pg, pel, dataS, ic, mode )
    cntPobj = pobj.copy()
    pobj = nm.cumsum( pobj, dtype = nm.int32 )

    objs = nm.zeros( (pobj[-1],), dtype = nm.int32 )
    uid = nm.zeros( (sum( nObj ),), dtype = nm.int32 )
    cnt = nm.zeros( (sum( nObj ),), dtype = nm.int32 )
    iu = mu.neighbourList( objs, uid, cnt, \
                           pg, pel, pobj, dataS, objIn.uid, ic, perm, mode )[1]


    if (nm.sometrue( cntPobj[1:] - cnt )):
        print cntPobj[1:]
        print cnt
        print cntPobj[1:] - cnt
        raise AssertionError

    if (iu != objIn.nUnique):
        print iu, "==", objIn.nUnique
        raise AssertionError

    return( (pg, pel, pobj, objs, cnt, uid) )

##
# 05.01.2005
# 04.08.2005
# 31.10.2005
def dm_printNeighbourList( objs, objList, pauses = False ):
    pg = objList.pg
    pel = objList.pel
    pd = objList.pd
    data = objList.data
    cnt = objList.cnt
    for ig in range( len( pg ) - 1 ):
        print ig
        for iel in range( pg[ig+1] - pg[ig] ):
            print "  ", iel
            nEdge = pel[pg[ig]+iel+1] - pel[pg[ig]+iel]
            for ii in range( nEdge ):
                cptr = pel[pg[ig]+iel] + ii
                ptr = pd[cptr]
                print "    ", ii, ":", objList.uid[cptr]
                for jj in range( cnt[cptr] ):
                    print "      ", jj, data[ptr+jj], \
                          objs.data[data[ptr+jj]]

            if pauses:
                spause()

##
# 21.12.2005, c
def dm_markSurfaceFaces( fa, nfa ):
    """ flag: 0 .. inner, 1 .. triangle, 2 .. quadrangle"""
    flag = nm.zeros( (fa.data.shape[0],), nm.int32 )

    pd = nfa.pd
    data = nfa.data
    cnt = nfa.cnt
    for ii in xrange( len( cnt ) ):
        if (cnt[ii] == 1):
            ptr = pd[ii]
            if (fa.data[data[ptr],-1] == -1):
                flag[data[ptr]] = 1
            else:
                flag[data[ptr]] = 2
    return flag

##
# 17.07.2006, c
class Domain( Struct ):
    """
    Domain is divided into groups, whose purpose is to have homogeneous
    data shapes."""

    ##
    # 17.07.2006, c
    def fromMesh( mesh, componentDir ):

        read = Reader( componentDir )

        geomEls = OneTypeList( GeomElement )
        for desc in mesh.descs:
#            print desc
            
            if not geomEls.find( desc ):
                geomEls.append( read( GeomElement, desc ) )

        interps = {}
        for gel in geomEls:
            name = gel.interpolation

            if interps.has_key( name ):
                gel.interp = interps[name]
            else:
                gel.interp = read( fea.Interpolant, name )
                interps[name] = gel.interp

        obj = Domain( name = mesh.name,
                      mesh = mesh,
                      geomEls = geomEls,
                      geomInterps = interps )
        obj.matIdsToIGs = {}
        for ig, matId in enumerate( mesh.matIds ):
            obj.matIdsToIGs[matId[0]] = ig
#        obj.nNod = obj.mesh.nod0.shape[0]
        
        return obj
    fromMesh = staticmethod( fromMesh )

    def setupGroups( self ):

        for gel in self.geomEls:
            gel.setup()
        
        # Call after gel.setup()
        for ginterp in self.geomInterps.itervalues():
            gel = self.geomEls[ginterp.geometry]
            ginterp.setup( gel )

        nGr = len( self.mesh.conns )
        nNod, dim = self.mesh.nod0.shape
        self.shape = Struct( nGr = len( self.mesh.conns ), nEl = 0,
                             nNod = nNod, dim = dim )
        self.groups = {}
        for ii in range( self.shape.nGr ):
            gel = self.geomEls[self.mesh.descs[ii]] # Shortcut.
            conn = self.mesh.conns[ii]
            vertices = nm.unique1d( conn )

            nVertex = vertices.shape[0]
            nEl, nEP = conn.shape
            nEdge = gel.data['v'].nEdge           
            nEdgeTotal = nEdge * nEl

            if gel.dim == 3:
                nFace = gel.data['v'].nFace
                nFaceTotal = nFace * nEl
            else:
                nFace = nFaceTotal = 0

            shape = Struct( nVertex = nVertex, nEl = nEl, nEP = nEP,
                            nEdge = nEdge, nEdgeTotal = nEdgeTotal,
                            nFace = nFace, nFaceTotal = nFaceTotal )
            self.groups[ii] = Struct( ig = ii,
                                      vertices = vertices,
                                      conn = conn,
                                      gel = gel,
                                      shape = shape )
            self.shape.nEl += nEl

    ##
    # c: 22.11.2007, r: 25.03.2008
    def iterGroups( self, igs = None ):
        if igs is None:
            for ig in xrange( self.shape.nGr ): # sorted by ig.
                yield self.groups[ig]
        else:
            for ig in igs:
                yield ig, self.groups[ig]

    ##
    # c: 25.03.2008, r: 25.03.2008
    def getCellOffsets( self ):
        offs = {}
        off = 0
        for group in self.iterGroups():
            ig = group.ig
            offs[ig] = off
            off += group.shape.nEl
        return offs

    ##
    # 22.08.2006, c
    def getMeshCoors( self ):
        return self.mesh.nod0[:,:-1]

    ##
    # 30.08.2007, c
    def getConns( self ):
        return self.mesh.conns

    
    ##
    # 08.06.2006, c
    # 17.07.2006
    def fixElementOrientation( self ):

        coors = self.mesh.nod0
        for ii, group in self.groups.iteritems():

            ori, conn = group.gel.orientation, group.conn

            itry = 0
            while itry < 2:
                flag = -nm.ones( conn.shape[0], dtype = nm.int32 )

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                mu.orientElements( flag, conn, coors,
                                   ori.vRoots, ori.vVecs,
                                   ori.swapFrom, ori.swapTo )
    #            print flag
                if nm.alltrue( flag == 0 ):
                    if itry > 0: output( '...corrected' )
                    itry = -1
                    break

                output( 'warning: bad element orienation, trying to correct...' )
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError, "elements cannot be oriented! (%d, %s)"\
                      % (ii, self.mesh.descs[ii] )
            elif flag[0] == -1:
                output( 'warning: element orienation not checked' )

    ##
    # 19.07.2006, c
    def getOrientation( self, ii, mode = 'edges' ):
        group = self.groups[ii]
        gd = group.gel.data['v']

        if mode == 'edges':
            ori = nm.zeros( (group.shape.nEl, gd.nEdge), nm.int32 )
            mu.orientEdges( ori, group.conn, gd.edges );
        elif mode == 'faces':
            output( 'orient faces' )
            raise NotImplementedError
        else:
            output( mode )
            raise ValueError
        
        return ori

    def hasFaces( self ):
        return sum( [group.shape.nFace
                     for group in self.iterGroups()] ) > 0
        

    ##
    # c: 04.08.2005, r: 20.02.2008
    def setupNeighbourLists( self, createEdgeList = 1, createFaceList = 1 ):
        modeNames = ['edges', 'faces']
        
        isFace = self.hasFaces()
        flags = [createEdgeList, createFaceList and isFace]
        sortCols = [[3,4], [3,4,5,6]]

        for mode in range( 2 ):
            if flags[mode]:
                output( 'setting up domain %s...' % modeNames[mode] )

                tt = time.clock()
                obj = Struct()
                obj.data, obj.gptr, obj.eptr \
                          = dm_createList( self.groups,
                                           self.mesh.nod0.shape[0], mode, 1 )

#                print "t = ", time.clock() - tt
                ii = nm.arange( obj.data.shape[0], dtype = nm.int32 );
                ii.shape = (ii.shape[0], 1)
                aux = nm.concatenate( (obj.data, ii), 1 ).copy()
    ##             print aux.flags['CONTIGUOUS']
                mu.sortRows( aux, nm.array( sortCols[mode], nm.int32 ) )
    ##             print "->", aux

#                print "t = ", time.clock() - tt
                obj.perm = perm = aux[:,-1].copy()
                obj.dataS = aux[:,:-1].copy()
                aux = nm.arange( perm.shape[0], dtype = nm.int32 )
                obj.permI = nm.zeros_like( obj.perm )
                obj.permI[perm] = aux

##                 print perm
##                 print obj.permI

                ic = nm.where( nm.sum( nm.absolute( \
                    obj.dataS[1:,3:] - obj.dataS[:-1,3:] ), 1 ), 0, 1 )
		ic = nm.asarray( ic, dtype = nm.int32 )
##                 print ic, len( ic )
                obj.nData = len( ic ) - 1
                obj.nUnique = len( ic ) - nm.sum( ic ) - 1
                obj.uniqueList = nm.asarray( nm.where( ic[:-1] == 0 )[0],
					     dtype = nm.int32 )
#               print "t = ", time.clock() - tt

                assert len( obj.uniqueList ) == obj.nUnique
#                print obj.nUnique, obj.uniqueList, obj.uniqueList.shape
                ii = nm.cumsum( ic[:-1] == 0, dtype = nm.int32 )
                obj.uid = ii.copy()
                obj.uid[0], obj.uid[1:] = 0, ii[:-1]
                obj.uidI = obj.uid[obj.permI[:-2]]
##                 print obj
##                 debug()

                nobj = Struct()
                nobj.pg, nobj.pel, nobj.pd, nobj.data, nobj.cnt, nobj.uid \
                         = dm_neighbourList( obj, self.groups, ic, perm, mode )
#                print "t = ", time.clock() - tt


                if mode == 0:
#                    dm_printNeighbourList( obj, nobj, pauses = True )
                    self.ed, self.ned = obj, nobj
                else:
                    self.fa, self.nfa = obj, nobj

                output( '...done in %.2f s' % (time.clock() - tt) )
        if not isFace:
            self.fa, self.nfa = None, None

    ##
    # 19.07.2006
    # 24.08.2006
    def getNeighbourLists( self, forceFaces = False ):
        if forceFaces and not self.fa:
            return self.ed, self.ned, self.ed, self.ned
        else:
            return self.ed, self.ned, self.fa, self.nfa

    ##
    # c: 31.10.2005, r: 20.02.2008
    def createRegions( self, regionDefs, funmod = None ):
        from sfe.fem.parseReg import createBNF, visitStack, printStack,\
             ParseException

        output( 'creating regions...' )
        tt = time.clock()
        regions = OneTypeList( Region )

        ##
        # 14.06.2006, c
        # 15.06.2006
        # 19.02.2007
        # 02.03.2007
        # 02.05.2007
        # 30.05.2007
        # 05.06.2007
        def regionLeaf( domain, rdef, funmod ):
            def _regionLeaf( level, op ):

                token, details = op['token'], op['orig']
                if token != 'KW_Region':
                    parseDef = token + '<' + ' '.join( details ) + '>'
##                     conns = [group.conn for group in domain.groups.itervalues()]
##                     vertexGroups = [group.vertices
##                                     for group in domain.groups.itervalues()]
                    region = Region( 'leaf', rdef, domain, parseDef )

                if token == 'KW_Region':
                    details = details[1][2:]
                    aux = regions.find( details )
                    if not aux:
                        raise ValueError, 'region %s does not exist' % details
                    else:
                        if rdef[:4] == 'copy':
                            region = aux.copy()
                        else:
                            region = aux

                elif token == 'KW_All':
                    region.setVertices( nm.arange( domain.mesh.nod0.shape[0],
                                                   dtype = nm.int32 ) )
                elif token == 'E_NIR':
                    where = details[2]
                    
                    if where[0] == '[':
                        out = nm.array( eval( where ), dtype = nm.int32 )
                        assert nm.amin( out ) >= 0
                        assert nm.amax( out ) < domain.mesh.nod0.shape[0]
                    else:
                        x = domain.mesh.nod0[:,0]
                        y = domain.mesh.nod0[:,1]
                        if domain.mesh.dim == 3:
                            z = domain.mesh.nod0[:,2]
                        else:
                            z = None
                        coorDict = {'x' : x, 'y' : y, 'z': z}
                        
                        out = nm.where( eval( where, {}, coorDict ) )[0]
                    region.setVertices( out )
                    
                elif token == 'E_NOS':

                    if domain.fa: # 3D.
                        fa, nfa = domain.fa, domain.nfa
                    else:
                        fa, nfa = domain.ed, domain.ned
                        
                    flag = dm_markSurfaceFaces( fa, nfa )
                    ii = nm.where( flag > 0 )[0]
                    aux = la.unique1d( fa.data[ii,3:].ravel() )
                    if aux[0] == -1: # Triangular faces have -1 as 4. point.
                        aux = aux[1:]
                    region.canCells = False
                    region.setVertices( aux )

                elif token == 'E_NBF':
                    where = details[2]
                    
                    x = domain.mesh.nod0[:,0]
                    if domain.shape.dim > 1:
                        y = domain.mesh.nod0[:,1]
                        if domain.shape.dim > 2:
                            z = domain.mesh.nod0[:,2]
                        else:
                            z = None
                    else:
                        y = None
                    aux = {'x' : x, 'y' : y, 'z': z}
                        
                    fun = 'funmod.' + where
#                    print fun
                    out = nm.where( eval( fun, {'funmod' : funmod}, aux ) )[0]

                    region.setVertices( out )

                elif token == 'E_EOG':

                    group = int( details[3] )

                    ig = domain.matIdsToIGs[group]
                    group = domain.groups[ig]
                    region.setFromGroup( ig, group.vertices, group.shape.nEl )

                elif token == 'E_ONIR':
                    aux = regions[details[3][2:]]
                    region.setVertices( aux.allVertices[0:1] )

                elif token == 'E_NI':
                    region.setVertices( nm.array( [int( details[1] )],
                                                  dtype = nm.int32 ) )

                else:
                    raise NotImplementedError, token
                return region
            
            return _regionLeaf

        ##
        # 14.06.2006, c
        # 15.06.2006
        def regionOp( domain, rdef, funmod ):
            def _regionOp( level, op, item1, item2 ):

                token = op['token']
                if token == 'OA_SubN':
                    return item1.subN( item2 )
                elif token == 'OA_SubE':
                    return item1.subE( item2 )
                elif token == 'OA_AddN':
                    return item1.addN( item2 )
                elif token == 'OA_AddE':
                    return item1.addE( item2 )
                elif token == 'OA_IntersectN':
                    return item1.intersectN( item2 )
                elif token == 'OA_IntersectE':
                    return item1.intersectE( item2 )
                else:
                    raise NotImplementedError, token
            return _regionOp

        stack = []
        bnf = createBNF( stack )

        ##
        # Sort region definitions by dependencies.
        depends = re.compile( 'r\.([a-zA-Z_0-9]+)' ).search
        graph = {}
        nameToSortName = {}
        for sortName, rdef in regionDefs.iteritems():
            name, sel = rdef.name, rdef.select
#            print sortName, name, sel
            if nameToSortName.has_key( name ):
                raise 'region %s/%s already defined!' % (sortName, name)
            nameToSortName[name] = sortName

            if not graph.has_key( name ):
                graph[name] = [0]

            while len( sel ):
                aux = depends( sel )
                if aux:
                    graph[name].append( aux.group()[2:] )
                    sel = sel[aux.end():]
                else:
                    sel = ''
#        print graph

        sortedRegions = sortByDependency( graph )
#        print sortedRegions
        
        ##
        # Define regions.
        for name in sortedRegions:
            sortName = nameToSortName[name]
            rdef = regionDefs[sortName]

            stack[:] = []
            try:
                out = bnf.parseString( rdef.select )
            except ParseException:
                print 'parsing failed:', rdef
                raise

#            printStack( copy( stack ) )

            region = visitStack( stack, regionOp( self, rdef.select, funmod ),
                                 regionLeaf( self, rdef.select, funmod ) )
            if hasattr( rdef, 'forbid' ):
                fb = re.compile( '^group +\d+(\s+\d+)*$' ).match( rdef.forbid )
                if fb:
                    groups = rdef.forbid[5:].strip().split()
                    forbid = [int( ii ) for ii in groups]
                else:
                    raise SyntaxError, 'bad forbid: %s' % rdef.forbid
                forbiddenIgs = [self.matIdsToIGs[matId] for matId in forbid]
##                 print forbiddenIgs
##                 pause()
                region.deleteGroups( forbiddenIgs )
            if hasattr( rdef, 'canCells' ):
                region.switchCells( rdef.canCells )
            region.completeDescription( self.ed, self.fa )

            region.typeName = region.name
            region.name = rdef.name
            region.sortName = sortName
            
            output( ' ', region.typeName, region.name, region.sortName )
#            print region.definition
#            print region.parseDef
            regions.append( region )

        # Sort by definition name.
        regions.sort( cmp = lambda i1, i2: cmp( i1.sortName, i2.sortName ) )
        self.regions = regions
        output( '...done in %.2f s' % (time.clock() - tt) )

        return regions

    ##
    # 26.07.2007, c
    def getDiameter( self ):
        bbox = self.mesh.getBoundingBox()
        return (bbox[1,:] - bbox[0,:]).max()

    ##
    # 31.07.2007, c
    def getElementDiameters( self, ig, cells, vg, mode, square = True ):
        group = self.groups[ig]
        diameters = nm.empty( (len( cells ), 1, 1, 1), dtype = nm.float64 )
        if vg is None:
            diameters.fill( 1.0 )
        else:
            vg.getElementDiameters( diameters, group.gel.data['v'].edges,
                                    self.getMeshCoors().copy(), group.conn,
                                    cells, mode )
        if square:
            out = diameters.squeeze()
        else:
            out = nm.sqrt( diameters.squeeze() )

        return out
            
    ##
    # 29.08.2007, re-c from 00.01.18
    def surfaceFaces( self ):

        if not self.fa:
            print "no faces defined!"
            raise ValueError

        fa = Struct()
        fa.data, fa.gptr, fa.eptr \
                 = dm_createList( self.groups, self.mesh.nod0.shape[0], 1, 0 )

        flag = dm_markSurfaceFaces( fa, self.nfa )

        surfFaces = []
        itri = nm.where( flag == 1 )[0]
        if itri.size:
            surfFaces.append( fa.data[itri,3:6] )
            itet = nm.where( flag == 2 )[0]
        if itet.size:
            surfFaces.append( fa.data[itet,3:7] )

        isurf = nm.where( flag >= 1 )[0]
        if isurf.size:
            lst = fa.data[isurf,0:3]

        return lst, surfFaces
