from sfepy.base.base import *

##
# 02.12.2004, c
# 22.06.2005
order = {'m' : 0, 'f' : 1, 'e' : 2, 'i' : 3, 's' : 4}
def eld_cmpNode( objA, objB ):
    ret = cmp( order[objA[0]], order[objB[0]] )
    if (ret == 0):
        ret = cmp( objA[1:], objB[1:] )
    return( ret )

##
# 22.06.2005, c
# 19.07.2005
# 23.08.2006
def gel_translate( coors, edges, faces ):
    dim = len( coors.values()[0] )
    out = Struct()

    out.ids = coors.keys()
    out.ids.sort( eld_cmpNode )

    out.coors = nm.array( [coors[ii] for ii in out.ids], nm.float64 )

    out.trans = dict.fromkeys( out.ids )
    for key in out.trans.iterkeys():
        out.trans[key] = out.ids.index( key )
    out.trans[''] = -1

    ##
    def trEdge( edgesIn ):
        nEdge = len( edgesIn )
        edges = nm.zeros( (nEdge, 2), nm.int32 )
        for ii, edge in enumerate( edgesIn ):
            edges[ii,:] = [out.trans[edge[0]], out.trans[edge[1]]]
        return (nEdge, edges)
    ##
    def trFace( facesIn ):
        nFace = len( facesIn )
        nFP = len( facesIn[0] )

        faces = nm.zeros( (nFace, nFP), nm.int32 )
        for ii, face in enumerate( facesIn ):
            for ifa in range( nFP ):
                faces[ii,ifa] = out.trans[face[ifa]]
        return (nFace, nFP, faces)

    out.nEdge, out.edges = trEdge( edges )
    if dim == 3:
        out.nFace, out.nFP, out.faces = trFace( faces )
        # Assign edges to a face (in order).
        indx = {3: [[0, 1], [1, 2], [2, 0]],
                4: [[0, 1], [1, 2], [2, 3], [3, 0]]}
        out.edgesOfFaces = []
        se = [set( edge ) for edge in out.edges]
        iis = indx[out.nFP]
        for face in out.faces:
            eof = []
            for ii in iis:
                edge = set( face[ii] )
                ie = se.index( edge )
                eof.append( ie )
            out.edgesOfFaces.append( eof )
        out.edgesOfFaces = nm.array( out.edgesOfFaces )
##         print out.edgesOfFaces
    else:
        out.nFace, out.faces = out.nEdge, out.edges
        out.edgesOfFaces = nm.arange( out.nEdge )[:,nm.newaxis]
    out.dim = dim
    out.nVertex = len( out.ids )

    return out

##
# 08.06.2006, c
def gel_setupOrientation( vecsTuple ):
    cycle = range( len( vecsTuple ) / 4 )

    roots = nm.array( [vecsTuple[4*ii] for ii in cycle], dtype = nm.int32 )
    vecs = nm.array( [vecsTuple[4*ii+1] for ii in cycle],
                     dtype = nm.int32, ndmin = 2 )
    swapFrom = nm.array( [vecsTuple[4*ii+2] for ii in cycle],
                         dtype = nm.int32, ndmin = 2 )
    swapTo = nm.array( [vecsTuple[4*ii+3] for ii in cycle],
                       dtype = nm.int32, ndmin = 2 )

    return roots, vecs, swapFrom, swapTo


##
# Geometric element.
# 29.11.2004, c
# 30.11.2004
# 01.12.2004
# 07.12.2004
# 08.12.2004
# 09.12.2004
class GeomElement( Struct ):
    # Class methods.
    def cmpEdge( objA, objB ):
        a = objA[3:].sort()
        b = objB[3:].sort()
        
        return( cmp( a, b ) )
    cmpEdge = classmethod( cmpEdge )

    ##
    # 30.03.2005
    # 17.06.2005
    # 07.12.2005
    def __init__( self, **kwargs ):
        Struct.__init__( self, **kwargs )
        self.setupDone = 0
        
    ##
    # 17.01.2005
    # 30.03.2005
    # 22.06.2005
    # 19.07.2005
    # 08.06.2006
    # 02.08.2006
    def setup( self ):
        if (self.setupDone): return

        self.dim = len( self.vCoors.values()[0] )

        self.data = {}
        self.data['v'] = gel_translate( self.vCoors, self.vEdges, self.vFaces )
        for key in self.sEdges.keys():
            if self.dim == 2:
                self.sFaces = {key : None}
            self.data[key] = gel_translate( self.sCoors,
                                            self.sEdges[key],
                                            self.sFaces[key] )

        ori = dictToStruct( self.orientation, flag = (1,) )
        ori.vRoots, ori.vVecs, ori.swapFrom, ori.swapTo =\
                    gel_setupOrientation( ori.vVecs )
        self.orientation = ori
        self.setupDone = 1
