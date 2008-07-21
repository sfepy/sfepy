from sfepy.base.base import *

##
# 02.12.2004, c
# 22.06.2005
order = {'m' : 0, 'f' : 1, 'e' : 2, 'i' : 3, 's' : 4}
def eld_cmp_node( obj_a, obj_b ):
    ret = cmp( order[obj_a[0]], order[obj_b[0]] )
    if (ret == 0):
        ret = cmp( obj_a[1:], obj_b[1:] )
    return( ret )

##
# 22.06.2005, c
# 19.07.2005
# 23.08.2006
def gel_translate( coors, edges, faces ):
    dim = len( coors.values()[0] )
    out = Struct()

    out.ids = coors.keys()
    out.ids.sort( eld_cmp_node )

    out.coors = nm.array( [coors[ii] for ii in out.ids], nm.float64 )

    out.trans = dict.fromkeys( out.ids )
    for key in out.trans.iterkeys():
        out.trans[key] = out.ids.index( key )
    out.trans[''] = -1

    ##
    def tr_edge( edges_in ):
        n_edge = len( edges_in )
        edges = nm.zeros( (n_edge, 2), nm.int32 )
        for ii, edge in enumerate( edges_in ):
            edges[ii,:] = [out.trans[edge[0]], out.trans[edge[1]]]
        return (n_edge, edges)
    ##
    def tr_face( faces_in ):
        n_face = len( faces_in )
        n_fp = len( faces_in[0] )

        faces = nm.zeros( (n_face, n_fp), nm.int32 )
        for ii, face in enumerate( faces_in ):
            for ifa in range( n_fp ):
                faces[ii,ifa] = out.trans[face[ifa]]
        return (n_face, n_fp, faces)

    out.n_edge, out.edges = tr_edge( edges )
    if dim == 3:
        out.n_face, out.n_fp, out.faces = tr_face( faces )
        # Assign edges to a face (in order).
        indx = {3: [[0, 1], [1, 2], [2, 0]],
                4: [[0, 1], [1, 2], [2, 3], [3, 0]]}
        out.edges_of_faces = []
        se = [set( edge ) for edge in out.edges]
        iis = indx[out.n_fp]
        for face in out.faces:
            eof = []
            for ii in iis:
                edge = set( face[ii] )
                ie = se.index( edge )
                eof.append( ie )
            out.edges_of_faces.append( eof )
        out.edges_of_faces = nm.array( out.edges_of_faces )
##         print out.edges_of_faces
    else:
        out.n_face, out.faces = out.n_edge, out.edges
        out.edges_of_faces = nm.arange( out.n_edge )[:,nm.newaxis]
    out.dim = dim
    out.n_vertex = len( out.ids )

    return out

##
# 08.06.2006, c
def gel_setup_orientation( vecs_tuple ):
    cycle = range( len( vecs_tuple ) / 4 )

    roots = nm.array( [vecs_tuple[4*ii] for ii in cycle], dtype = nm.int32 )
    vecs = nm.array( [vecs_tuple[4*ii+1] for ii in cycle],
                     dtype = nm.int32, ndmin = 2 )
    swap_from = nm.array( [vecs_tuple[4*ii+2] for ii in cycle],
                         dtype = nm.int32, ndmin = 2 )
    swap_to = nm.array( [vecs_tuple[4*ii+3] for ii in cycle],
                       dtype = nm.int32, ndmin = 2 )

    return roots, vecs, swap_from, swap_to


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
    def cmp_edge( obj_a, obj_b ):
        a = obj_a[3:].sort()
        b = obj_b[3:].sort()
        
        return( cmp( a, b ) )
    cmp_edge = classmethod( cmp_edge )

    ##
    # 30.03.2005
    # 17.06.2005
    # 07.12.2005
    def __init__( self, **kwargs ):
        Struct.__init__( self, **kwargs )
        self.setup_done = 0
        
    ##
    # 17.01.2005
    # 30.03.2005
    # 22.06.2005
    # 19.07.2005
    # 08.06.2006
    # 02.08.2006
    def setup( self ):
        if (self.setup_done): return

        self.dim = len( self.v_coors.values()[0] )

        self.data = {}
        self.data['v'] = gel_translate( self.v_coors, self.v_edges, self.v_faces )
        for key in self.s_edges.keys():
            if self.dim == 2:
                self.s_faces = {key : None}
            self.data[key] = gel_translate( self.s_coors,
                                            self.s_edges[key],
                                            self.s_faces[key] )

        ori = dict_to_struct( self.orientation, flag = (1,) )
        ori.v_roots, ori.v_vecs, ori.swap_from, ori.swap_to =\
                    gel_setup_orientation( ori.v_vecs )
        self.orientation = ori
        self.setup_done = 1
