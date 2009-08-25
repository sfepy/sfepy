##
# 24.11.2004, c
from sfepy.base.base import *
import sfepy.base.la as la
from geomElement import GeomElement

#from core import *

##
# 19.07.2005, c
# 20.07.2005
# 02.08.2005
# 21.09.2005
# 22.09.2005
# 27.11.2005
# 16.12.2005
def ap_nodegen_lagrange( ao, key, nod, gel, want_bar_coors = 0 ):
    ##
    def gen_edges( nodes, nts, iseq, nt, edges, ao ):
        for ii, edge in enumerate( edges ):
            n1 = nodes[edge[0],:].copy()
            n2 = nodes[edge[1],:].copy()
            for ie in range( ao - 1 ):
                c2 = ie + 1
                c1 = ao - c2
                nts[iseq] = [nt, ii]
                aux = [int( tmp ) for tmp in delta * (c1 * n1 + c2 * n2)]
                nodes[iseq,:] = aux
                iseq += 1
        return iseq

    ##
    def gen_faces( nodes, nts, iseq, nt, faces, ao ):
        for ii, face in enumerate( faces ):
            n1 = nodes[face[0],:].copy()
            n2 = nodes[face[1],:].copy()
            n3 = nodes[face[2],:].copy()
            for i1 in range( ao - 2 ):
                for i2 in range( ao - 2 - i1 ):
                    c3 = i1 + 1
                    c2 = i2 + 1
                    c1 = ao - c3 - c2
                    nts[iseq] = [nt, ii]
                    aux = [int( tmp ) for tmp
                           in delta * (c1 * n1 + c2 * n2 + c3 * n3)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

    ##
    def gen_bubble( nodes, nts, iseq, nt, ao ):
        n1 = nodes[0,:].copy()
        n2 = nodes[1,:].copy()
        n3 = nodes[2,:].copy()
        n4 = nodes[3,:].copy()
        for i1 in range( ao - 3 ):
            for i2 in range( ao - 3 ):
                for i3 in range( ao - 3 - i1 - i2 ):
                    c4 = i1 + 1
                    c3 = i2 + 1
                    c2 = i3 + 1
                    c1 = ao - c4 - c3 - c2
                    nts[iseq] = [nt, 0]
                    aux = [int( tmp ) for tmp
                           in delta * (c1 * n1 + c2 * n2 + c3 * n3 + c4 * n4)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq
                    
    ##
    fac = lambda n : reduce( lambda a, b : a * (b + 1), xrange( n ), 1 )

    gd = gel.data[key]
    dim = gd.dim
    n_v = gd.n_vertex

    n_nod = fac( ao + dim ) / (fac( ao ) * fac( dim ))
##     print n_nod, gd
    nodes = nm.zeros( (n_nod, n_v), nm.int32 )
    nts = nm.zeros( (n_nod, 2), nm.int32 )

    if ao == 0:
        nts[0,:] = [3, 0]
        nodes[0,:] = nm.zeros( (n_v,), nm.int32 )
    else:
        iseq = 0
        delta = 1.0 / float( ao )
        # Vertex nodes.
        nts[0:n_v,0] = 0
        nts[0:n_v,1] = nm.arange( n_v, dtype = nm.int32 )
        aux = ao * nm.identity( n_v, dtype = nm.int32 )
        nodes[iseq:iseq+n_v,:] = aux
        iseq += n_v

        if dim == 1:
            iseq = gen_edges( nodes, nts, iseq, 3, [[0, 1]], ao )
        elif dim == 2:
            iseq = gen_edges( nodes, nts, iseq, 1, gd.edges, ao )
            iseq = gen_faces( nodes, nts, iseq, 3, [[0, 1, 2]], ao )
        elif dim == 3:
            iseq = gen_edges( nodes, nts, iseq, 1, gd.edges, ao )
            iseq = gen_faces( nodes, nts, iseq, 2, gd.faces, ao )
            iseq = gen_bubble( nodes, nts, iseq, 3, ao )
        else:
            raise NotImplementedError
        
##     print nm.concatenate( (nts, nodes), 1 )

    # Check orders.
    naos = nm.sum( nodes, 1 )
    if not nm.alltrue( naos == ao ):
        print nodes
        raise AssertionError

    if want_bar_coors:
        ##
        # Barycentric coordinates.
        if ao == 0:
            tmp = nm.ones( (n_nod, n_v), nm.int32 )
            bar_coors = nm.dot( tmp, gd.coors ) / n_v
        else:
            bar_coors = nm.dot( nodes, gd.coors ) / ao

        return (nts, nodes, bar_coors)
    else:
        return (nts, nodes)

##
# 23.11.2005, c
# 16.12.2005
def ap_nodegen_lagrange_pxb( ao, key, nod, gel, want_bar_coors = 0 ):
    
    nts, nodes = ap_nodegen_lagrange( ao, key, nod, gel, 0 )
    shape = [nts.shape[0] + 1, 2]
    nts = nm.resize( nts, shape )
    nts[-1,:] = [3, 0]

    shape = [nodes.shape[0] + 1, nodes.shape[1]]
    nodes = nm.resize( nodes, shape )
    nodes[-1,0] = -1
    nodes[-1,1:] = 0
    
    if want_bar_coors:
        gd = gel.data[key]
        bar_coors = nm.dot( nodes, gd.coors ) / ao
        n_v = gel.data[key].n_vertex
        tmp = nm.ones( (n_v,), nm.int32 )
        bar_coors[-1,:] =  nm.dot( tmp, gd.coors ) / n_v
        return (nts, nodes, bar_coors)
    else:
        return (nts, nodes)

##
# 06.12.2005, c
# 07.12.2005
# 02.08.2006
def ap_nodegen_lagrange_tensor111( ao, key, nod, gel, want_bar_coors = 0 ):

    ##
    def gen_edges( nodes, nts, iseq, nt, edges, ao ):
        delta = 1.0 / float( ao )
        for ii, edge in enumerate( edges ):
            n1 = nodes[edge[0],:].copy()
            n2 = nodes[edge[1],:].copy()
            for ie in range( ao - 1 ):
                c2 = ie + 1
                c1 = ao - c2
                nts[iseq] = [nt, ii]
                aux = [int( tmp ) for tmp in delta * (c1 * n1 + c2 * n2)]
                nodes[iseq,:] = aux
                iseq += 1
        return iseq

    ##
    def gen_faces( nodes, nts, iseq, nt, faces, ao ):
        delta = 1.0 / (float( ao ) ** 2)
        for ii, face in enumerate( faces ):
            n1 = nodes[face[0],:].copy()
            n2 = nodes[face[1],:].copy()
            n3 = nodes[face[2],:].copy()
            n4 = nodes[face[3],:].copy()
            for i1 in range( ao - 1 ):
                for i2 in range( ao - 1 ):
                    c4 = i1 + 1
                    c3 = i2 + 1
                    c2 = ao - c4
                    c1 = ao - c3
                    nts[iseq] = [2, ii]
                    aux = [int( tmp ) for tmp
                           in delta * (c1 * c2 * n1 + c2 * c3 * n2
                                       + c3 * c4 * n3 + c4 * c1 * n4)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq
    ##
    def gen_bubble( nodes, nts, iseq, nt, ao ):
        delta = 1.0 / (float( ao ) ** 3)
        n1 = nodes[0,:].copy()
        n2 = nodes[1,:].copy()
        n3 = nodes[2,:].copy()
        n4 = nodes[3,:].copy()
        n5 = nodes[4,:].copy()
        n6 = nodes[5,:].copy()
        n7 = nodes[6,:].copy()
        n8 = nodes[7,:].copy()
        for i1 in range( ao - 1 ):
            for i2 in range( ao - 1 ):
                for i3 in range( ao - 1 ):
                    c6 = i1 + 1
                    c5 = i2 + 1
                    c4 = i3 + 1
                    c3 = ao - c6
                    c2 = ao - c5
                    c1 = ao - c4
                    nts[iseq] = [nt, 0]
                    aux = [int( tmp ) for tmp
                           in delta * (c1 * c2 * c3 * n1 + c4 * c2 * c3 * n2
                                       + c5 * c4 * c3 * n3 + c1 * c3 * c5 * n4
                                       + c1 * c2 * c6 * n5 + c4 * c2 * c6 * n6
                                       + c5 * c4 * c6 * n7 + c1 * c6 * c5 * n8)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

    # Requires fixed vertex numbering!
    vertex_maps = {3 : [[0, 0, 0],
                      [ao, 0, 0],
                      [ao, ao, 0],
                      [0, ao, 0],
                      [0, 0, ao],
                      [ao, 0, ao],
                      [ao, ao, ao],
                      [0, ao, ao]],
                  2 : [[0, 0],
                      [ao, 0],
                      [ao, ao],
                      [0, ao]],
                  1 : [0,
                       ao] }
    gd = gel.data[key]
    dim = gd.dim
    vertex_map = vertex_maps[dim]

    ##
    # Make tensor product of 1D nodes.
    n_nod = (ao + 1) ** dim
    n_v = gd.n_vertex
    nodes = nm.zeros( (n_nod, 2 * dim), nm.int32 )
    nts = nm.zeros( (n_nod, 2), nm.int32 )

    if ao == 0:
        nts[0,:] = [3, 0]
        nodes[0,:] = nm.zeros( (n_nod,), nm.int32 )
    else:
        iseq = 0
        ##
        # Vertex nodes.
        nts[0:n_v,0] = 0
        nts[0:n_v,1] = nm.arange( n_v, dtype = nm.int32 )
        aux = ao * nm.identity( n_v, dtype = nm.int32 )
        if dim == 3:
            for ii in range( n_v ):
                i1, i2, i3 = vertex_map[ii]
                nodes[iseq,:] = [ao - i1, i1, ao - i2, i2, ao - i3, i3]
                iseq += 1
        elif dim == 2:
            for ii in range( n_v ):
                i1, i2 = vertex_map[ii]
                nodes[iseq,:] = [ao - i1, i1, ao - i2, i2]
                iseq += 1
        else:
            for ii in range( n_v ):
                i1 = vertex_map[ii]
                nodes[iseq,:] = [ao - i1, i1]
                iseq += 1

        if dim == 1:
            iseq = gen_edges( nodes, nts, iseq, 3, [[0, 1]], ao )
        elif dim == 2:
            iseq = gen_edges( nodes, nts, iseq, 1, gd.edges, ao )
            iseq = gen_faces( nodes, nts, iseq, 3, [[0, 1, 2, 3]], ao )
        elif dim == 3:
            iseq = gen_edges( nodes, nts, iseq, 1, gd.edges, ao )
            iseq = gen_faces( nodes, nts, iseq, 2, gd.faces, ao )
            iseq = gen_bubble( nodes, nts, iseq, 3, ao )
        else:
            raise NotImplementedError
        
##     print nm.concatenate( (nts, nodes), 1 )

##     print nodes
##     print nts
##     pause()

    # Check orders.
    naos = nm.sum( nodes, 1 )
    if not nm.alltrue( naos == (ao * dim) ):
        print nm.concatenate( (nodes, naos[:,nm.newaxis]), 1 )
        raise AssertionError

    if want_bar_coors:
        ##
        # Barycentric coordinates.
        if ao == 0:
            tmp = nm.ones( (n_nod, n_v), nm.int32 )
            bar_coors = nm.dot( tmp, gd.coors ) / n_v

        else:
            c_min = nm.amin( gd.coors, 0 )
            c_max = nm.amax( gd.coors, 0 )
    ##         print gd.coors
    ##         print c_min, c_max

            cr = nm.arange( 2 * dim )
            bar_coors = (nodes[:,cr[::2]] * c_min
                         + nodes[:,cr[1::2]] * c_max) / ao
##             print bar_coors
        
        return (nts, nodes, bar_coors)
    else:
        return (nts, nodes)

##
# 20.07.2005, c
# 21.07.2005
# 12.10.2005
# 05.06.2006
def ap_bfgen_lagrange( *args ):

    class Generator( Struct ):
        def __init__( self, v_coors, bf, nodes, var_set ):
            if bf.grad > 1:
                raise NotImplementedError

            self.nodes = nodes
            self.n_v = self.nodes.shape[1]
            self.ao = nm.sum( nodes[0,:] )
            mtx = nm.ones( (self.n_v, self.n_v), nm.float64 )
            mtx[0:self.n_v-1,:] = nm.transpose( v_coors )
            self.mtx_i = nla.inv( mtx )
            self.rhs = nm.ones( (self.n_v,), nm.float64 )
            
        # 01.09.2007
        def __call__( self, coors, node, iv = None,
                      suppress_errors = False, eps = 1e-15 ):

            # Barycentric coordinates.
            bc = nm.zeros( (self.n_v, coors.shape[0]), nm.float64 )
            for ic, coor in enumerate( coors ):
                self.rhs[0:self.n_v-1] = coor
                bc[:,ic] = nm.dot( self.mtx_i, self.rhs )
                error = False
                for ii, val in enumerate( bc[:,ic] ):
                    if val < 0.:
                        if val > (-eps):
                            bc[ii,ic] = 0.
                        else:
                            error = True
                    if val > 1.:
                        if val < (1 + eps):
                            bc[ii,ic] = 1.
                        else:
                            error = True
                    if error and not suppress_errors:
                        msg = 'quadrature point outside of element!'
                        msg += '\nic: %s coor: %s node: %s iv: %s'\
                               % (ic, coor, node, iv)
                        msg += '\nbc: %s' % bc
                        raise AssertionError(msg)

            if iv is None:
                val = nm.ones( (coors.shape[0],), nm.float64 )
                for i1 in range( self.n_v ):
                    for i2 in range( node[i1] ):
                        val *= (self.ao * bc[i1,:] - i2) / (i2 + 1.0)
            else:
                val = nm.zeros( (coors.shape[0],), nm.float64 )
                for ii in range( self.n_v ):
                    vv = nm.ones( (coors.shape[0],), nm.float64 )
                    for i1 in range( self.n_v ):
                        if i1 == ii: continue
                        for i2 in range( node[i1] ):
                            vv *= (self.ao * bc[i1,:] - i2) / (i2 + 1.0)

                    dval = nm.zeros( (coors.shape[0],), nm.float64 )
                    for i1 in range( node[ii] ):
                        dd = nm.ones( (coors.shape[0],), nm.float64 )
                        for i2 in range( node[ii] ):
                            if i1 == i2: continue
                            dd *= (self.ao * bc[ii,:] - i2) / (i2 + 1.0)
                        dval += dd * self.ao / (i1 + 1.0)

                    val += vv * dval * self.mtx_i[ii,iv]

            return val

    return Generator( *args )
    
##
# 23.11.2005, c
# 19.12.2005
def ap_bfgen_lagrange_pxb( *args ):

    class Generator( Struct ):
        def __init__( self, v_coors, bf, nodes, var_set ):
            if bf.grad > 1:
                raise NotImplementedError

            self.bf_gen = ap_bfgen_lagrange( v_coors, bf, nodes, var_set )
            self.n_f = nodes.shape[0] - 1

            # Make a 'hypercubic' (cubic in 2D) node.
            self.bnode = nm.array( [nodes[-1,:]], nm.int32 )
            self.bnode[:] = 1
            self.bf_gen_bubble = ap_bfgen_lagrange( v_coors, bf, self.bnode, \
                                                  var_set )

        # 01.09.2007
        def __call__( self, coors, node, iv = None, suppress_errors = False ):
            valb = self.bf_gen_bubble( coors, self.bnode[0], iv )
            if nm.sum( node ) >= 0:
                bf = self.bf_gen( coors, node, iv,
                                 suppress_errors = suppress_errors )
                return bf - valb / float( self.n_f )
            else:
                return valb
            
    return Generator( *args )

##
# 07.12.2005, c
def ap_bfgen_lagrange_tensor111( *args ):

    class Generator( Struct ):
        def __init__( self, v_coors, bf, nodes, var_set ):

            if bf.grad > 1:
                raise NotImplementedError

            dim = v_coors.shape[1]
            c_min = nm.amin( v_coors, 0 )
            c_max = nm.amax( v_coors, 0 )
##             print dim, c_min, c_max

            ##
            # Make a fake 1D gel.
            gd = Struct( n_vertex = 2, dim = 1 )
            gel = GeomElement( data = {'aux' : gd} )

            self.bf_gens = []
            for ii in range( dim ):
                vc = nm.array( [[c_min[ii]], [c_max[ii]]] )
                ao = nm.sum( nodes[0,2*ii:2*ii+2] )
##                 print vc, ao
                aux, nod1 = ap_nodegen_lagrange( ao, 'aux', None, gel, 0 )
##                 print nod1
##                 pause()
                bf_gen = ap_bfgen_lagrange( vc, bf, nod1, var_set )
                self.bf_gens.append( bf_gen )

        # 01.09.2007
        def __call__( self, coors, node, iv = None, suppress_errors = False ):

##             print coors, node, iv
            dim = coors.shape[1]
            if iv is None:
                val = 1.0
                for ii in range( dim ):
                    val *= self.bf_gens[ii]( coors[:,ii],
                                            node[2*ii:2*ii+2], iv,
                                            suppress_errors = suppress_errors )
                return val
            else:
                val = 1.0
                for ii in range( dim ):
                    if ii == iv:
                        val *= self.bf_gens[ii]( coors[:,ii],
                                                node[2*ii:2*ii+2], 0,
                                                suppress_errors = suppress_errors )
                    else:
                        val *= self.bf_gens[ii]( coors[:,ii],
                                                node[2*ii:2*ii+2], None,
                                                suppress_errors = suppress_errors )
                return val

            
    return Generator( *args )

##
# (volume ('v', 'b'), surface ('s*')) generator
ap_node_generators = {
    'Lagrange' : (ap_nodegen_lagrange, ap_nodegen_lagrange),
    'LagrangeP1B' : (ap_nodegen_lagrange_pxb, ap_nodegen_lagrange),
    'LagrangeP2B' : (ap_nodegen_lagrange_pxb, ap_nodegen_lagrange),
    'LagrangeTensor111' : (ap_nodegen_lagrange_tensor111,
                           ap_nodegen_lagrange_tensor111)
}

ap_bf_generators = {
    'Lagrange' : (ap_bfgen_lagrange, ap_bfgen_lagrange),
    'LagrangeP1B' : (ap_bfgen_lagrange_pxb, ap_bfgen_lagrange),
    'LagrangeP2B' : (ap_bfgen_lagrange_pxb, ap_bfgen_lagrange),
    'LagrangeTensor111' : (ap_bfgen_lagrange_tensor111,
                           ap_bfgen_lagrange_tensor111)
}
