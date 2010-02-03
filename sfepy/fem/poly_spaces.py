from sfepy.base.base import *
from sfepy.base.la import cycle

class LagrangeNodes(Struct):

    @staticmethod
    def append_edges(nodes, nts, iseq, nt, edges, order):
        delta = 1.0 / float(order)

        for ii, edge in enumerate(edges):
            n1 = nodes[edge[0],:].copy()
            n2 = nodes[edge[1],:].copy()
            for ie in range(order - 1):
                c2 = ie + 1
                c1 = order - c2
                nts[iseq] = [nt, ii]
                aux = [int(tmp) for tmp in delta * (c1 * n1 + c2 * n2)]
                nodes[iseq,:] = aux
                iseq += 1
        return iseq

    @staticmethod
    def append_faces(nodes, nts, iseq, nt, faces, order):
        delta = 1.0 / float(order)

        for ii, face in enumerate(faces):
            n1 = nodes[face[0],:].copy()
            n2 = nodes[face[1],:].copy()
            n3 = nodes[face[2],:].copy()
            for i1 in range(order - 2):
                for i2 in range(order - 2 - i1):
                    c3 = i1 + 1
                    c2 = i2 + 1
                    c1 = order - c3 - c2
                    nts[iseq] = [nt, ii]
                    aux = [int(tmp) for tmp
                           in delta * (c1 * n1 + c2 * n2 + c3 * n3)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

    @staticmethod
    def append_bubbles(nodes, nts, iseq, nt, order):
        delta = 1.0 / float(order)

        n1 = nodes[0,:].copy()
        n2 = nodes[1,:].copy()
        n3 = nodes[2,:].copy()
        n4 = nodes[3,:].copy()
        for i1 in range(order - 3):
            for i2 in range(order - 3):
                for i3 in range(order - 3 - i1 - i2):
                    c4 = i1 + 1
                    c3 = i2 + 1
                    c2 = i3 + 1
                    c1 = order - c4 - c3 - c2
                    nts[iseq] = [nt, 0]
                    aux = [int(tmp) for tmp
                           in delta * (c1 * n1 + c2 * n2 + c3 * n3 + c4 * n4)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

    @staticmethod
    def append_tp_edges(nodes, nts, iseq, nt, edges, ao):
        delta = 1.0 / float(ao)
        for ii, edge in enumerate(edges):
            n1 = nodes[edge[0],:].copy()
            n2 = nodes[edge[1],:].copy()
            for ie in range(ao - 1):
                c2 = ie + 1
                c1 = ao - c2
                nts[iseq] = [nt, ii]
                aux = [int(tmp) for tmp in delta * (c1 * n1 + c2 * n2)]
                nodes[iseq,:] = aux
                iseq += 1
        return iseq

    @staticmethod
    def append_tp_faces(nodes, nts, iseq, nt, faces, ao):
        delta = 1.0 / (float(ao) ** 2)
        for ii, face in enumerate( faces ):
            n1 = nodes[face[0],:].copy()
            n2 = nodes[face[1],:].copy()
            n3 = nodes[face[2],:].copy()
            n4 = nodes[face[3],:].copy()
            for i1 in range(ao - 1):
                for i2 in range(ao - 1):
                    c4 = i1 + 1
                    c3 = i2 + 1
                    c2 = ao - c4
                    c1 = ao - c3
                    nts[iseq] = [2, ii]
                    aux = [int(tmp) for tmp
                           in delta * (c1 * c2 * n1 + c2 * c3 * n2
                                       + c3 * c4 * n3 + c4 * c1 * n4)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

    @staticmethod
    def append_tp_bubbles(nodes, nts, iseq, nt, ao):
        delta = 1.0 / (float(ao) ** 3)
        n1 = nodes[0,:].copy()
        n2 = nodes[1,:].copy()
        n3 = nodes[2,:].copy()
        n4 = nodes[3,:].copy()
        n5 = nodes[4,:].copy()
        n6 = nodes[5,:].copy()
        n7 = nodes[6,:].copy()
        n8 = nodes[7,:].copy()
        for i1 in range(ao - 1):
            for i2 in range(ao - 1):
                for i3 in range(ao - 1):
                    c6 = i1 + 1
                    c5 = i2 + 1
                    c4 = i3 + 1
                    c3 = ao - c6
                    c2 = ao - c5
                    c1 = ao - c4
                    nts[iseq] = [nt, 0]
                    aux = [int(tmp) for tmp
                           in delta * (c1 * c2 * c3 * n1 + c4 * c2 * c3 * n2
                                       + c5 * c4 * c3 * n3 + c1 * c3 * c5 * n4
                                       + c1 * c2 * c6 * n5 + c4 * c2 * c6 * n6
                                       + c5 * c4 * c6 * n7 + c1 * c6 * c5 * n8)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

                    
class PolySpace(Struct):
    _all = None

    keys = {
        (1, 2) : 'simplex',
        (3, 4) : 'simplex',
        (3, 4) : 'simplex',
        (2, 4) : 'tensor_product',
        (3, 8) : 'tensor_product',
    }

    @staticmethod
    def any_from_args(name, geometry, order, base='lagrange',
                      force_bubble=False):
        if PolySpace._all is None:
            PolySpace._all = find_subclasses(globals(),
                                                 [PolySpace])
        table = PolySpace._all

        key = '%s_%s' % (base, PolySpace.keys[(geometry.dim,
                                            geometry.n_vertex)])
        if force_bubble:
            key += '_bubble'

        return table[key](name, geometry, order)

    def __init__(self, name, geometry, order):
        self.name = name
        self.geometry = geometry
        self.order = order

        self.bbox = nm.vstack((geometry.coors.min(0), geometry.coors.max(0)))

    def eval_base(coors, diff=None,
                  suppress_errors=False, eps=1e-15):
        """Implemented in subclasses."""

class LagrangeSimplexPolySpace(PolySpace):
    name = 'lagrange_simplex'

    def __init__(self, name, geometry, order):
        PolySpace.__init__(self, name, geometry, order)

        n_v = geometry.n_vertex

        mtx = nm.ones((n_v, n_v), nm.float64)
        mtx[0:n_v-1,:] = nm.transpose(geometry.coors)
        self.mtx_i = nla.inv(mtx)
        self.rhs = nm.ones((n_v,), nm.float64)

        self.nodes, self.nts, self.node_coors = self._define_nodes()
        self.n_nod = self.nodes.shape[0]

    def _define_nodes(self):
        # Factorial.
        fac = lambda n : reduce(lambda a, b : a * (b + 1), xrange(n), 1)

        geometry = self.geometry
        n_v, dim = geometry.n_vertex, geometry.dim
        order = self.order

        n_nod = fac(order + dim) / (fac(order) * fac(dim))
        ## print n_nod, gd
        nodes = nm.zeros((n_nod, n_v), nm.int32)
        nts = nm.zeros((n_nod, 2), nm.int32)

        if order == 0:
            nts[0,:] = [3, 0]
            nodes[0,:] = nm.zeros((n_v,), nm.int32)

        else:
            iseq = 0
            delta = 1.0 / float(order)

            # Vertex nodes.
            nts[0:n_v,0] = 0
            nts[0:n_v,1] = nm.arange(n_v, dtype = nm.int32)
            aux = order * nm.identity(n_v, dtype = nm.int32)
            nodes[iseq:iseq+n_v,:] = aux
            iseq += n_v

            if dim == 1:
                iseq = LagrangeNodes.append_edges(nodes, nts, iseq, 3,
                                                  [[0, 1]], order)
            elif dim == 2:
                iseq = LagrangeNodes.append_edges(nodes, nts, iseq, 1,
                                                  geometry.edges, order)
                iseq = LagrangeNodes.append_faces(nodes, nts, iseq, 3,
                                                  [[0, 1, 2]], order)
            elif dim == 3:
                iseq = LagrangeNodes.append_edges(nodes, nts, iseq, 1,
                                                  geometry.edges, order)
                iseq = LagrangeNodes.append_faces(nodes, nts, iseq, 2,
                                                  geometry.faces, order)
                iseq = LagrangeNodes.append_bubbles(nodes, nts, iseq, 3,
                                                    order)
            else:
                raise NotImplementedError

        ## print nm.concatenate((nts, nodes), 1)

        # Check orders.
        orders = nm.sum(nodes, 1)
        if not nm.all(orders == order):
            raise AssertionError('wrong orders! (%d == all of %s)'
                                 % (order, orders))

        # Coordinates of the nodes.
        if order == 0:
            tmp = nm.ones((n_nod, n_v), nm.int32)
            node_coors = nm.dot(tmp, geometry.coors) / n_v

        else:
            node_coors = nm.dot(nodes, geometry.coors) / order

        return nodes, nts, node_coors

    def eval_base(self, coors, diff=False,
                  suppress_errors=False, eps=1e-15):
        from extmods.fem import eval_lagrange_simplex

        if diff:
            bdim = self.geometry.dim
        else:
            bdim = 1
        
        base = nm.empty((coors.shape[0], bdim, self.n_nod), dtype=nm.float64)

        # Work array.
        bc = nm.zeros((self.geometry.n_vertex, coors.shape[0]), nm.float64)

        eval_lagrange_simplex(base, coors, self.nodes, self.order, diff,
                              self.mtx_i, bc, suppress_errors, eps)

        return base

class LagrangeSimplexBPolySpace(LagrangeSimplexPolySpace):
    name = 'lagrange_simplex_bubble'

    def __init__(self, name, geometry, order):
        LagrangeSimplexPolySpace.__init__(self, name, geometry, order)

        nodes, nts, node_coors = self.nodes, self.nts, self.node_coors

        shape = [nts.shape[0] + 1, 2]
        nts = nm.resize(nts, shape)
        nts[-1,:] = [3, 0]

        shape = [nodes.shape[0] + 1, nodes.shape[1]]
        nodes = nm.resize(nodes, shape)
        nodes[-1,0] = -1
        nodes[-1,1:] = 0

        n_v = self.geometry.n_vertex
        tmp = nm.ones((n_v,), nm.int32)

        node_coors = nm.vstack((node_coors,
                                nm.dot(tmp, self.geometry.coors) / n_v))

        self.nodes, self.nts, self.node_coors = nodes, nts, node_coors
        self.n_nod = self.nodes.shape[0]

class LagrangeTensorProductPolySpace(PolySpace):
    name = 'lagrange_tensor_product'

    def __init__(self, name, geometry, order):
        PolySpace.__init__(self, name, geometry, order)

        g1d = Struct(n_vertex = 2,
                     dim = 1,
                     coors = self.bbox[:,0:1])
        self.ps1d = LagrangeSimplexPolySpace('P_aux', g1d, order)

        self.nodes, self.nts, self.node_coors = self._define_nodes()
        self.n_nod = self.nodes.shape[0]

    def _define_nodes(self):
        geometry = self.geometry
        order = self.order

        n_v, dim = geometry.n_vertex, geometry.dim

        # Requires fixed vertex numbering!
        vertex_maps = {3 : [[0, 0, 0],
                            [1, 0, 0],
                            [1, 1, 0],
                            [0, 1, 0],
                            [0, 0, 1],
                            [1, 0, 1],
                            [1, 1, 1],
                            [0, 1, 1]],
                       2 : [[0, 0],
                            [1, 0],
                            [1, 1],
                            [0, 1]],
                       1 : [[0,
                             1]]}
        vertex_map = order * nm.array(vertex_maps[dim], dtype=nm.int32)

        n_nod = (order + 1) ** dim
        nodes = nm.zeros((n_nod, 2 * dim), nm.int32)
        nts = nm.zeros((n_nod, 2), nm.int32)

        if order == 0:
            nts[0,:] = [3, 0]
            nodes[0,:] = nm.zeros((n_nod,), nm.int32)

        else:
            iseq = 0

            # Vertex nodes.
            nts[0:n_v,0] = 0
            nts[0:n_v,1] = nm.arange( n_v, dtype = nm.int32 )
            aux = order * nm.identity( n_v, dtype = nm.int32 )
            if dim == 3:
                for ii in range(n_v):
                    i1, i2, i3 = vertex_map[ii]
                    nodes[iseq,:] = [order - i1, i1,
                                     order - i2, i2,
                                     order - i3, i3]
                    iseq += 1
            elif dim == 2:
                for ii in range(n_v):
                    i1, i2 = vertex_map[ii]
                    nodes[iseq,:] = [order - i1, i1, order - i2, i2]
                    iseq += 1
            else:
                for ii in range(n_v):
                    i1 = vertex_map[ii]
                    nodes[iseq,:] = [order - i1, i1]
                    iseq += 1

            if dim == 1:
                iseq = LagrangeNodes.append_tp_edges(nodes, nts, iseq, 3,
                                                     [[0, 1]], order)
            elif dim == 2:
                iseq = LagrangeNodes.append_tp_edges(nodes, nts, iseq, 1,
                                                     geometry.edges, order)
                iseq = LagrangeNodes.append_tp_faces(nodes, nts, iseq, 3,
                                                     [[0, 1, 2, 3]], order)
            elif dim == 3:
                iseq = LagrangeNodes.append_tp_edges(nodes, nts, iseq, 1,
                                                     geometry.edges, order)
                iseq = LagrangeNodes.append_tp_faces(nodes, nts, iseq, 2,
                                                     geometry.faces, order)
                iseq = LagrangeNodes.append_tp_bubbles(nodes, nts, iseq, 3,
                                                       order)
            else:
                raise NotImplementedError

        # Check orders.
        orders = nm.sum(nodes, 1)
        if not nm.all(orders == order * dim):
            raise AssertionError('wrong orders! (%d == all of %s)'
                                 % (order * dim, orders))

        # Coordinates of the nodes.
        if order == 0:
            tmp = nm.ones((n_nod, n_v), nm.int32)
            node_coors = nm.dot(tmp, geometry.coors) / n_v

        else:
            c_min, c_max = self.bbox[:,0]
            
            cr = nm.arange(2 * dim)
            node_coors = (nodes[:,cr[::2]] * c_min
                          + nodes[:,cr[1::2]] * c_max) / order

        return nodes, nts, node_coors

    def _eval_base_debug(self, coors, diff=False,
                         suppress_errors=False, eps=1e-15):
        """Python version of eval_base()."""
        dim = self.geometry.dim

        ev = self.ps1d.eval_base
        
        if diff:
            base = nm.ones((coors.shape[0], dim, self.n_nod), dtype=nm.float64)

            for ii in xrange(dim):
                self.ps1d.nodes = self.nodes[:,2*ii:2*ii+2].copy()
                self.ps1d.n_nod = self.n_nod

                for iv in xrange(dim):
                    if ii == iv:
                        base[:,iv:iv+1,:] *= ev(coors[:,ii:ii+1].copy(),
                                                diff=True,
                                                suppress_errors=suppress_errors,
                                                eps=eps)
                
                    else:
                        base[:,iv:iv+1,:] *= ev(coors[:,ii:ii+1].copy(),
                                                diff=False,
                                                suppress_errors=suppress_errors,
                                                eps=eps)

        else:
            base = nm.ones((coors.shape[0], 1, self.n_nod), dtype=nm.float64)

            for ii in xrange(dim):
                self.ps1d.nodes = self.nodes[:,2*ii:2*ii+2].copy()
                self.ps1d.n_nod = self.n_nod
                
                base *= ev(coors[:,ii:ii+1].copy(),
                           diff=diff,
                           suppress_errors=suppress_errors,
                           eps=eps)

        return base

    def eval_base(self, coors, diff=False,
                  suppress_errors=False, eps=1e-15):
        from extmods.fem import eval_lagrange_tensor_product as ev
        ps1d = self.ps1d

        dim = self.geometry.dim
        if diff:
            bdim = dim
        else:
            bdim = 1

        base = nm.empty((coors.shape[0], bdim, self.n_nod), dtype=nm.float64)

        # Work arrays.
        bc = nm.zeros((ps1d.geometry.n_vertex, coors.shape[0]), nm.float64)
        base1d = nm.ones((coors.shape[0], 1, self.n_nod), dtype=nm.float64)

        ev(base, coors, self.nodes, self.order, diff, ps1d.mtx_i, bc, base1d,
           suppress_errors, eps)
                
        return base
