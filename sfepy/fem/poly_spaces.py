import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import find_subclasses, Struct

class LagrangeNodes(Struct):
    """Helper class for defining nodes of Lagrange elements."""

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
                aux = [int(round(tmp)) for tmp in delta * (c1 * n1 + c2 * n2)]
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
                    aux = [int(round(tmp)) for tmp
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
                    aux = [int(round(tmp)) for tmp
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
                aux = [int(round(tmp)) for tmp in delta * (c1 * n1 + c2 * n2)]
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
                    aux = [int(round(tmp)) for tmp
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
                    aux = [int(round(tmp)) for tmp
                           in delta * (c1 * c2 * c3 * n1 + c4 * c2 * c3 * n2
                                       + c5 * c4 * c3 * n3 + c1 * c3 * c5 * n4
                                       + c1 * c2 * c6 * n5 + c4 * c2 * c6 * n6
                                       + c5 * c4 * c6 * n7 + c1 * c6 * c5 * n8)]
                    nodes[iseq,:] = aux
                    iseq += 1
        return iseq

class NodeDescription(Struct):
    """
    Describe FE nodes defined on different parts of a reference element.
    """

    def _describe_facets(self, ii):
        nts = self.node_types[ii]
        ik = nm.where(nts[1:,1] > nts[:-1,1])[0]

        if len(ik) == 0:
            ifacets = None
            n_dof = 0

        else:
            ii = ii.astype(nm.int32)

            ik = nm.r_[0, ik + 1, nts.shape[0]]
            ifacets = [ii[ik[ir] : ik[ir+1]] for ir in range(len(ik) - 1)]
            n_dof = len(ifacets[0])

        return ifacets, n_dof

    def _describe_other(self, ii):
        if len(ii):
            return ii, len(ii)

        else:
            return None, 0

    def _get_facet_nodes(self, ifacets, nodes):
        if ifacets is None:
            return None

        else:
            return [nodes[ii] for ii in ifacets]

    def _get_nodes(self, ii, nodes):
        if ii is None:
            return None

        else:
            return nodes[ii]

    def __init__(self, node_types, nodes):
        self.node_types = node_types

        # Vertex nodes.
        ii = nm.where(node_types[:,0] == 0)[0]
        self.vertex, self.n_vertex_nod = self._describe_other(ii)
        self.vertex_nodes = self._get_nodes(self.vertex, nodes)

        # Edge nodes.
        ii = nm.where(node_types[:,0] == 1)[0]
        self.edge, self.n_edge_nod = self._describe_facets(ii)
        self.edge_nodes = self._get_facet_nodes(self.edge, nodes)

        # Face nodes.
        ii = nm.where(node_types[:,0] == 2)[0]
        self.face, self.n_face_nod = self._describe_facets(ii)
        self.face_nodes = self._get_facet_nodes(self.face, nodes)

        # Bubble nodes.
        ii = nm.where(node_types[:,0] == 3)[0]
        self.bubble, self.n_bubble_nod = self._describe_other(ii)
        self.bubble_nodes = self._get_nodes(self.bubble, nodes)

    def has_extra_nodes(self):
        """
        Return True if the element has some edge, face or bubble nodes.
        """
        return (self.n_edge_nod + self.n_face_nod + self.n_bubble_nod) > 0

class PolySpace(Struct):
    """Abstract polynomial space class."""
    _all = None

    keys = {
        (1, 2) : 'simplex',
        (2, 3) : 'simplex',
        (3, 4) : 'simplex',
        (2, 4) : 'tensor_product',
        (3, 8) : 'tensor_product',
    }

    @staticmethod
    def any_from_args(name, geometry, order, base='lagrange',
                      force_bubble=False):
        """
        Construct a particular polynomial space classes according to the
        arguments passed in.
        """
        if name is None:
            name = PolySpace.suggest_name(geometry, order, base, force_bubble)

        if PolySpace._all is None:
            PolySpace._all = find_subclasses(globals(),
                                                 [PolySpace])
        table = PolySpace._all

        key = '%s_%s' % (base, PolySpace.keys[(geometry.dim,
                                               geometry.n_vertex)])
        if force_bubble:
            key += '_bubble'

        return table[key](name, geometry, order)

    @staticmethod
    def suggest_name(geometry, order, base='lagrange',
                     force_bubble=False):
        """
        Suggest the polynomial space name given its constructor parameters.
        """
        aux = geometry.get_interpolation_name()[:-1]
        if force_bubble:
            return aux + ('%dB' % order)
        else:
            return aux + ('%d' % order)

    def __init__(self, name, geometry, order):
        self.name = name
        self.geometry = geometry
        self.order = order

        self.bbox = nm.vstack((geometry.coors.min(0), geometry.coors.max(0)))

    def eval_base(self, coors, diff=False,
                  suppress_errors=False, eps=1e-15):
        """
        Evaluate the basis in points given by coordinates. The real work is
        done in _eval_base() implemented in subclasses.

        Parameters
        ----------
        coors : array_like
            The coordinates of points where the basis is evaluated. See Notes.
        diff : bool
            If True, return the first derivative.
        suppress_errors : bool
            If True, do not report points outside the reference domain.
        eps : float
            Accuracy for comparing coordinates.

        Returns
        -------
        base : array
            The basis (shape (n_coor, 1, n_base)) or its derivative (shape
            (n_coor, dim, n_base)) evaluated in the given points.

        Notes
        -----
        If coors.ndim == 3, several point sets are assumed, with equal number
        of points in each of them. This is the case, for example, of the
        values of the volume base functions on the element facets. The indexing
        (of bf_b(g)) is then (ifa,iqp,:,n_ep), so that the facet can be set in
        C using FMF_SetCell.
        """
        coors = nm.asarray(coors)
        if not coors.ndim in (2, 3):
            raise ValueError('coordinates must have 2 or 3 dimensions! (%d)'
                             % coors.ndim)

        if (coors.ndim == 2):
            base = self._eval_base(coors, diff=diff,
                                   suppress_errors=suppress_errors,
                                   eps=eps)

        else: # Several point sets.
            if diff:
                bdim = self.geometry.dim
            else:
                bdim = 1

            base = nm.empty((coors.shape[0], coors.shape[1],
                             bdim, self.n_nod), dtype=nm.float64)

            for ii, _coors in enumerate(coors):
                base[ii] = self._eval_base(_coors, diff=diff,
                                           suppress_errors=suppress_errors,
                                           eps=eps)

        return base

    def get_mtx_i(self):
        return self.mtx_i

    def describe_nodes(self):
        return NodeDescription(self.nts, self.nodes)

class LagrangeSimplexPolySpace(PolySpace):
    """Lagrange polynomial space on a simplex domain."""
    name = 'lagrange_simplex'

    def __init__(self, name, geometry, order):
        PolySpace.__init__(self, name, geometry, order)

        n_v = geometry.n_vertex

        mtx = nm.ones((n_v, n_v), nm.float64)
        mtx[0:n_v-1,:] = nm.transpose(geometry.coors)
        self.mtx_i = nm.ascontiguousarray(nla.inv(mtx))
        self.rhs = nm.ones((n_v,), nm.float64)

        self.nodes, self.nts, node_coors = self._define_nodes()
        self.node_coors = nm.ascontiguousarray(node_coors)
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

    def _eval_base(self, coors, diff=False,
                   suppress_errors=False, eps=1e-15):
        """See PolySpace.eval_base()."""
        from extmods.bases import eval_lagrange_simplex

        base = eval_lagrange_simplex(coors, self.mtx_i, self.nodes,
                                     self.order, diff,
                                     eps, not suppress_errors)

        return base

class LagrangeSimplexBPolySpace(LagrangeSimplexPolySpace):
    """Lagrange polynomial space with forced bubble function on a simplex
    domain."""
    name = 'lagrange_simplex_bubble'

    def __init__(self, name, geometry, order):
        LagrangeSimplexPolySpace.__init__(self, name, geometry, order)

        nodes, nts, node_coors = self.nodes, self.nts, self.node_coors

        shape = [nts.shape[0] + 1, 2]
        nts = nm.resize(nts, shape)
        nts[-1,:] = [3, 0]

        shape = [nodes.shape[0] + 1, nodes.shape[1]]
        nodes = nm.resize(nodes, shape)
        # Make a 'hypercubic' (cubic in 2D) node.
        nodes[-1,:] = 1

        n_v = self.geometry.n_vertex
        tmp = nm.ones((n_v,), nm.int32)

        node_coors = nm.vstack((node_coors,
                                nm.dot(tmp, self.geometry.coors) / n_v))

        self.nodes, self.nts = nodes, nts
        self.node_coors = nm.ascontiguousarray(node_coors)

        self.bnode = nodes[-1:,:]

        self.n_nod = self.nodes.shape[0]

    def _eval_base(self, coors, diff=False,
                   suppress_errors=False, eps=1e-15):
        """See PolySpace.eval_base()."""
        from extmods.bases import eval_lagrange_simplex

        base = eval_lagrange_simplex(coors, self.mtx_i, self.nodes[:-1],
                                     self.order, diff,
                                     eps, not suppress_errors)
        bubble = eval_lagrange_simplex(coors, self.mtx_i, self.bnode,
                                       int(self.bnode.sum()), diff,
                                       eps, not suppress_errors)

        base -= bubble / (self.n_nod - 1)
        base = nm.ascontiguousarray(nm.dstack((base, bubble)))

        return base

class LagrangeTensorProductPolySpace(PolySpace):
    """Lagrange polynomial space on a tensor product domain."""
    name = 'lagrange_tensor_product'

    def __init__(self, name, geometry, order):
        PolySpace.__init__(self, name, geometry, order)

        g1d = Struct(n_vertex = 2,
                     dim = 1,
                     coors = self.bbox[:,0:1])
        self.ps1d = LagrangeSimplexPolySpace('P_aux', g1d, order)

        self.nodes, self.nts, node_coors = self._define_nodes()
        self.node_coors = nm.ascontiguousarray(node_coors)
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

    def _eval_base(self, coors, diff=False,
                   suppress_errors=False, eps=1e-15):
        """See PolySpace.eval_base()."""
        from extmods.bases import eval_lagrange_tensor_product as ev
        ps1d = self.ps1d

        base = ev(coors, ps1d.mtx_i, self.nodes, self.order, diff,
                  eps, not suppress_errors)

        return base

    def get_mtx_i(self):
        return self.ps1d.mtx_i
