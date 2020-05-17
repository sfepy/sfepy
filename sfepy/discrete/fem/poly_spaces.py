from __future__ import absolute_import
import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import find_subclasses, assert_, Struct
from sfepy.linalg import combine, insert_strided_axis
from six.moves import range
from functools import reduce

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
               1 : [[0],
                    [1]],
               0 : [[0]]}

def transform_basis(transform, bf):
    """
    Transform a basis `bf` using `transform` array of matrices.
    """
    if bf.ndim == 3:
        nbf = nm.einsum('cij,qdj->cqdi', transform, bf, order='C')

    elif bf.ndim == 4:
        if bf.shape[0] == 1:
            nbf = nm.einsum('cij,qdj->cqdi', transform, bf[0], order='C')

        else:
            nbf = nm.einsum('cij,cqdj->cqdi', transform, bf, order='C')

    # Note: the 2nd derivatives are not supported here.
    # Workaround for NumPy 1.14.0 - order is ignored(?)
    nbf = nm.ascontiguousarray(nbf)

    return nbf

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
                    nts[iseq] = [nt, ii]
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
        (0, 1) : 'simplex',
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
            # TODO how to get poly_spaces.py import LegendreSimplexPolySpace,
            #  LegendreTensorProductPolySpace from sfepy.discrete.dg.dg_poly_spaces
            #  to have them available here in globals? There is a cyclic import in
            #  the way.
            PolySpace._all = find_subclasses(globals(), [PolySpace])
        table = PolySpace._all

        key = '%s_%s' % (base, PolySpace.keys[(geometry.dim,
                                               geometry.n_vertex)])
        if (geometry.name == '1_2') and (key not in table):
            key = '%s_%s' % (base, 'tensor_product')

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

    def eval_base(self, coors, diff=0, ori=None, force_axis=False,
                  transform=None, suppress_errors=False, eps=1e-15):
        """
        Evaluate the basis or its first or second derivatives in points given
        by coordinates. The real work is done in _eval_base() implemented in
        subclasses.

        Note that the second derivative code is a work-in-progress and only
        `coors` and `transform` arguments are used.

        Parameters
        ----------
        coors : array_like
            The coordinates of points where the basis is evaluated. See Notes.
        diff : 0, 1 or 2
            If nonzero, return the given derivative.
        ori : array_like, optional
            Optional orientation of element facets for per element basis.
        force_axis : bool
            If True, force the resulting array shape to have one more axis even
            when `ori` is None.
        transform : array_like, optional
            The basis transform array.
        suppress_errors : bool
            If True, do not report points outside the reference domain.
        eps : float
            Accuracy for comparing coordinates.

        Returns
        -------
        base : array
            The basis (shape (n_coor, 1, n_base)) or its first derivative
            (shape (n_coor, dim, n_base)) or its second derivative (shape
            (n_coor, dim, dim, n_base)) evaluated in the given points. An
            additional axis is pre-pended of length n_cell, if `ori` is given,
            or of length 1, if `force_axis` is True.

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
            base = self._eval_base(coors, diff=diff, ori=ori,
                                   suppress_errors=suppress_errors,
                                   eps=eps)

            if (base.ndim == 3) and force_axis:
                base = base[None, ...]

            if not base.flags['C_CONTIGUOUS']:
                base = nm.ascontiguousarray(base)

        else: # Several point sets.
            if diff:
                bdim = self.geometry.dim
            else:
                bdim = 1

            base = nm.empty((coors.shape[0], coors.shape[1],
                             bdim, self.n_nod), dtype=nm.float64)

            for ii, _coors in enumerate(coors):
                base[ii] = self._eval_base(_coors, diff=diff, ori=ori,
                                           suppress_errors=suppress_errors,
                                           eps=eps)

        if transform is not None:
            base = transform_basis(transform, base)

        return base

    def get_mtx_i(self):
        return self.mtx_i

    def describe_nodes(self):
        return NodeDescription(self.nts, self.nodes)

class LagrangePolySpace(PolySpace):

    def create_context(self, cmesh, eps, check_errors, i_max, newton_eps,
                       tdim=None):
        from sfepy.discrete.fem.extmods.bases import CLagrangeContext

        ref_coors = self.geometry.coors

        if cmesh is not None:
            mesh_coors = cmesh.coors

            conn = cmesh.get_conn(cmesh.tdim, 0)
            mesh_conn = conn.indices.reshape(cmesh.n_el, -1).astype(nm.int32)

            if tdim is None:
                tdim = cmesh.tdim

        else:
            mesh_coors = mesh_conn = None

        if tdim is None:
            raise ValueError('supply either cmesh or tdim!')

        ctx = CLagrangeContext(order=self.order,
                               tdim=tdim,
                               nodes=self.nodes,
                               ref_coors=ref_coors,
                               mesh_coors=mesh_coors,
                               mesh_conn=mesh_conn,
                               mtx_i=self.get_mtx_i(),
                               eps=eps,
                               check_errors=check_errors,
                               i_max=i_max,
                               newton_eps=newton_eps)

        return ctx

    def _eval_base(self, coors, diff=0, ori=None,
                   suppress_errors=False, eps=1e-15):
        """
        See :func:`PolySpace.eval_base()`.
        """
        if diff == 2:
            base = self._eval_hessian(coors)

        else:
            base = self.eval_ctx.evaluate(coors, diff=diff,
                                          eps=eps,
                                          check_errors=not suppress_errors)
        return base

class LagrangeSimplexPolySpace(LagrangePolySpace):
    """Lagrange polynomial space on a simplex domain."""
    name = 'lagrange_simplex'

    def __init__(self, name, geometry, order, init_context=True):
        PolySpace.__init__(self, name, geometry, order)

        n_v = geometry.n_vertex

        mtx = nm.ones((n_v, n_v), nm.float64)
        mtx[0:n_v-1,:] = nm.transpose(geometry.coors)
        self.mtx_i = nm.ascontiguousarray(nla.inv(mtx))
        self.rhs = nm.ones((n_v,), nm.float64)

        self.nodes, self.nts, node_coors = self._define_nodes()
        self.node_coors = nm.ascontiguousarray(node_coors)
        self.n_nod = self.nodes.shape[0]

        if init_context:
            self.eval_ctx = self.create_context(None, 0, 1e-15, 100, 1e-8,
                                                tdim=n_v - 1)

        else:
            self.eval_ctx = None

    def _define_nodes(self):
        # Factorial.
        fac = lambda n : reduce(lambda a, b : a * (b + 1), range(n), 1)

        geometry = self.geometry
        n_v, dim = geometry.n_vertex, geometry.dim
        order = self.order

        n_nod = fac(order + dim) // (fac(order) * fac(dim))
        ## print n_nod, gd
        nodes = nm.zeros((n_nod, n_v), nm.int32)
        nts = nm.zeros((n_nod, 2), nm.int32)

        if order == 0:
            nts[0,:] = [3, 0]
            nodes[0,:] = nm.zeros((n_v,), nm.int32)

        else:
            iseq = 0

            # Vertex nodes.
            nts[0:n_v,0] = 0
            nts[0:n_v,1] = nm.arange(n_v, dtype = nm.int32)
            aux = order * nm.identity(n_v, dtype = nm.int32)
            nodes[iseq:iseq+n_v,:] = aux
            iseq += n_v

            if dim == 0:
                pass

            elif dim == 1:
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

    def _eval_hessian(self, coors):
        """
        Evaluate the second derivatives of the basis.
        """
        def get_bc(coor):
            rhs = nm.concatenate((coor, [1]))
            bc = nm.dot(self.mtx_i, rhs)

            return bc

        def get_val(bc, node, omit=[]):
            val = nm.ones(1, nm.float64)
            for i1 in range(bc.shape[0]):
                if i1 in omit: continue

                for i2 in range(node[i1]):
                    val *= (self.order * bc[i1] - i2) / (i2 + 1.0)

            return val

        def get_der(bc1, node1, omit=[]):
            val = nm.zeros(1, nm.float64)
            for i1 in range(node1):
                if i1 in omit: continue

                aux = nm.ones(1, nm.float64)
                for i2 in range(node1):
                    if (i1 == i2) or (i2 in omit): continue
                    aux *= (self.order * bc1 - i2) / (i2 + 1.0)

                val += aux * self.order / (i1 + 1.0)

            return val

        n_v = self.mtx_i.shape[0]
        dim = n_v - 1

        mi = self.mtx_i[:, :dim]
        bfgg = nm.zeros((coors.shape[0], dim, dim, self.n_nod),
                        dtype=nm.float64)

        for ic, coor in enumerate(coors):
            bc = get_bc(coor)

            for ii, node in enumerate(self.nodes):
                for ig1, bc1 in enumerate(bc): # 1. derivative w.r.t. bc1.
                    for ig2, bc2 in enumerate(bc): # 2. derivative w.r.t. bc2.
                        if ig1 == ig2:
                            val = get_val(bc, node, omit=[ig1])

                            vv = 0.0
                            for i1 in range(node[ig1]):
                                aux = get_der(bc2, node[ig2], omit=[i1])
                                vv += aux * self.order / (i1 + 1.0)

                            val *= vv

                        else:
                            val = get_val(bc, node, omit=[ig1, ig2])
                            val *= get_der(bc1, node[ig1])
                            val *= get_der(bc2, node[ig2])

                        bfgg[ic, :, :, ii] += val * mi[ig1] * mi[ig2][:, None]

        return bfgg

class LagrangeSimplexBPolySpace(LagrangeSimplexPolySpace):
    """Lagrange polynomial space with forced bubble function on a simplex
    domain."""
    name = 'lagrange_simplex_bubble'

    def __init__(self, name, geometry, order, init_context=True):
        LagrangeSimplexPolySpace.__init__(self, name, geometry, order,
                                          init_context=False)

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

        if init_context:
            self.eval_ctx = self.create_context(None, 0, 1e-15, 100, 1e-8,
                                                tdim=n_v - 1)

        else:
            self.eval_ctx = None

    def create_context(self, *args, **kwargs):
        ctx = LagrangePolySpace.create_context(self, *args, **kwargs)
        ctx.is_bubble = 1

        return ctx

class LagrangeTensorProductPolySpace(LagrangePolySpace):
    """Lagrange polynomial space on a tensor product domain."""
    name = 'lagrange_tensor_product'

    def __init__(self, name, geometry, order, init_context=True):
        PolySpace.__init__(self, name, geometry, order)

        g1d = Struct(n_vertex = 2,
                     dim = 1,
                     coors = self.bbox[:,0:1].copy())
        self.ps1d = LagrangeSimplexPolySpace('P_aux', g1d, order,
                                             init_context=False)

        self.nodes, self.nts, node_coors = self._define_nodes()
        self.node_coors = nm.ascontiguousarray(node_coors)
        self.n_nod = self.nodes.shape[0]

        if init_context:
            tdim = int(nm.sqrt(geometry.n_vertex))
            self.eval_ctx = self.create_context(None, 0, 1e-15, 100, 1e-8,
                                                tdim=tdim)

        else:
            self.eval_ctx = None

    def _define_nodes(self):
        geometry = self.geometry
        order = self.order

        n_v, dim = geometry.n_vertex, geometry.dim

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
            order * nm.identity( n_v, dtype = nm.int32 )
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
                    i1 = vertex_map[ii][0]
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

    def _eval_base_debug(self, coors, diff=False, ori=None,
                         suppress_errors=False, eps=1e-15):
        """Python version of eval_base()."""
        dim = self.geometry.dim

        ev = self.ps1d.eval_base

        if diff:
            base = nm.ones((coors.shape[0], dim, self.n_nod), dtype=nm.float64)

            for ii in range(dim):
                self.ps1d.nodes = self.nodes[:,2*ii:2*ii+2].copy()
                self.ps1d.n_nod = self.n_nod

                for iv in range(dim):
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

            for ii in range(dim):
                self.ps1d.nodes = self.nodes[:,2*ii:2*ii+2].copy()
                self.ps1d.n_nod = self.n_nod

                base *= ev(coors[:,ii:ii+1].copy(),
                           diff=diff,
                           suppress_errors=suppress_errors,
                           eps=eps)

        return base

    def _eval_hessian(self, coors):
        """
        Evaluate the second derivatives of the basis.
        """
        evh = self.ps1d.eval_base

        dim = self.geometry.dim
        bfgg = nm.zeros((coors.shape[0], dim, dim, self.n_nod),
                        dtype=nm.float64)

        v0s = []
        v1s = []
        v2s = []
        for ii in range(dim):
            self.ps1d.nodes = self.nodes[:,2*ii:2*ii+2].copy()
            self.ps1d.n_nod = self.n_nod
            ev = self.ps1d.create_context(None, 0, 1e-15, 100, 1e-8,
                                          tdim=1).evaluate

            v0s.append(ev(coors[:, ii:ii+1].copy())[:, 0, :])
            v1s.append(ev(coors[:, ii:ii+1].copy(), diff=1)[:, 0, :])
            v2s.append(evh(coors[:, ii:ii+1], diff=2)[:, 0, 0, :])

        for ir in range(dim):
            vv = v2s[ir] # Destroys v2s!
            for ik in range(dim):
                if ik == ir: continue
                vv *= v0s[ik]

            bfgg[:, ir, ir, :] = vv

            for ic in range(dim):
                if ic == ir: continue
                val = v1s[ir] * v1s[ic]
                for ik in range(dim):
                    if (ik == ir) or (ik == ic): continue
                    val *= v0s[ik]

                bfgg[:, ir, ic, :] += val

        return bfgg

    def get_mtx_i(self):
        return self.ps1d.mtx_i

class LobattoTensorProductPolySpace(PolySpace):
    """
    Hierarchical polynomial space using Lobatto functions.

    Each row of the `nodes` attribute defines indices of Lobatto functions that
    need to be multiplied together to evaluate the corresponding shape
    function. This defines the ordering of basis functions on the reference
    element.
    """
    name = 'lobatto_tensor_product'

    def __init__(self, name, geometry, order):
        PolySpace.__init__(self, name, geometry, order)

        aux = self._define_nodes()
        self.nodes, self.nts, node_coors, self.face_axes, self.sfnodes = aux
        self.node_coors = nm.ascontiguousarray(node_coors)
        self.n_nod = self.nodes.shape[0]

        aux = nm.where(self.nodes > 0, self.nodes, 1)
        self.node_orders = nm.prod(aux, axis=1)
        self.edge_indx = nm.where(self.nts[:, 0] == 1)[0]
        self.face_indx = nm.where(self.nts[:, 0] == 2)[0]

        self.face_axes_nodes = self._get_face_axes_nodes(self.face_axes)

    def _get_counts(self):
        order = self.order
        dim = self.geometry.dim

        n_nod = (order + 1) ** dim
        n_per_edge = (order - 1)
        n_per_face = (order - 1) ** (dim - 1)
        n_bubble = (order - 1) ** dim

        return n_nod, n_per_edge, n_per_face, n_bubble

    def _define_nodes(self):
        geometry = self.geometry
        order = self.order

        n_v, dim = geometry.n_vertex, geometry.dim

        n_nod, n_per_edge, n_per_face, n_bubble = self._get_counts()

        nodes = nm.zeros((n_nod, dim), nm.int32)
        nts = nm.zeros((n_nod, 2), nm.int32)

        # Vertex nodes.
        nts[0:n_v, 0] = 0
        nts[0:n_v, 1] = nm.arange(n_v, dtype=nm.int32)
        nodes[0:n_v] = nm.array(vertex_maps[dim], dtype=nm.int32)
        ii = n_v

        # Edge nodes.
        if (dim > 1) and (n_per_edge > 0):
            ik = nm.arange(2, order + 1, dtype=nm.int32)
            zo = nm.zeros((n_per_edge, 2), dtype=nm.int32)
            zo[:, 1] = 1
            for ie, edge in enumerate(geometry.edges):
                n1, n2 = nodes[edge]
                ifix = nm.where(n1 == n2)[0]
                irun = nm.where(n1 != n2)[0][0]
                ic = n1[ifix]

                nodes[ii:ii + n_per_edge, ifix] = zo[:, ic]
                nodes[ii:ii + n_per_edge, irun] = ik
                nts[ii:ii + n_per_edge] = [[1, ie]]
                ii += n_per_edge

        # 3D face nodes.
        face_axes = []
        sfnodes = None
        if (dim == 3) and (n_per_face > 0):
            n_face = len(geometry.faces)
            sfnodes = nm.zeros((n_per_face * n_face, dim), nm.int32)
            ii0 = ii

            ik = nm.arange(2, order + 1, dtype=nm.int32)
            zo = nm.zeros((n_per_face, 2), dtype=nm.int32)
            zo[:, 1] = 1

            for ifa, face in enumerate(geometry.faces):
                ns = nodes[face]

                diff = nm.diff(ns, axis=0)
                asum = nm.abs(diff).sum(axis=0)
                ifix = nm.where(asum == 0)[0][0]
                ic = ns[0, ifix]
                irun1 = nm.where(asum == 2)[0][0]
                irun2 = nm.where(asum == 1)[0][0]

                iy, ix = nm.meshgrid(ik, ik)

                nodes[ii:ii + n_per_face, ifix] = zo[:, ic]
                nodes[ii:ii + n_per_face, irun1] = ix.ravel()
                nodes[ii:ii + n_per_face, irun2] = iy.ravel()
                nts[ii:ii + n_per_face] = [[2, ifa]]

                ij = ii - ii0
                sfnodes[ij:ij + n_per_face, ifix] = zo[:, ic]
                sfnodes[ij:ij + n_per_face, irun1] = iy.ravel()
                sfnodes[ij:ij + n_per_face, irun2] = ix.ravel()

                face_axes.append([irun1, irun2])

                ii += n_per_face

        face_axes = nm.array(face_axes)

        # Bubble nodes.
        if n_bubble > 0:
            ik = nm.arange(2, order + 1, dtype=nm.int32)
            nodes[ii:] = nm.array([aux for aux in combine([ik] * dim)])
            nts[ii:ii + n_bubble] = [[3, 0]]
            ii += n_bubble

        assert_(ii == n_nod)

        # Coordinates of the "nodes". All nodes on a facet have the same
        # coordinates - the centre of the facet.
        c_min, c_max = self.bbox[:, 0]

        node_coors = nm.zeros(nodes.shape, dtype=nm.float64)
        node_coors[:n_v] = nodes[:n_v]

        if (dim > 1) and (n_per_edge > 0):
            ie = nm.where(nts[:, 0] == 1)[0]
            node_coors[ie] = node_coors[geometry.edges[nts[ie, 1]]].mean(1)

        if (dim == 3) and (n_per_face > 0):
            ifa = nm.where(nts[:, 0] == 2)[0]
            node_coors[ifa] = node_coors[geometry.faces[nts[ifa, 1]]].mean(1)

        if n_bubble > 0:
            ib = nm.where(nts[:, 0] == 3)[0]
            node_coors[ib] = node_coors[geometry.conn].mean(0)

        return nodes, nts, node_coors, face_axes, sfnodes

    def _get_face_axes_nodes(self, face_axes):
        if not len(face_axes): return None

        nodes = self.nodes[self.face_indx]
        n_per_face = self._get_counts()[2]
        anodes = nm.tile(nodes[:n_per_face, face_axes[0]], (6, 1))

        return anodes

    def _eval_base(self, coors, diff=False, ori=None,
                   suppress_errors=False, eps=1e-15):
        """
        See PolySpace.eval_base().
        """
        from .extmods.lobatto_bases import eval_lobatto_tensor_product as ev
        c_min, c_max = self.bbox[:, 0]

        base = ev(coors, self.nodes, c_min, c_max, self.order, diff)

        if ori is not None:
            ebase = nm.tile(base, (ori.shape[0], 1, 1, 1))

            if self.edge_indx.shape[0]:
                # Orient edge functions.
                ie, ii = nm.where(ori[:, self.edge_indx] == 1)
                ii = self.edge_indx[ii]
                ebase[ie, :, :, ii] *= -1.0

            if self.face_indx.shape[0]:
                # Orient face functions.
                fori = ori[:, self.face_indx]

                # ... normal axis order
                ie, ii = nm.where((fori == 1) | (fori == 2))
                ii = self.face_indx[ii]
                ebase[ie, :, :, ii] *= -1.0

                # ... swapped axis order
                sbase = ev(coors, self.sfnodes, c_min, c_max, self.order, diff)
                sbase = insert_strided_axis(sbase, 0, ori.shape[0])

                # ...overwrite with swapped axes basis.
                ie, ii = nm.where(fori >= 4)
                ii2 = self.face_indx[ii]
                ebase[ie, :, :, ii2] = sbase[ie, :, :, ii]

                # ...deal with orientation.
                ie, ii = nm.where((fori == 5) | (fori == 6))
                ii = self.face_indx[ii]
                ebase[ie, :, :, ii] *= -1.0

            base = ebase

        return base
