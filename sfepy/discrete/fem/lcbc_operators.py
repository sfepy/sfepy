"""
Operators for enforcing linear combination boundary conditions in nodal FEM
setting.
"""
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import output, assert_, Container, Struct
from sfepy.discrete.common.dof_info import DofInfo, expand_nodes_to_equations
from sfepy.discrete.fem.utils import (compute_nodal_normals,
                                      compute_nodal_edge_dirs)

class LCBCOperator(Struct):
    """
    Base class for LCBC operators.
    """

    def treat_pbcs(self, dofs, master):
        """
        Treat dofs with periodic BC.
        """
        master = nm.intersect1d(dofs, master)
        if len(master) == 0: return

        # Remove rows with master DOFs.
        remove = nm.searchsorted(nm.sort(dofs), master)
        keep = nm.setdiff1d(nm.arange(len(dofs)), remove)

        mtx = self.mtx[keep]

        # Remove empty columns, update new DOF count.
        mtx = mtx.tocsc()
        indptr = nm.unique(mtx.indptr)
        self.mtx = sp.csc_matrix((mtx.data, mtx.indices, indptr),
                                 shape=(mtx.shape[0], indptr.shape[0] - 1))
        self.n_dof = self.mtx.shape[1]

class RigidOperator(LCBCOperator):
    """
    Transformation matrix operator for rigid LCBCs.
    """

    def __init__(self, name, nodes, field, dof_names, all_dof_names):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        coors = field.get_coor(nodes)
        n_nod, dim = coors.shape

        mtx_e = nm.tile(nm.eye(dim, dtype=nm.float64), (n_nod, 1))

        if dim == 2:
            mtx_r = nm.empty((dim * n_nod, 1), dtype=nm.float64)
            mtx_r[0::dim,0] = -coors[:,1]
            mtx_r[1::dim,0] = coors[:,0]
            n_rigid_dof = 3

        elif dim == 3:
            mtx_r = nm.zeros((dim * n_nod, dim), dtype=nm.float64)
            mtx_r[0::dim,1] = coors[:,2]
            mtx_r[0::dim,2] = -coors[:,1]
            mtx_r[1::dim,0] = -coors[:,2]
            mtx_r[1::dim,2] = coors[:,0]
            mtx_r[2::dim,0] = coors[:,1]
            mtx_r[2::dim,1] = -coors[:,0]
            n_rigid_dof = 6

        else:
            msg = 'dimension in [2, 3]: %d' % dim
            raise ValueError(msg)

        self.n_dof = n_rigid_dof
        self.mtx = nm.hstack((mtx_r, mtx_e))

        # Strip unconstrained dofs.
        aux = dim * nm.arange(n_nod)
        indx = [aux + all_dof_names.index(dof) for dof in dof_names]
        indx = nm.array(indx).T.ravel()

        self.mtx = self.mtx[indx]

def _save_vectors(filename, vectors, region, mesh, data_name):
    """
    Save vectors defined in region nodes as vector data in mesh vertices.
    """
    nv = nm.zeros_like(mesh.coors)
    nmax = region.vertices.shape[0]
    nv[region.vertices] = vectors[:nmax]

    out = {data_name : Struct(name='output_data', mode='vertex', data=nv)}
    mesh.write(filename, out=out, io='auto')

class NoPenetrationOperator(LCBCOperator):
    """
    Transformation matrix operator for no-penetration LCBCs.
    """
    def __init__(self, name, nodes, region, field, dof_names, filename=None):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        dim = region.dim
        assert_(len(dof_names) == dim)

        normals = compute_nodal_normals(nodes, region, field)

        if filename is not None:
            _save_vectors(filename, normals, region, field.domain.mesh, 'n')

        ii = nm.abs(normals).argmax(1)
        n_nod, dim = normals.shape

        irs = set(range(dim))

        data = []
        rows = []
        cols = []
        for idim in xrange(dim):
            ic = nm.where(ii == idim)[0]
            if len(ic) == 0: continue

            ir = list(irs.difference([idim]))
            nn = nm.empty((len(ic), dim - 1), dtype=nm.float64)
            for ik, il in enumerate(ir):
                nn[:,ik] = - normals[ic,il] / normals[ic,idim]

            irn = dim * ic + idim
            ics = [(dim - 1) * ic + ik for ik in xrange(dim - 1)]
            for ik in xrange(dim - 1):
                rows.append(irn)
                cols.append(ics[ik])
                data.append(nn[:,ik])

            ones = nm.ones((nn.shape[0],), dtype=nm.float64)
            for ik, il in enumerate(ir):
                rows.append(dim * ic + il)
                cols.append(ics[ik])
                data.append(ones)

        rows = nm.concatenate(rows)
        cols = nm.concatenate(cols)
        data = nm.concatenate(data)

        n_np_dof = n_nod * (dim - 1)
        mtx = sp.coo_matrix((data, (rows, cols)), shape=(n_nod * dim, n_np_dof))

        self.n_dof = n_np_dof
        self.mtx = mtx.tocsr()

class NormalDirectionOperator(LCBCOperator):
    """
    Transformation matrix operator for normal direction LCBCs.

    The substitution (in 3D) is:

    .. math::
        [u_1, u_2, u_3]^T = [n_1, n_2, n_3]^T w

    The new DOF is :math:`w`.
    """
    def __init__(self, name, nodes, region, field, dof_names, filename=None):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        dim = region.dim
        assert_(len(dof_names) == dim)

        vectors = self.get_vectors(nodes, region, field, filename=filename)

        n_nod, dim = vectors.shape

        data = vectors.ravel()
        rows = nm.arange(data.shape[0])
        cols = nm.repeat(nm.arange(n_nod), dim)

        mtx = sp.coo_matrix((data, (rows, cols)), shape=(n_nod * dim, n_nod))

        self.n_dof = n_nod
        self.mtx = mtx.tocsr()

    def get_vectors(self, nodes, region, field, filename=None):
        normals = compute_nodal_normals(nodes, region, field)

        if filename is not None:
            _save_vectors(filename, normals, region, field.domain.mesh, 'n')

        return normals

class EdgeDirectionOperator(NormalDirectionOperator):
    """
    Transformation matrix operator for edges direction LCBCs.

    The substitution (in 3D) is:

    .. math::
        [u_1, u_2, u_3]^T = [d_1, d_2, d_3]^T w,

    where :math:`\ul{d}` is an edge direction vector averaged into a node. The
    new DOF is :math:`w`.
    """
    def get_vectors(self, nodes, region, field, filename=None):
        edirs = compute_nodal_edge_dirs(nodes, region, field)

        if filename is not None:
            _save_vectors(filename, edirs, region, field.domain.mesh, 'e')

        return edirs

class IntegralMeanValueOperator(LCBCOperator):
    """
    Transformation matrix operator for integral mean value LCBCs.
    All node DOFs are sumed to the new one.
    """
    def __init__(self, name, nodes, region, field, dof_names):
        Struct.__init__(self, name=name, nodes=nodes, dof_names=dof_names)

        dpn = len(dof_names)
        n_nod = nodes.shape[0]

        data = nm.ones((n_nod * dpn,))
        rows = nm.arange(data.shape[0])
        cols = nm.zeros((data.shape[0],))

        mtx = sp.coo_matrix((data, (rows, cols)), shape=(n_nod * dpn, dpn))

        self.n_dof = dpn
        self.mtx = mtx.tocsr()

class LCBCOperators(Container):
    """
    Container holding instances of LCBCOperator subclasses for a single
    variable.

    Parameters
    ----------
    name : str
        The object name.
    eq_map : EquationMap instance
        The equation mapping of the variable.
    offset : int
        The offset added to markers distinguishing the individual LCBCs.
    """
    def __init__(self, name, eq_map, offset):
        Container.__init__(self, name=name, eq_map=eq_map, offset=offset)

        self.eq_lcbc = nm.zeros((self.eq_map.n_eq,), dtype=nm.int32)
        self.markers = []
        self.n_transformed_dof = []
        self.n_op = 0
        self.ics = None

    def add_from_bc(self, bc, field):
        """
        Create a new LCBC operator described by `bc`, and add it to the
        container.

        Parameters
        ----------
        bc : LinearCombinationBC instance
            The LCBC condition description.
        field : Field instance
            The field of the variable.
        """
        region = bc.region
        dofs, kind = bc.dofs

        nmaster = field.get_dofs_in_region(region, merge=True)

        if kind == 'rigid':
            op = RigidOperator('%d_rigid' % len(self),
                               nmaster, field, dofs, self.eq_map.dof_names)

        elif kind == 'no_penetration':
            filename = bc.get('filename', None)
            op = NoPenetrationOperator('%d_no_penetration' % len(self),
                                       nmaster, region, field, dofs,
                                       filename=filename)

        elif kind == 'normal_direction':
            filename = bc.get('filename', None)
            op = NormalDirectionOperator('%d_normal_direction' % len(self),
                                         nmaster, region, field, dofs,
                                         filename=filename)

        elif kind == 'edge_direction':
            filename = bc.get('filename', None)
            op = EdgeDirectionOperator('%d_edge_direction' % len(self),
                                       nmaster, region, field, dofs,
                                       filename=filename)

        elif kind == 'integral_mean_value':
            op = IntegralMeanValueOperator('%d_integral_mean_value' % len(self),
                                           nmaster, region, field, dofs)

        self.append(op)

    def append(self, op):
        Container.append(self, op)

        eq = self.eq_map.eq
        dofs = expand_nodes_to_equations(op.nodes, op.dof_names,
                                         self.eq_map.dof_names)
        meq = eq[dofs]
        assert_(nm.all(meq >= 0))

        if self.eq_map.n_epbc:
            op.treat_pbcs(dofs, self.eq_map.master)

        self.markers.append(self.offset + self.n_op + 1)
        self.eq_lcbc[meq] = self.markers[-1]

        self.n_transformed_dof.append(op.n_dof)
        self.n_op = len(self)

    def finalize(self):
        """
        Call this after all LCBCs of the variable have been added.

        Initializes the global column indices.
        """
        self.ics = nm.cumsum(nm.r_[0, self.n_transformed_dof])

def make_global_lcbc_operator(lcbc_ops, adi, new_only=False):
    """
    Assemble all LCBC operators into a single matrix.

    Returns
    -------
    mtx_lc : csr_matrix
        The global LCBC operator in the form of a CSR matrix.
    lcdi : DofInfo
        The global active LCBC-constrained DOF information.
    new_only : bool
        If True, the operator columns will contain only new DOFs.
    """
    n_dof = adi.ptr[-1]
    eq_lcbc = nm.zeros((n_dof,), dtype=nm.int32)

    n_dof_new = 0
    n_free = {}
    n_new = {}
    for var_name, lcbc_op in lcbc_ops.iteritems():
        if lcbc_op is None: continue

        indx = adi.indx[var_name]
        eq_lcbc[indx] = lcbc_op.eq_lcbc

        n_free[var_name] = len(nm.where(lcbc_op.eq_lcbc == 0)[0])
        n_new[var_name] = nm.sum(lcbc_op.n_transformed_dof)

        n_dof_new += n_new[var_name]

    if n_dof_new == 0:
        return None, None

    ii = nm.nonzero(eq_lcbc)[0]
    n_constrained = ii.shape[0]
    n_dof_free = n_dof - n_constrained
    n_dof_reduced = n_dof_free + n_dof_new
    output('dofs: total %d, free %d, constrained %d, new %d'\
            % (n_dof, n_dof_free, n_constrained, n_dof_new))
    output(' -> reduced %d' % (n_dof_reduced))

    lcdi = DofInfo('lcbc_active_state_dof_info')
    fdi = DofInfo('free_dof_info')
    ndi = DofInfo('new_dof_info')
    for var_name in adi.var_names:
        nf = n_free.get(var_name, adi.n_dof[var_name])
        nn = n_new.get(var_name, 0)
        fdi.append_raw(var_name, nf)
        ndi.append_raw(var_name, nn)
        lcdi.append_raw(var_name, nn + nf)

    assert_(lcdi.ptr[-1] == n_dof_reduced)

    rows = []
    cols = []
    data = []
    for var_name, lcbc_op in lcbc_ops.iteritems():
        if lcbc_op is None: continue

        if new_only:
            offset = ndi.indx[var_name].start

        else:
            offset = lcdi.indx[var_name].start + fdi.n_dof[var_name]

        for ii, op in enumerate(lcbc_op):
            indx = nm.where(eq_lcbc == lcbc_op.markers[ii])[0]
            icols = nm.arange(offset + lcbc_op.ics[ii],
                              offset + lcbc_op.ics[ii+1])

            if isinstance(op.mtx, sp.spmatrix):
                lr, lc, lv = sp.find(op.mtx)
                rows.append(indx[lr])
                cols.append(icols[lc])
                data.append(lv)

            else:
                irs, ics = nm.meshgrid(indx, icols)
                rows.append(irs.ravel())
                cols.append(ics.ravel())
                data.append(op.mtx.T.ravel())

    rows = nm.concatenate(rows)
    cols = nm.concatenate(cols)
    data = nm.concatenate(data)

    if new_only:
        mtx_lc = sp.coo_matrix((data, (rows, cols)),
                               shape=(n_dof, n_dof_new))

    else:
        mtx_lc = sp.coo_matrix((data, (rows, cols)),
                               shape=(n_dof, n_dof_reduced))

        ir = nm.where(eq_lcbc == 0)[0]

        ic = nm.empty((n_dof_free,), dtype=nm.int32)
        for var_name in adi.var_names:
            ii = nm.arange(fdi.n_dof[var_name], dtype=nm.int32)
            ic[fdi.indx[var_name]] = lcdi.indx[var_name].start + ii

        mtx_lc2 = sp.coo_matrix((nm.ones((ir.shape[0],)), (ir, ic)),
                                shape=(n_dof, n_dof_reduced), dtype=nm.float64)


        mtx_lc = mtx_lc + mtx_lc2

    mtx_lc = mtx_lc.tocsr()

    return mtx_lc, lcdi
