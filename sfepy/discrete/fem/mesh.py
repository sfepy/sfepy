import time
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import Struct, invert_dict, get_default, output, assert_
from meshio import MeshIO, supported_cell_types

def find_map(x1, x2, eps=1e-8, allow_double=False, join=True):
    """
    Find a mapping between common coordinates in x1 and x2, such that
    x1[cmap[:,0]] == x2[cmap[:,1]]
    """
    off, dim = x1.shape
    ir = nm.zeros((off + x2.shape[0],), dtype=nm.int32)
    ir[off:] = off

    x1 = nm.round(x1.T / eps) * eps
    x2 = nm.round(x2.T / eps) * eps
    xx = nm.c_[x1, x2]

    keys = [xx[ii] for ii in range(dim)]
    iis = nm.lexsort(keys=keys)

    xs = xx.T[iis]
    xd = nm.sqrt(nm.sum(nm.diff(xs, axis=0)**2.0, axis=1))

    ii = nm.where(xd < eps)[0]
    off1, off2 = ir[iis][ii], ir[iis][ii+1]
    i1, i2 = iis[ii] - off1, iis[ii+1] - off2
    dns = nm.where(off1 == off2)[0]
    if dns.size:
        output('double node(s) in:')
        for dn in dns:
            if off1[dn] == 0:
                output('x1: %d %d -> %s %s' % (i1[dn], i2[dn],
                                               x1[:,i1[dn]], x1[:,i2[dn]]))
            else:
                output('x2: %d %d -> %s %s' % (i1[dn], i2[dn],
                                               x2[:,i1[dn]], x2[:,i2[dn]]))
        if not allow_double:
            raise ValueError('double node(s)! (see above)')

    if join:
        cmap = nm.c_[i1, i2]
        return cmap
    else:
        return i1, i2

def merge_mesh(x1, ngroups1, conns1, mat_ids1, x2, ngroups2, conns2, mat_ids2,
               cmap, eps=1e-8):
    """
    Merge two meshes in common coordinates found in x1, x2.

    Notes
    -----
    Assumes the same number and kind of element groups in both meshes!
    """
    n1 = x1.shape[0]
    n2 = x2.shape[0]

    err = nm.sum(nm.sum(nm.abs(x1[cmap[:,0],:-1] - x2[cmap[:,1],:-1])))
    if abs(err) > (10.0 * eps):
        raise ValueError('nonmatching meshes! (error: %e)' % err)

    mask = nm.ones((n2,), dtype=nm.int32)
    mask[cmap[:,1]] = 0
    remap = nm.cumsum(mask) + n1 - 1
    remap[cmap[:,1]] = cmap[:,0]

    i2 = nm.setdiff1d(nm.arange( n2, dtype=nm.int32), cmap[:,1])
    xx = nm.r_[x1, x2[i2]]
    ngroups = nm.r_[ngroups1, ngroups2[i2]]

    conns = []
    for ii, conn1 in enumerate(conns1):
        conn = nm.vstack((conn1, remap[conns2[ii]]))
        conns.append(conn)

    mat_ids = None
    if (mat_ids1 is not None) and (mat_ids2 is not None):
        mat_ids = []

        for ii, mm1 in enumerate(mat_ids1):
            mm = nm.concatenate((mm1, mat_ids2[ii]))
            mat_ids.append(mm)

    return xx, ngroups, conns, mat_ids

def fix_double_nodes(coor, ngroups, conns, eps):
    """
    Detect and attempt fixing double nodes in a mesh.

    The double nodes are nodes having the same coordinates
    w.r.t. precision given by `eps`.
    """
    n_nod, dim = coor.shape
    cmap = find_map(coor, nm.zeros((0,dim)), eps=eps, allow_double=True)
    if cmap.size:
        output('double nodes in input mesh!')
        output('trying to fix...')

        while cmap.size:
            # Just like in Variable.equation_mapping()...
            ii = nm.argsort(cmap[:,1])
            scmap = cmap[ii]

            eq = nm.arange(n_nod)
            eq[scmap[:,1]] = -1
            eqi = eq[eq >= 0]
            eq[eqi] = nm.arange(eqi.shape[0])
            remap = eq.copy()
            remap[scmap[:,1]] = eq[scmap[:,0]]
            output(coor.shape)
            coor = coor[eqi]
            ngroups = ngroups[eqi]
            output(coor.shape)
            ccs = []
            for conn in conns:
                ccs.append(remap[conn])
            conns = ccs
            cmap = find_map(coor, nm.zeros((0,dim)), eps=eps,
                            allow_double=True)
        output('...done')
    return coor, ngroups, conns

def get_min_edge_size(coor, conns):
    """
    Get the smallest edge length.
    """
    mes = 1e16
    for conn in conns:
        n_ep = conn.shape[1]
        for ir in range(n_ep):
            x1 = coor[conn[:,ir]]
            for ic in range(ir + 1, n_ep):
                x2 = coor[conn[:,ic]]
                aux = nm.sqrt(nm.sum((x2 - x1)**2.0, axis=1).min())
                mes = min(mes, aux)

    return mes

def get_min_vertex_distance(coor, guess):
    """Can miss the minimum, but is enough for our purposes."""
    # Sort by x.
    ix = nm.argsort(coor[:,0])
    scoor = coor[ix]

    mvd = 1e16

    # Get mvd in chunks potentially smaller than guess.
    n_coor = coor.shape[0]

    i0 = i1 = 0
    x0 = scoor[i0,0]
    while 1:
        while ((scoor[i1,0] - x0) < guess) and (i1 < (n_coor - 1)):
            i1 += 1

        ## print i0, i1, x0, scoor[i1,0]
        aim, aa1, aa2, aux = get_min_vertex_distance_naive(scoor[i0:i1+1])
        if aux < mvd:
            im, a1, a2 = aim, aa1 + i0, aa2 + i0
        mvd = min(mvd, aux)
        i0 = i1 = int(0.5 * (i1 + i0)) + 1
        ## i0 += 1
        x0 = scoor[i0,0]
        ## print '-', i0

        if i1 == n_coor - 1: break

    ## print im, ix[a1], ix[a2], a1, a2, scoor[a1], scoor[a2]

    return mvd

def get_min_vertex_distance_naive(coor):

    ii = nm.arange(coor.shape[0])
    i1, i2 = nm.meshgrid(ii, ii)
    i1 = i1.flatten()
    i2 = i2.flatten()

    ii = nm.where(i1 < i2)
    aux = coor[i1[ii]] - coor[i2[ii]]
    aux = nm.sum(aux**2.0, axis=1)

    im = aux.argmin()

    return im, i1[ii][im], i2[ii][im], nm.sqrt(aux[im])

def make_mesh(coor, ngroups, conns, mesh_in):
    """Create a mesh reusing mat_ids and descs of mesh_in."""
    mat_ids = []
    for ii, conn in enumerate(conns):
        mat_id = nm.empty((conn.shape[0],), dtype=nm.int32)
        mat_id.fill(mesh_in.mat_ids[ii][0])
        mat_ids.append(mat_id)

    mesh_out = Mesh.from_data('merged mesh', coor, ngroups, conns,
                              mat_ids, mesh_in.descs)
    return mesh_out

class Mesh(Struct):
    """
    The Mesh class is a light proxy to CMesh.

    Input and output is handled by the MeshIO class and subclasses.
    """

    @staticmethod
    def from_surface(surf_faces, mesh_in):
        """
        Create a mesh given a set of surface faces and the original mesh.
        """
        aux = nm.concatenate([faces.ravel() for faces in surf_faces])
        inod = nm.unique(aux)

        n_nod = len(inod)
        n_nod_m, dim = mesh_in.coors.shape

        aux = nm.arange(n_nod, dtype=nm.int32)
        remap = nm.zeros((n_nod_m,), nm.int32)
        remap[inod] = aux

        mesh = Mesh(mesh_in.name + "_surf")

        mesh.coors = mesh_in.coors[inod]
        mesh.ngroups = mesh_in.ngroups[inod]

        sfm = {3 : "2_3", 4 : "2_4"}
        mesh.conns = []
        mesh.descs = []
        mesh.mat_ids = []
        for ii, sf in enumerate(surf_faces):
            n_el, n_fp = sf.shape

            conn = remap[sf]
            mat_id = nm.empty((conn.shape[0],), dtype=nm.int32)
            mat_id.fill(ii)

            mesh.descs.append(sfm[n_fp])
            mesh.conns.append(conn)
            mesh.mat_ids.append(mat_id)

        mesh._set_shape_info()

        return mesh

    @staticmethod
    def from_file(filename=None, io='auto', prefix_dir=None,
                  omit_facets=False):
        """
        Read a mesh from a file.

        Parameters
        ----------
        filename : string or function or MeshIO instance or Mesh instance
            The name of file to read the mesh from. For convenience, a
            mesh creation function or a MeshIO instance or directly a Mesh
            instance can be passed in place of the file name.
        io : *MeshIO instance
            Passing *MeshIO instance has precedence over filename.
        prefix_dir : str
            If not None, the filename is relative to that directory.
        omit_facets : bool
            If True, do not read cells of lower dimension than the space
            dimension (faces and/or edges). Only some MeshIO subclasses
            support this!
        """
        if isinstance(filename, Mesh):
            return filename

        if io == 'auto':
            if filename is None:
                output('filename or io must be specified!')
                raise ValueError
            else:
                io = MeshIO.any_from_filename(filename, prefix_dir=prefix_dir)

        cell_types = ', '.join(supported_cell_types[io.format])
        output('reading mesh [%s] (%s)...' % (cell_types, io.filename))
        tt = time.clock()

        trunk = io.get_filename_trunk()
        mesh = Mesh(trunk)
        mesh = io.read(mesh, omit_facets=omit_facets)

        output('...done in %.2f s' % (time.clock() - tt))

        mesh._set_shape_info()

        return mesh

    @staticmethod
    def from_region(region, mesh_in, localize=False, is_surface=False):
        """
        Create a mesh corresponding to cells, or, if `is_surface` is True, to
        facets, of a given region.
        """
        cmesh_in = mesh_in.cmesh

        if not is_surface:
            if not region.shape.n_cell:
                raise ValueError('region %s has no cells!' % region.name)

            cmesh = cmesh_in.create_new(region.cells, region.tdim,
                                        localize=localize)

        else:
            if not region.shape.n_facet:
                raise ValueError('region %s has no facets!' % region.name)

            cmesh = cmesh_in.create_new(region.facets, region.tdim - 1,
                                        localize=localize)

        mesh = Mesh(mesh_in.name + "_reg", cmesh=cmesh)

        return mesh

    @staticmethod
    def from_data(name, coors, ngroups, conns, mat_ids, descs,
                  nodal_bcs=None):
        """
        Create a mesh from mesh IO data.
        """
        mesh = Mesh(name)
        mesh._set_io_data(coors=coors,
                          ngroups=ngroups,
                          conns=conns,
                          mat_ids=mat_ids,
                          descs=descs,
                          nodal_bcs=nodal_bcs)
        mesh._set_shape_info()
        return mesh

    def __init__(self, name='mesh', cmesh=None):
        """
        Create a Mesh.

        By default, the mesh is empty, see the `cmesh` argument.

        Parameters
        ----------
        name : str
            Object name.
        cmesh : CMesh, optional
            If given, use this as the cmesh.
        """
        Struct.__init__(self, name=name, nodal_bcs={}, io=None)
        if cmesh is not None:
            self.cmesh = cmesh
            self._collect_descs()
            self._set_shape_info()

    def copy(self, name=None):
        """
        Make a deep copy of the mesh.

        Parameters
        ----------
        name : str
            Name of the copied mesh.
        """
        cmesh = self.cmesh.create_new()
        return Mesh(name=name, cmesh=cmesh)

    def __add__(self, other):
        """
        Merge the two meshes, assuming they have the same number and kind of
        element groups.
        """
        cmap = find_map(self.coors, other.coors)
        aux = merge_mesh(self.coors, self.ngroups, self.conns, self.mat_ids,
                         other.coors, other.ngroups, other.conns, other.mat_ids,
                         cmap)
        coors, ngroups, conns, mat_ids = aux

        mesh = Mesh.from_data(self.name + ' + ' + other.name,
                              coors, ngroups, conns, mat_ids, self.descs)

        return mesh

    def _collect_descs(self):
        cmesh = self.cmesh
        i2k = invert_dict(cmesh.key_to_index)
        cts = nm.unique(cmesh.cell_types)
        self.descs = [i2k[ct] for ct in cts]

    def _set_shape_info(self):
        self.n_nod, self.dim = self.coors.shape
        self.n_el = self.cmesh.n_el
        self.dims = [int(ii[0]) for ii in self.descs]

    def _set_io_data(self, coors, ngroups, conns, mat_ids, descs,
                     nodal_bcs=None):
        """
        Set mesh data.

        Parameters
        ----------
        coors : array
            Coordinates of mesh nodes.
        ngroups : array
            Node groups.
        conns : list of arrays
            The array of mesh elements (connectivities) for each element group.
        mat_ids : list of arrays
            The array of material ids for each element group.
        descs: list of strings
            The element type for each element group.
        nodal_bcs : dict of arrays, optional
            The nodes defining regions for boundary conditions referred
            to by the dict keys in problem description files.
        """
        ac = nm.ascontiguousarray
        coors = ac(coors, dtype=nm.float64)

        if ngroups is None:
            ngroups = nm.zeros((coors.shape[0],), dtype=nm.int32)

        self.descs = descs
        self.nodal_bcs = get_default(nodal_bcs, {})

        from sfepy.discrete.fem.extmods.cmesh import CMesh
        self.cmesh = CMesh.from_data(coors, ac(ngroups),
                                     [ac(conn, dtype=nm.int32)
                                      for conn in conns],
                                     ac(nm.concatenate(mat_ids)), descs)

    def _get_io_data(self):
        """
        Return data to be used by `MeshIO`.
        """
        cmesh = self.cmesh
        conns, mat_ids = [], []
        for desc in self.descs:
            conn, cells = self.get_conn(desc, ret_cells=True)
            conns.append(conn)
            mat_ids.append(cmesh.cell_groups[cells])

        return cmesh.coors, cmesh.vertex_groups, conns, mat_ids, self.descs

    @property
    def coors(self):
        return self.cmesh.coors

    def write(self, filename=None, io=None, out=None, float_format=None,
              **kwargs):
        """
        Write mesh + optional results in `out` to a file.

        Parameters
        ----------
        filename : str, optional
            The file name. If None, the mesh name is used instead.
        io : MeshIO instance or 'auto', optional
            Passing 'auto' respects the extension of `filename`.
        out : dict, optional
            The output data attached to the mesh vertices and/or cells.
        float_format : str, optional
            The format string used to print floats in case of a text file
            format.
        **kwargs : dict, optional
            Additional arguments that can be passed to the `MeshIO` instance.
        """
        if filename is None:
            filename = self.name + '.mesh'

        if io is None:
            io = self.io
            if io is None:
                io = 'auto'

        if io == 'auto':
            io = MeshIO.any_from_filename(filename)

        io.set_float_format(float_format)
        io.write(filename, self, out, **kwargs)

    def get_bounding_box(self):
        return nm.vstack((nm.amin(self.coors, 0), nm.amax(self.coors, 0)))

    def get_conn(self, desc, ret_cells=False):
        """
        Get the rectangular cell-vertex connectivity corresponding to `desc`.
        If `ret_cells` is True, the corresponding cells are returned as well.
        """
        if desc not in self.descs:
            raise ValueError("'%s' not in %s!" % (desc, self.descs))

        cmesh = self.cmesh

        cells = nm.where(cmesh.cell_types == cmesh.key_to_index[desc])
        cells = cells[0].astype(nm.uint32)

        cdim = int(desc[0])
        conn = cmesh.get_incident(0, cells, cdim)
        conn = conn.reshape((cells.shape[0], -1)).astype(nm.int32)

        if ret_cells:
            return conn, cells

        else:
            return conn

    def transform_coors(self, mtx_t, ref_coors=None):
        """
        Transform coordinates of the mesh by the given transformation matrix.

        Parameters
        ----------
        mtx_t : array
           The transformation matrix `T` (2D array). It is applied
           depending on its shape:

           - `(dim, dim): x = T * x`
           - `(dim, dim + 1): x = T[:, :-1] * x + T[:, -1]`
        ref_coors : array, optional
           Alternative coordinates to use for the transformation instead
           of the mesh coordinates, with the same shape as `self.coors`.
        """
        if ref_coors is None:
            ref_coors = self.coors

        if mtx_t.shape[1] > self.coors.shape[1]:
            self.coors[:] = nm.dot(ref_coors, mtx_t[:,:-1].T) + mtx_t[:,-1]
        else:
            self.coors[:] = nm.dot(ref_coors, mtx_t.T)

    def create_conn_graph(self, verbose=True):
        """
        Create a graph of mesh connectivity.

        Returns
        -------
        graph : csr_matrix
            The mesh connectivity graph as a SciPy CSR matrix.
        """
        from extmods.cmesh import create_mesh_graph

        shape = (self.n_nod, self.n_nod)
        output('graph shape:', shape, verbose=verbose)
        if nm.prod(shape) == 0:
            output('no graph (zero size)!', verbose=verbose)
            return None

        output('assembling mesh graph...', verbose=verbose)
        tt = time.clock()

        nnz, prow, icol = create_mesh_graph(shape[0], shape[1],
                                            len(self.conns),
                                            self.conns, self.conns)
        output('...done in %.2f s' % (time.clock() - tt), verbose=verbose)
        output('graph nonzeros: %d (%.2e%% fill)' \
               % (nnz, float(nnz) / nm.prod(shape)))

        data = nm.ones((nnz,), dtype=nm.bool)
        graph = sp.csr_matrix((data, icol, prow), shape)

        return graph

    def explode_groups(self, eps, return_emap=False):
        """
        Explode the mesh element groups by `eps`, i.e. split group
        interface nodes and shrink each group towards its centre by
        `eps`.

        Parameters
        ----------
        eps : float in `[0.0, 1.0]`
            The group shrinking factor.
        return_emap : bool, optional
            If True, also return the mapping against original mesh
            coordinates that result in the exploded mesh coordinates.
            The mapping can be used to map mesh vertex data to the
            exploded mesh vertices.

        Returns
        -------
        mesh : Mesh
            The new mesh with exploded groups.
        emap : spmatrix, optional
            The maping for exploding vertex values. Only provided if
            `return_emap` is True.
        """
        assert_(0.0 <= eps <= 1.0)

        remap = nm.empty((self.n_nod,), dtype=nm.int32)
        offset = 0

        if return_emap:
            rows, cols = [], []

        coors = []
        ngroups = []
        conns = []
        mat_ids = []
        descs = []
        for ig, conn in enumerate(self.conns):
            nodes = nm.unique(conn)
            group_coors = self.coors[nodes]
            n_nod = group_coors.shape[0]

            centre = group_coors.sum(axis=0) / float(n_nod)
            vectors = group_coors - centre[None, :]
            new_coors = centre + (vectors * eps)

            remap[nodes] = nm.arange(n_nod, dtype=nm.int32) + offset
            new_conn = remap[conn]

            coors.append(new_coors)
            ngroups.append(self.ngroups[nodes])

            conns.append(new_conn)
            mat_ids.append(self.mat_ids[ig])
            descs.append(self.descs[ig])

            offset += n_nod

            if return_emap:
                cols.append(nodes)
                rows.append(remap[nodes])

        coors = nm.concatenate(coors, axis=0)
        ngroups = nm.concatenate(ngroups, axis=0)

        mesh = Mesh.from_data('exploded_' + self.name,
                              coors, ngroups, conns, mat_ids, descs)

        if return_emap:
            rows = nm.concatenate(rows)
            cols = nm.concatenate(cols)
            data = nm.ones(rows.shape[0], dtype=nm.float64)

            emap = sp.coo_matrix((data, (rows, cols)),
                                 shape=(mesh.n_nod, self.n_nod))

            return mesh, emap

        else:
            return mesh
