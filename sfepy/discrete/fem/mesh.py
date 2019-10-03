from __future__ import absolute_import
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import Struct, invert_dict, get_default, output,\
     assert_, is_sequence
from sfepy.base.timing import Timer
from .meshio import MeshIO
import six
from scipy.spatial import cKDTree

eps = 1e-9


def set_accuracy(eps):
    globals()['eps'] = eps

def find_map(x1, x2, allow_double=False, join=True):
    """
    Find a mapping between common coordinates in x1 and x2, such that
    x1[cmap[:,0]] == x2[cmap[:,1]]
    """
    kdtree = cKDTree(nm.vstack([x1, x2]))
    cmap = kdtree.query_pairs(eps, output_type='ndarray')

    dns1 = nm.where(cmap[:, 1] < x1.shape[0])[0]
    dns2 = nm.where(cmap[:, 0] >= x1.shape[0])[0]

    if (dns1.size + dns2.size):
        output('double node(s) in:')
        for dn in dns1:
            idxs = cmap[dn, :]
            output('x1: %d %d -> %s %s' % (idxs[0], idxs[1],
                                           x1[idxs[0], :], x1[idxs[1], :]))
        for dn in dns2:
            idxs = cmap[dn, :]
            output('x2: %d %d -> %s %s' % (idxs[0], idxs[1],
                                           x2[idxs[0], :], x2[idxs[1], :]))

        if not allow_double:
            raise ValueError('double node(s)! (see above)')

    cmap[:, 1] -= x1.shape[0]

    return cmap if join else (cmap[:, 0], cmap[:, 1])

def merge_mesh(x1, ngroups1, conn1, mat_ids1, x2, ngroups2, conn2, mat_ids2,
               cmap):
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

    conn = nm.vstack((conn1, remap[conn2]))

    mat_ids = None
    if (mat_ids1 is not None) and (mat_ids2 is not None):
        mat_ids = nm.concatenate((mat_ids1, mat_ids2))

    return xx, ngroups, conn, mat_ids

def fix_double_nodes(coor, ngroups, conns):
    """
    Detect and attempt fixing double nodes in a mesh.

    The double nodes are nodes having the same coordinates
    w.r.t. precision given by `eps`.
    """
    n_nod, dim = coor.shape
    cmap = find_map(coor, nm.zeros((0,dim)), allow_double=True)
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
            cmap = find_map(coor, nm.zeros((0,dim)), allow_double=True)
        output('...done')
    return coor, ngroups, conns

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
    def from_file(filename=None, io='auto', prefix_dir=None,
                  omit_facets=False, file_format=None):
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
                io = MeshIO.any_from_filename(filename, prefix_dir=prefix_dir,
                                              file_format=file_format)

        output('reading mesh (%s)...' % io.filename)
        timer = Timer(start=True)

        trunk = io.get_filename_trunk()
        mesh = Mesh(trunk)
        mesh = io.read(mesh, omit_facets=omit_facets)

        # FIXME - hot fix for reading 1D meshes
        if len(mesh.descs) == 1 and mesh.descs[0] == "1_2":
            output("forcing 1D")
            data = list(mesh._get_io_data(cell_dim_only=1))
            data[0] = data[0][:, :1]
            mesh = Mesh.from_data(mesh.name, *data)

        output('...done in %.2f s' % timer.stop())

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

        if localize:
            remap = nm.empty(mesh_in.cmesh.n_coor, dtype=nm.int32)
            remap[:] = -1
            remap[region.vertices] = nm.arange(region.vertices.shape[0])

            mesh.nodal_bcs = {}
            for key, val in six.iteritems(mesh_in.nodal_bcs):
                new_val = remap[val]
                mesh.nodal_bcs[key] = new_val[new_val >= 0]

        else:
            mesh.nodal_bcs = mesh_in.nodal_bcs.copy()

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
        if name is None:
            name = self.name

        cmesh = self.cmesh.create_new()
        return Mesh(name=name, cmesh=cmesh)

    def __add__(self, other):
        """
        Merge the two meshes, assuming they have the same kind of the single
        element group.
        """
        cmap = find_map(self.coors, other.coors)
        desc = self.descs[0]
        aux = merge_mesh(self.coors, self.cmesh.vertex_groups,
                         self.get_conn(desc), self.cmesh.cell_groups,
                         other.coors, other.cmesh.vertex_groups,
                         other.get_conn(desc), other.cmesh.cell_groups,
                         cmap)
        coors, ngroups, conn, mat_ids = aux

        mesh = Mesh.from_data(self.name + ' + ' + other.name,
                              coors, ngroups, [conn], [mat_ids], [desc])

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

        from sfepy.discrete.common.extmods.cmesh import CMesh
        self.cmesh = CMesh.from_data(coors, ac(ngroups),
                                     [ac(conn, dtype=nm.int32)
                                      for conn in conns],
                                     ac(nm.concatenate(mat_ids)), descs)

    def _get_io_data(self, cell_dim_only=None):
        """
        Return data to be used by `MeshIO`.
        """
        cmesh = self.cmesh
        conns, mat_ids = [], []
        if cell_dim_only is not None:
            if not is_sequence(cell_dim_only):
                cell_dim_only = [cell_dim_only]
            descs = [ii for ii in self.descs if int(ii[0]) in cell_dim_only]
        else:
            descs = self.descs
        for desc in descs:
            conn, cells = self.get_conn(desc, ret_cells=True)
            conns.append(conn)
            mat_ids.append(cmesh.cell_groups[cells])

        return cmesh.coors, cmesh.vertex_groups, conns, mat_ids, descs

    @property
    def coors(self):
        return self.cmesh.coors

    def write(self, filename=None, io=None, out=None, float_format=None,
              file_format=None, **kwargs):
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
            io = MeshIO.any_from_filename(filename, file_format=file_format,
                                          mode='w')

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
        from sfepy.discrete.common.extmods.cmesh import create_mesh_graph

        shape = (self.n_nod, self.n_nod)
        output('graph shape:', shape, verbose=verbose)
        if nm.prod(shape) == 0:
            output('no graph (zero size)!', verbose=verbose)
            return None

        output('assembling mesh graph...', verbose=verbose)
        timer = Timer(start=True)

        conn = self.get_conn(self.descs[0])
        nnz, prow, icol = create_mesh_graph(shape[0], shape[1],
                                            1, [conn], [conn])
        output('...done in %.2f s' % timer.stop(), verbose=verbose)
        output('graph nonzeros: %d (%.2e%% fill)' \
               % (nnz, float(nnz) / nm.prod(shape)), verbose=verbose)

        data = nm.ones((nnz,), dtype=nm.bool)
        graph = sp.csr_matrix((data, icol, prow), shape)

        return graph
