import time
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import Struct, get_default, output, assert_
from meshio import MeshIO

def make_point_cells(indx, dim):
    conn = nm.zeros((indx.shape[0], dim + 1), dtype=nm.int32)
    for ii in range(0, dim + 1):
        conn[:,ii] = indx
    return conn

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

def make_inverse_connectivity(conns, n_nod, ret_offsets=True):
    """
    For each mesh node referenced in the connectivity conns, make a list of
    elements it belongs to.
    """
    from itertools import chain

    iconn = [[] for ii in xrange(n_nod)]
    n_els = [0] * n_nod
    for ig, conn in enumerate(conns):
        for iel, row in enumerate(conn):
            for node in row:
                iconn[node].extend([ig, iel])
                n_els[node] += 1

    n_els = nm.array(n_els, dtype=nm.int32)
    iconn = nm.fromiter(chain(*iconn), nm.int32)

    if ret_offsets:
        offsets = nm.cumsum(nm.r_[0, n_els], dtype=nm.int32)
        return offsets, iconn

    else:
        return n_els, iconn

class Mesh(Struct):
    """
    Contains the FEM mesh together with all utilities related to it.

    Input and output is handled by the MeshIO class and subclasses.
    The Mesh class only contains the real mesh - nodes, connectivity,
    regions, plus methods for doing operations on this mesh.

    Example of creating and working with a mesh::

        In [1]: from sfepy.fem import Mesh
        In [2]: m = Mesh.from_file("meshes/3d/cylinder.vtk")
        sfepy: reading mesh (meshes/3d/cylinder.vtk)...
        sfepy: ...done in 0.04 s

        In [3]: m.coors
        Out[3]:
        array([[  1.00000000e-01,   2.00000000e-02,  -1.22460635e-18],
               [  1.00000000e-01,   1.80193774e-02,   8.67767478e-03],
               [  1.00000000e-01,   1.24697960e-02,   1.56366296e-02],
               ...,
               [  8.00298527e-02,   5.21598617e-03,  -9.77772215e-05],
               [  7.02544004e-02,   3.61610291e-04,  -1.16903153e-04],
               [  3.19633596e-02,  -1.00335972e-02,   9.60460305e-03]])

        In [4]: m.ngroups
        Out[4]: array([0, 0, 0, ..., 0, 0, 0])

        In [5]: m.conns
        Out[5]:
        [array([[ 28,  60,  45,  29],
               [ 28,  60,  57,  45],
               [ 28,  57,  27,  45],
               ...,
               [353, 343, 260, 296],
               [353, 139, 181, 140],
               [353, 295, 139, 140]])]

        In [6]: m.mat_ids
        Out[6]: [array([6, 6, 6, ..., 6, 6, 6])]

        In [7]: m.descs
        Out[7]: ['3_4']

        In [8]: m
        Out[8]: Mesh:meshes/3d/cylinder

        In [9]: print m
        Mesh:meshes/3d/cylinder
          conns:
            [array([[ 28,  60,  45,  29],
                   [ 28,  60,  57,  45],
                   [ 28,  57,  27,  45],
                   ...,
                   [353, 343, 260, 296],
                   [353, 139, 181, 140],
                   [353, 295, 139, 140]])]
          coors:
            [[  1.00000000e-01   2.00000000e-02  -1.22460635e-18]
             [  1.00000000e-01   1.80193774e-02   8.67767478e-03]
             [  1.00000000e-01   1.24697960e-02   1.56366296e-02]
             ...,
             [  8.00298527e-02   5.21598617e-03  -9.77772215e-05]
             [  7.02544004e-02   3.61610291e-04  -1.16903153e-04]
             [  3.19633596e-02  -1.00335972e-02   9.60460305e-03]]
          descs:
            ['3_4']
          dim:
            3
          el_offsets:
            [   0 1348]
          io:
            None
          mat_ids:
            [array([6, 6, 6, ..., 6, 6, 6])]
          n_e_ps:
            [4]
          n_el:
            1348
          n_els:
            [1348]
          n_nod:
            354
          name:
            meshes/3d/cylinder
          ngroups:
            [0 0 0 ..., 0 0 0]
          setup_done:
            0

    The Mesh().coors is an array of node coordinates and Mesh().conns is the
    list of elements of each type (see Mesh().desc), so for example if you want
    to know the coordinates of the nodes of the fifth finite element of the
    type 3_4 do::

        In [10]: m.descs
        Out[10]: ['3_4']

    So now you know that the finite elements of the type 3_4 are in a.conns[0]::

        In [11]: m.coors[m.conns[0][4]]
        Out[11]:
        array([[  1.00000000e-01,   1.80193774e-02,  -8.67767478e-03],
               [  1.00000000e-01,   1.32888539e-02,  -4.35893200e-04],
               [  1.00000000e-01,   2.00000000e-02,  -1.22460635e-18],
               [  9.22857574e-02,   1.95180454e-02,  -4.36416134e-03]])

    The element ids are of the form "<dimension>_<number of nodes>", i.e.:

    - 2_2 ... line
    - 2_3 ... triangle
    - 2_4 ... quadrangle
    - 3_2 ... line
    - 3_4 ... tetrahedron
    - 3_8 ... hexahedron

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

        output('reading mesh (%s)...' % (io.filename))
        tt = time.clock()

        trunk = io.get_filename_trunk()
        mesh = Mesh(trunk)
        mesh = io.read(mesh, omit_facets=omit_facets)

        output('...done in %.2f s' % (time.clock() - tt))

        mesh._set_shape_info()

        return mesh

    @staticmethod
    def from_region(region, mesh_in, save_edges=False, save_faces=False,
                    localize=False, is_surface=False):
        """
        Create a mesh corresponding to a given region.
        """
        mesh = Mesh(mesh_in.name + "_reg")
        mesh.coors = mesh_in.coors.copy()
        mesh.ngroups = mesh_in.ngroups.copy()

        mesh.conns = []
        mesh.descs = []
        mesh.mat_ids = []

        if not is_surface:
            if region.has_cells():
                for ig in region.igs:
                    mesh.descs.append(mesh_in.descs[ig])
                    els = region.get_cells(ig)
                    mesh.mat_ids.append(mesh_in.mat_ids[ig][els,:].copy())
                    mesh.conns.append(mesh_in.conns[ig][els,:].copy())

            if save_edges:
                cmesh = region.domain.cmesh
                for ig in region.igs:
                    edges = region.get_edges(ig)
                    if not edges.size: continue

                    verts = cmesh.get_incident(0, edges, 1)
                    verts.shape = (verts.shape[0] / 2, 2)

                    mesh.descs.append('1_2')
                    mesh.conns.append(verts)

                    mat_ids = nm.repeat(ig, verts.shape[0])
                    mesh.mat_ids.append(mat_ids)

            if save_faces:
                mesh._append_region_faces(region)

            if save_edges or save_faces:
                mesh.descs.append('1_1')
                mesh.mat_ids.append(-nm.ones_like(region.vertices))
                mesh.conns.append(region.vertices[:, None])

        else:
            mesh._append_region_faces(region, force_faces=True)

        mesh._set_shape_info()

        if localize:
            mesh.localize(region.vertices)

        return mesh

    @staticmethod
    def from_data(name, coors, ngroups, conns, mat_ids, descs,
                  igs=None, nodal_bcs=None):
        """
        Create a mesh from mesh data.
        """
        if igs is None:
            igs = range(len(conns))
        mesh = Mesh(name)
        mesh._set_data(coors=coors,
                       ngroups=ngroups,
                       conns=[conns[ig] for ig in igs],
                       mat_ids=[mat_ids[ig] for ig in igs],
                       descs=[descs[ig] for ig in igs],
                       nodal_bcs=nodal_bcs)
        mesh._set_shape_info()
        return mesh

    def __init__(self, name='mesh', filename=None,
                 prefix_dir=None, **kwargs):
        """Create a Mesh.

        Parameters
        ----------
        name : str
            Object name.
        filename : str
            Loads a mesh from the specified file, if not None.
        prefix_dir : str
            If not None, the filename is relative to that directory.
        """
        Struct.__init__(self, name=name, **kwargs)
        self.nodal_bcs = {}

        if filename is None:
            self.io = None
            self.setup_done = 0

        else:
            io = MeshIO.any_from_filename(filename, prefix_dir=prefix_dir)
            output('reading mesh (%s)...' % (io.filename))
            tt = time.clock()
            io.read(self)
            output('...done in %.2f s' % (time.clock() - tt))
            self._set_shape_info()

    def copy(self, name=None):
        """Make a deep copy of self.

        Parameters
        ----------
        name : str
            Name of the copied mesh.
        """
        return Struct.copy(self, deep=True, name=name)

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

    def _set_shape_info(self):
        self.n_nod, self.dim = self.coors.shape
        self.n_els = nm.array([conn.shape[0] for conn in self.conns])
        self.n_e_ps = nm.array([conn.shape[1] for conn in self.conns])
        self.el_offsets = nm.cumsum(nm.r_[0, self.n_els])
        self.n_el = nm.sum(self.n_els)
        self.dims = [int(ii[0]) for ii in self.descs]

    def _set_data(self, coors, ngroups, conns, mat_ids, descs, nodal_bcs=None):
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
        self.coors = nm.ascontiguousarray(coors)

        if ngroups is None:
            self.ngroups = nm.zeros((self.coors.shape[0],), dtype=nm.int32)

        else:
            self.ngroups = nm.ascontiguousarray(ngroups)

        self.conns = [nm.asarray(conn, dtype=nm.int32) for conn in conns]
        self.mat_ids = [nm.asarray(mat_id, dtype=nm.int32)
                        for mat_id in mat_ids]
        self.descs = descs
        self.nodal_bcs = get_default(nodal_bcs, {})

    def _append_region_faces(self, region, force_faces=False):
        dim = self.coors.shape[1]
        if (not force_faces) and (dim == 2): return

        cmesh = region.domain.cmesh
        for ig in region.igs:
            faces = region.get_facets(ig)
            if not faces.size: continue

            verts, offs = cmesh.get_incident(0, faces, cmesh.dim - 1,
                                             ret_offsets=True)
            n_fp = offs[1] - offs[0]
            verts.shape = (verts.shape[0] / n_fp, n_fp)

            self.descs.append('%d_%d' % (dim - 1, n_fp))
            self.conns.append(verts)

            mat_ids = nm.repeat(ig, verts.shape[0])
            self.mat_ids.append(mat_ids)

    def write(self, filename=None, io=None,
              coors=None, igs=None, out=None, float_format=None, **kwargs):
        """
        Write mesh + optional results in `out` to a file.

        Parameters
        ----------
        filename : str, optional
            The file name. If None, the mesh name is used instead.
        io : MeshIO instance or 'auto', optional
            Passing 'auto' respects the extension of `filename`.
        coors : array, optional
            The coordinates that can be used instead of the mesh coordinates.
        igs : array_like, optional
            Passing a list of group ids selects only those groups for writing.
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

        if coors is None:
            coors = self.coors

        if igs is None:
            igs = range(len(self.conns))

        aux_mesh = Mesh.from_data(self.name, coors, self.ngroups,
                                  self.conns, self.mat_ids, self.descs,
                                  igs=igs, nodal_bcs=self.nodal_bcs)
        io.set_float_format(float_format)
        io.write(filename, aux_mesh, out, **kwargs)

    def get_bounding_box(self):
        return nm.vstack((nm.amin(self.coors, 0), nm.amax(self.coors, 0)))

    def get_element_coors(self, ig=None):
        """
        Get the coordinates of vertices elements in group `ig`.

        Parameters
        ----------
        ig : int, optional
            The element group. If None, the coordinates for all groups
            are returned, filled with zeros at places of missing
            vertices, i.e. where elements having less then the full number
            of vertices (`n_ep_max`) are.

        Returns
        -------
        coors : array
            The coordinates in an array of shape `(n_el, n_ep_max, dim)`.
        """
        cc = self.coors
        n_ep_max = self.n_e_ps.max()

        coors = nm.empty((self.n_el, n_ep_max, self.dim), dtype=cc.dtype)
        for ig, conn in enumerate(self.conns):
            i1, i2 = self.el_offsets[ig], self.el_offsets[ig + 1]
            coors[i1:i2, :conn.shape[1], :] = cc[conn]

        return coors

    def localize(self, inod):
        """
        Strips nodes not in inod and remaps connectivities.
        Omits elements where remap[conn] contains -1...
        """
        remap = nm.empty((self.n_nod,), dtype=nm.int32)
        remap.fill(-1)
        remap[inod] = nm.arange(inod.shape[0], dtype=nm.int32)

        self.coors = self.coors[inod]
        self.ngroups = self.ngroups[inod]
        conns = []
        mat_ids = []
        for ig, conn in enumerate(self.conns):
            if conn.shape[0] == 0:
                continue

            aux = remap[conn]
            ii = nm.unique(nm.where(aux == -1)[0])
            ii = nm.setdiff1d(nm.arange(conn.shape[0], dtype=nm.int32), ii)
            conns.append(aux[ii])
            mat_ids.append(self.mat_ids[ig][ii])
        self.conns = conns
        self.mat_ids = mat_ids

        self._set_shape_info()

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
