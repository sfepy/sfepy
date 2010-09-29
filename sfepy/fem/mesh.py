import time
import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import Struct, output, assert_
from meshio import MeshIO

##
# 28.05.2007, c
def make_point_cells( indx, dim ):
    conn = nm.zeros( (indx.shape[0], dim + 1), dtype = nm.int32 )
    for ii in range( 0, dim + 1 ):
        conn[:,ii] = indx
    return conn

##
# 23.05.2007, updated from matlab version, r: 05.05.2008
def find_map( x1, x2, eps = 1e-8, allow_double = False, join = True ):
    """
    Find a mapping between common coordinates in x1 and x2, such that
    x1[cmap[:,0]] == x2[cmap[:,1]]
    """
    off, dim = x1.shape
    ir = nm.zeros( (off + x2.shape[0],), dtype = nm.int32 )
    ir[off:] = off

    x1 = nm.round( x1.T / eps ) * eps
    x2 = nm.round( x2.T / eps ) * eps
    xx = nm.c_[x1, x2]

    keys = [xx[ii] for ii in range( dim )]
    iis = nm.lexsort( keys = keys )

    xs = xx.T[iis]
##     import scipy.io as io
##     io.write_array( 'sss1', x1.T )
##     io.write_array( 'sss2', x2.T )
##     io.write_array( 'sss', xs, precision = 16 )
##     pause()
    xd = nm.sqrt( nm.sum( nm.diff( xs, axis = 0 )**2.0, axis = 1 ) )

    ii = nm.where( xd < eps )[0]
    off1, off2 = ir[iis][ii], ir[iis][ii+1]
    i1, i2 = iis[ii] - off1, iis[ii+1] - off2
    dns = nm.where( off1 == off2 )[0]
    if dns.size:
        print 'double node(s) in:'
        for dn in dns:
            if off1[dn] == 0:
                print 'x1: %d %d -> %s %s' % (i1[dn], i2[dn],
                                              x1[:,i1[dn]], x1[:,i2[dn]])
            else:
                print 'x2: %d %d -> %s %s' % (i1[dn], i2[dn],
                                              x2[:,i1[dn]], x2[:,i2[dn]])
        if not allow_double:
            raise ValueError

    if join:
        cmap = nm.c_[i1, i2]
        return cmap
    else:
        return i1, i2

def merge_mesh( x1, ngroups1, conns1, x2, ngroups2, conns2, cmap, eps = 1e-8 ):
    """Merge two meshes in common coordinates found in x1, x2."""
    n1 = x1.shape[0]
    n2 = x2.shape[0]

    err = nm.sum( nm.sum( nm.abs( x1[cmap[:,0],:-1] - x2[cmap[:,1],:-1] ) ) )
    if abs( err ) > (10.0 * eps):
        print 'nonmatching meshes!', err
        raise ValueError

    mask = nm.ones( (n2,), dtype = nm.int32 )
    mask[cmap[:,1]] = 0
#    print mask, nm.cumsum( mask )
    remap = nm.cumsum( mask ) + n1 - 1
    remap[cmap[:,1]] = cmap[:,0]
#    print remap

    i2 = nm.setdiff1d( nm.arange(  n2, dtype = nm.int32 ),
                       cmap[:,1] )
    xx = nm.r_[x1, x2[i2]]
    ngroups = nm.r_[ngroups1, ngroups2[i2]]

    conns = []
    for ii in xrange( len( conns1 ) ):
        conn = nm.vstack( (conns1[ii], remap[conns2[ii]]) )
        conns.append( conn )
    
    return xx, ngroups, conns

def make_mesh( coor, ngroups, conns, mesh_in ):
    """Create a mesh reusing mat_ids and descs of mesh_in."""
    mat_ids = []
    for ii, conn in enumerate( conns ):
        mat_id = nm.empty( (conn.shape[0],), dtype = nm.int32 )
        mat_id.fill( mesh_in.mat_ids[ii][0] )
        mat_ids.append( mat_id )
        
    mesh_out = Mesh.from_data( 'merged mesh', coor, ngroups, conns,
                               mat_ids, mesh_in.descs )
    return mesh_out

def make_inverse_connectivity(conns, n_nod, ret_offsets=True):
    """
    For each mesh node referenced in the connectivity conns, make a list of
    elements it belongs to.
    """
    from itertools import chain

    iconn = [[] for ii in xrange( n_nod )]
    n_els = [0] * n_nod
    for ig, conn in enumerate( conns ):
        for iel, row in enumerate( conn ):
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

##
# Mesh.
# 13.12.2004, c
# 02.01.2005
class Mesh( Struct ):
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

    def from_surface( surf_faces, mesh_in ):
        """
        Create a mesh given a set of surface faces and the original mesh.
        """
        aux = nm.concatenate([faces.ravel() for faces in surf_faces])
        inod = nm.unique1d(aux)

        n_nod = len( inod )
        n_nod_m, dim = mesh_in.coors.shape

        aux = nm.arange( n_nod, dtype=nm.int32 )
        remap = nm.zeros( (n_nod_m,), nm.int32 )
        remap[inod] = aux

        mesh = Mesh( mesh_in.name + "_surf" )

        mesh.coors = mesh_in.coors[inod]
        mesh.ngroups = mesh_in.ngroups[inod]

        sfm = {3 : "2_3", 4 : "2_4"}
        mesh.conns = []
        mesh.descs = []
        mesh.mat_ids = []
        for ii, sf in enumerate( surf_faces ):
            n_el, n_fp = sf.shape

            conn = remap[sf]
            mat_id = nm.empty( (conn.shape[0],), dtype = nm.int32 )
            mat_id.fill( ii )

            mesh.descs.append( sfm[n_fp] )
            mesh.conns.append( conn )
            mesh.mat_ids.append( mat_id )

        mesh._set_shape_info()
        
        return mesh
    from_surface = staticmethod( from_surface )

    @staticmethod
    def from_file(filename=None, io='auto', prefix_dir=None):
        """
        Read a mesh from a file.

        Parameters
        ----------
        filename : string like
            The filename.

        io : *MeshIO instance
            Passing *MeshIO instance has precedence over filename.

        prefix_dir : str
            If not None, the filename is relative to that directory.
        """
        if io == 'auto':
            if filename is None:
                output( 'filename or io must be specified!' )
                raise ValueError
            else:
                io = MeshIO.any_from_filename(filename, prefix_dir=prefix_dir)

        output('reading mesh (%s)...' % (io.filename))
        tt = time.clock()

        trunk = io.get_filename_trunk()
        mesh = Mesh(trunk)
        mesh = io.read(mesh)

        output('...done in %.2f s' % (time.clock() - tt))

        mesh._set_shape_info()

        return mesh

    @staticmethod
    def from_region(region, mesh_in, save_edges=False, save_faces=False,
                    localize=False, is_surface=False):
        """
        Create a mesh corresponding to a given region.
        """
        mesh = Mesh( mesh_in.name + "_reg" )
        mesh.coors = mesh_in.coors.copy()
        mesh.ngroups = mesh_in.ngroups.copy()

        mesh.conns = []
        mesh.descs = []
        mesh.mat_ids = []

        if not is_surface:

            if region.has_cells():
                for ig in region.igs:
                    mesh.descs.append( mesh_in.descs[ig] )
                    els = region.get_cells( ig )
                    mesh.mat_ids.append( mesh_in.mat_ids[ig][els,:].copy() )
                    mesh.conns.append( mesh_in.conns[ig][els,:].copy() )

            if save_edges:
                ed = region.domain.ed
                for ig in region.igs:
                    edges = region.get_edges( ig )
                    mesh.descs.append( '1_2' )
                    mesh.mat_ids.append( ed.data[edges,0] + 1 )
                    mesh.conns.append( ed.data[edges,-2:].copy() )

            if save_faces:
                mesh._append_region_faces(region)

            if save_edges or save_faces:
                mesh.descs.append( {2 : '2_3', 3 : '3_4'}[mesh_in.dim] )
                mesh.mat_ids.append( -nm.ones_like( region.all_vertices ) )
                mesh.conns.append(make_point_cells(region.all_vertices,
                                                   mesh_in.dim))

        else:
            mesh._append_region_faces(region)

        mesh._set_shape_info()

        if localize:
            mesh.localize( region.all_vertices )

        return mesh

    ##
    # c: 02.01.2008, r: 02.01.2008
    def from_region_and_field( region, field ):
        mesh = Mesh.from_region(region, field.domain.mesh)
        mesh.name = mesh.name + '_field'

        nodes = region.get_field_nodes( field, merge = True )

        aux = field.get_extra_nodes_as_simplices( nodes )
        mesh.coors = field.aps.coors
        mesh.ngroups = nm.zeros( (mesh.coors.shape[0],), dtype = nm.int32 )
        mesh.descs.append( aux[0] )
        mesh.mat_ids.append( aux[1] )
        mesh.conns.append( aux[2] )

        mesh.localize( nodes )

        return mesh
    from_region_and_field = staticmethod( from_region_and_field )

    def from_data( name, coors, ngroups, conns, mat_ids, descs, igs = None ):
        """
        Create a mesh from mesh data.
        """
        if igs is None:
            igs = range( len( conns ) )
        mesh = Mesh(name)
        mesh._set_data(coors = coors,
                       ngroups = ngroups,
                       conns = [conns[ig] for ig in igs],
                       mat_ids = [mat_ids[ig] for ig in igs],
                       descs = [descs[ig] for ig in igs])
        mesh._set_shape_info()
        return mesh
    from_data = staticmethod( from_data )
        

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

        if filename is None:
            self.io = None
            self.setup_done = 0

        else:
            io = MeshIO.any_from_filename(filename, prefix_dir=prefix_dir)
            output( 'reading mesh (%s)...' % (io.filename) )
            tt = time.clock()
            io.read(self)
            output( '...done in %.2f s' % (time.clock() - tt) )
            self._set_shape_info()
            
    def copy(self, name=None):
        """Make a deep copy of self.

        Parameters
        ----------
        name : str
            Name of the copied mesh.
        """
        return Struct.copy(self, deep=True, name=name)

    ##
    # 04.08.2006, c
    # 29.09.2006
    def _set_shape_info( self ):
        self.n_nod, self.dim = self.coors.shape
        self.n_els = nm.array( [conn.shape[0] for conn in self.conns] )
        self.n_e_ps = nm.array( [conn.shape[1] for conn in self.conns] )
        self.el_offsets = nm.cumsum( nm.r_[0, self.n_els] )
        self.n_el = nm.sum( self.n_els )

    def _set_data(self, coors, ngroups, conns, mat_ids, descs):
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
        """
        self.coors = nm.ascontiguousarray(coors)

        if ngroups is None:
            self.ngroups = nm.zeros((self.coors.shape[0],), dtype=nm.int32)

        else:
            self.ngroups = nm.ascontiguousarray(ngroups)

        self.conns = [nm.asarray(conn, dtype=nm.int32) for conn in conns]
        self.mat_ids = [nm.asarray(mat_id, dtype=nm.int32) for mat_id in mat_ids]
        self.descs = descs

    def _append_region_faces(self, region):
        from sfepy.fem.domain import dm_create_list

        domain = region.domain
        data, _, _ = dm_create_list(domain.groups, domain.mesh.n_nod, 1, 0)

        for ig in region.igs:
            faces = region.get_faces(ig)
            fdata = data[faces]

            i3 = nm.where(fdata[:,-1] == -1)[0]
            i4 = nm.where(fdata[:,-1] != -1)[0]

            if i3.size:
                self.descs.append('2_3')
                self.mat_ids.append(fdata[i3,0] + 1)
                self.conns.append(fdata[i3,-4:-1].copy())

            if i4.size:
                self.descs.append('2_4')
                self.mat_ids.append(fdata[i4,0] + 1)
                self.conns.append(fdata[i4,-4:].copy())

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
            io = MeshIO.any_from_filename( filename )

        if coors is None:
            coors = self.coors

        if igs is None:
            igs = range( len( self.conns ) )

        aux_mesh = Mesh.from_data( self.name, coors, self.ngroups,
                                   self.conns, self.mat_ids, self.descs, igs )
        io.set_float_format( float_format )
        io.write( filename, aux_mesh, out, **kwargs )

    ##
    # 23.05.2007, c
    def get_bounding_box( self ):
        return nm.vstack( (nm.amin( self.coors, 0 ), nm.amax( self.coors, 0 )) )

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
            aux = remap[conn]
            ii = nm.unique1d(nm.where(aux == -1)[0])
            ii = nm.setdiff1d(nm.arange(conn.shape[0], dtype=nm.int32), ii)
            conns.append(aux[ii])
            mat_ids.append(self.mat_ids[ig][ii])
        self.conns = conns
        self.mat_ids = mat_ids

        self._set_shape_info()

    ##
    # c: 18.01.2008, r: 18.01.2008
    def transform_coors( self, mtx_t, ref_coors = None ):
        """x = T * x."""
        if ref_coors is None:
            ref_coors = self.coors

        self.coors[:] = nm.dot( ref_coors, mtx_t.T )

    def create_conn_graph(self, verbose=True):
        """
        Create a graph of mesh connectivity.

        Returns
        -------
        graph : csr_matrix
            The mesh connectivity graph as a SciPy CSR matrix.    
        """
        from extmods.fem import raw_graph

        shape = (self.n_nod, self.n_nod)
        output('graph shape:', shape, verbose=verbose)
        if nm.prod(shape) == 0:
            output('no graph (zero size)!', verbose=verbose)
            return None

        output('assembling mesh graph...', verbose=verbose)
        tt = time.clock()

        ret, prow, icol = raw_graph(int(shape[0]), int(shape[1]),
                                    len(self.conns), self.conns, self.conns )
        output('...done in %.2f s' % (time.clock() - tt), verbose=verbose)
        nnz = prow[-1]
        output('graph nonzeros: %d (%.2e%% fill)' \
               % (nnz, float(nnz) / nm.prod(shape)))
	
        data = nm.ones((nnz,), dtype=nm.bool)
        graph = sp.csr_matrix((data, icol, prow), shape)

        return graph
