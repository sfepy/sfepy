from sfepy.base.base import *
import sfepy.base.la as la
from meshio import MeshIO

import os.path as op

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

##
# 23.05.2007, updated from matlab version
# 24.05.2007
# 25.05.2007
# 28.05.2007
def merge_mesh( x1, conns1, x2, conns2, cmap, eps = 1e-8 ):
    nc = cmap.shape[0]
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

    conns = []
    for ii in xrange( len( conns1 ) ):
        conn = nm.vstack( (conns1[ii], remap[conns2[ii]]) )
        conns.append( conn )
    
    return xx, conns

##
# 24.05.2007, c
def make_mesh( coor, conns, mesh_in ):
    nod0 = nm.c_[coor, nm.zeros( (coor.shape[0],), dtype = nm.float64 )]

    mat_ids = []
    for ii, conn in enumerate( conns ):
        mat_id = nm.empty( (conn.shape[0],), dtype = nm.int32 )
        mat_id.fill( mesh_in.mat_ids[ii][0] )
        mat_ids.append( mat_id )
        
    mesh_out = Mesh.from_data( 'merged mesh', nod0, conns,
                             mat_ids, mesh_in.descs )
    return mesh_out

def make_inverse_connectivity( conns, n_nod, combine_groups = False ):
    """
    For each mesh node referenced in the connectivity conns, make a list of
    elements it belongs to. If combine_groups is True, elements are referenced
    by (ig, iel), otherwise by iel only.
    """
    if combine_groups:
        iconn = [[] for ii in xrange( n_nod )]
        for ig, conn in enumerate( conns ):
            for iel, row in enumerate( conn ):
                for node in row:
                    iconn[node].append( (ig, iel) )
        return iconn

    else:
        iconns = []
        for ig, conn in enumerate( conns ):
            iconn = [[] for ii in xrange( n_nod )]
            for iel, row in enumerate( conn ):
                for node in row:
                    iconn[node].append( iel )
            iconns.append( iconn )
        return iconns

def find_closest_nodes( x1, x2, num = 1 ):
    """
    For the point x2 find num closest points in x1. Naive algorithm!
    """
    dist = la.norm_l2_along_axis( x1 - x2 )
    ii = dist.argsort()

    out = ii[:num]
    if num == 1:
        out = out[0]

    return out
    
##
# 30.08.2007, c
# 31.08.2007
# 05.09.2007
def find_refinement( coor, conns, fcoor, fconns, eps, check_refined = True ):
    ##
    # Find mesh vertices in the refined mesh.
    vmap = find_map( coor, fcoor )
#    print vmap
    if vmap.shape[0] != coor.shape[0]:
        print 'nonconforming meshes!'
        print vmap.shape, coor.shape, fcoor.shape
        raise ValueError
    ii = nm.argsort( vmap[:,0] )
    vmap = vmap[ii,1]

    ##
    # Inverted connectivity for the refined mesh (= eonlist in mesh_graph()!).
    # Elements are numbered locally per group.
    print fcoor.shape[0]
    iconns = make_inverse_connectivity( fconns, fcoor.shape[0] )
##    print iconns

    ##
    # Create element -> fine elements map.
    dim = coor.shape[1]
    emaps = []
    for ig, conn in enumerate( conns ):
        emap = []
        iconn = iconns[ig]
        fconn = fconns[ig]

        n_ep = fconn.shape[1]
        for iel, row in enumerate( conn ):
            fels = []
#            print '*** iel:', iel, row

            seeds = vmap[row]
            while seeds.size:
#                print 'seeds:', seeds

                new_seeds = []
                for fnode in seeds:
#                    print fnode
                    pfels = [fel for fel in iconn[fnode] if not fel in fels]
                    if not pfels:
                        continue
                    pfels = nm.array( pfels )
#                    print pfels
                    fnodes = nm.array( fconn[pfels,:] )
#                    print fnodes
##                 print fcoor[fnodes,:]

                    # Use average element coordinate (simplex -> always convex).
                    eflag = la.points_in_simplex( fcoor[fnodes,:].sum( 1 ) / n_ep,
                                                coor[row,:], eps )
#                    print eflag
                    if eflag.any():
                        fels.append( pfels[eflag] )
                        new_seeds.append( fnodes[eflag] )
#                    print pfels[eflag], '->', fels

#                print 'new:', new_seeds
                seeds = nm.setdiff1d( nm.asarray( new_seeds ).ravel(), seeds )
#                pause()
            emap.append( nm.array( fels ).ravel() )
##             print '->>', iel, fels
##             pause()
##         print '->>>', emap
        if check_refined:
            n_ref = 2**dim
            for row in emap:
                assert_( len( row ) == n_ref )
        emaps.append( emap )
##         pause()

    iemaps = []
    n_f_els = []
    for ig, emap in enumerate( emaps ):
        iemap = nm.empty( (fconns[ig].shape[0],), dtype = nm.int32 )
        iemap.fill( -1 )
        n_f_el = nm.empty( (len( emap ),), dtype = nm.int32 )
        for iel, row in enumerate( emap ):
            for fiel in row:
                iemap[fiel] = iel
            n_f_el[iel] = len( row )
        iemaps.append( iemap )
        n_f_els.append( n_f_el )

    return emaps, iemaps, n_f_els

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

    Example of creating and working with a mesh:

    >>> from sfepy.fem import Mesh
    >>> m = Mesh.from_file("database/simple.vtk")
    sfepy: reading mesh (database/simple.vtk)...
    sfepy: ...done in 0.06 s
    >>> m.nod0
    array([[  1.00000000e-01,   2.00000000e-02,  -1.22460635e-18,
              0.00000000e+00],
           [  1.00000000e-01,   1.80193774e-02,   8.67767478e-03,
              0.00000000e+00],
           [  1.00000000e-01,   1.24697960e-02,   1.56366296e-02,
              0.00000000e+00],
           ..., 
           [  8.00298527e-02,   5.21598617e-03,  -9.77772215e-05,
              0.00000000e+00],
           [  7.02544004e-02,   3.61610291e-04,  -1.16903153e-04,
              0.00000000e+00],
           [  3.19633596e-02,  -1.00335972e-02,   9.60460305e-03,
              0.00000000e+00]])
    >>> m.conns
    [array([[ 28,  60,  45,  29],
           [ 28,  60,  57,  45],
           [ 28,  57,  27,  45],
           ..., 
           [353, 343, 260, 296],
           [353, 139, 181, 140],
           [353, 295, 139, 140]])]
    >>> m.mat_ids
    [array([6, 6, 6, ..., 6, 6, 6])]
    >>> m.descs
    ['3_4']
    >>> m
    Mesh:database/simple
    >>> print m
    Mesh:database/simple
      n_e_ps:
        [4]
      dim:
        3
      n_el:
        1348
      name:
        database/simple
      descs:
        ['3_4']
      n_nod:
        354
      mat_ids:
        [array([6, 6, 6, ..., 6, 6, 6])]
      n_els:
        [1348]
      nod0:
        [[  1.00000000e-01   2.00000000e-02  -1.22460635e-18   0.00000000e+00]
         [  1.00000000e-01   1.80193774e-02   8.67767478e-03   0.00000000e+00]
         [  1.00000000e-01   1.24697960e-02   1.56366296e-02   0.00000000e+00]
         ..., 
         [  8.00298527e-02   5.21598617e-03  -9.77772215e-05   0.00000000e+00]
         [  7.02544004e-02   3.61610291e-04  -1.16903153e-04   0.00000000e+00]
         [  3.19633596e-02  -1.00335972e-02   9.60460305e-03   0.00000000e+00]]
      io:
        None
      conns:
        [array([[ 28,  60,  45,  29],
               [ 28,  60,  57,  45],
               [ 28,  57,  27,  45],
               ..., 
               [353, 343, 260, 296],
               [353, 139, 181, 140],
               [353, 295, 139, 140]])]
      setup_done:
        0
      el_offsets:
        [   0 1348]


    The Mesh().nod0 is an array of nodes and Mesh().conns is the list of
    elements of each type (see Mesh().desc), so for example if you want to know
    the coordinates of the nodes of the fifth finite element of the type 3_4 do:

    In [1]: a.descs
    Out[1]: ['3_4']

    So now you know that the finite elements of the type 3_4 are in a.conns[0]:

    In [2]: [a.nod0[n] for n in a.conns[0][4]]
    Out[2]:
    [array([ 3.0877856 , -4.40913864, -1.58148163,  0.        ]),
     array([ 3.28954489, -4.53265378, -2.07926241,  0.        ]),
     array([ 3.00343981, -4.09445003, -2.14632505,  0.        ]),
     array([ 3.51217117, -4.05946689, -1.68843294,  0.        ])]

    The element ids are of the form "<dimension>_<number of nodes>", i.e.:

    2_2 ... line
    2_3 ... triangle
    2_4 ... quadrangle
    3_2 ... line
    3_4 ... tetrahedron
    3_8 ... hexahedron

    """
    ##
    # 19.01.2005, c
    # 03.03.2005
    # 08.03.2005
    # 05.10.2005
    # 04.08.2006
    # 29.08.2007
    def from_surface( surf_faces, mesh_in ):

        inod = la.as_unique_set( surf_faces )
        n_nod = len( inod )
        n_nod_m, n_col = mesh_in.nod0.shape

        aux = nm.arange( n_nod )
        remap = nm.zeros( (n_nod_m,), nm.int32 )
        remap[inod] = aux

        mesh = Mesh( mesh_in.name + "_surf" )

        mesh.nod0 = mesh_in.nod0[inod,:n_col]

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

    ##
    # c: 25.01.2006, r: 23.06.2008
    def from_file( filename = None, io = 'auto' ):
        """passing *MeshIO instance has precedence over filename"""
        if io == 'auto':
            if filename is None:
                output( 'filename or io must be specified!' )
                raise ValueError
            else:
                io = MeshIO.any_from_filename( filename )
                if isinstance( filename, file ):
                    trunk = 'from_descriptor'
                else:
                    trunk = op.splitext( filename )[0]
        else:
            trunk = io.filename

        output( 'reading mesh (%s)...' % (filename) )
        tt = time.clock()
        mesh = Mesh( trunk )
        mesh = io.read( mesh )
        output( '...done in %.2f s' % (time.clock() - tt) )
        mesh._set_shape_info()
        return mesh
    from_file = staticmethod( from_file )

    ##
    # c: 17.02.2006, r: 28.04.2008
    def from_region( region, mesh_in, ed = None, fa = None, localize = None ):
        mesh = Mesh( mesh_in.name + "_reg" )
        mesh.nod0 = mesh_in.nod0.copy()
        
        mesh.conns = []
        mesh.descs = []
        mesh.mat_ids = []
        if region.has_cells():
            for ig in region.igs:
                mesh.descs.append( mesh_in.descs[ig] )
                els = region.get_cells( ig )
                mesh.mat_ids.append( mesh_in.mat_ids[ig][els,:].copy() )
                mesh.conns.append( mesh_in.conns[ig][els,:].copy() )

        if ed is not None:
            for ig in region.igs:
                edges = region.get_edges( ig )
                mesh.descs.append( '1_2' )
                mesh.mat_ids.append( ed.data[edges,0] + 1 )
                mesh.conns.append( ed.data[edges,-2:].copy() )

        if fa is not None:
            for ig in region.igs:
                faces = region.get_faces( ig )
                fdata = fa.data[faces]
                i3 = nm.where( fdata[:,-1] == -1 )[0]
                i4 = nm.where( fdata[:,-1] != -1 )[0]
                if i3.size:
                    mesh.descs.append( '2_3' )
                    mesh.mat_ids.append( fdata[i3,0] + 1 )
                    mesh.conns.append( fdata[i3,-4:-1].copy() )
                if i4.size:
                    mesh.descs.append( '2_4' )
                    mesh.mat_ids.append( fdata[i4,0] + 1 )
                    mesh.conns.append( fdata[i4,-4:].copy() )

        if (ed is not None) or (fa is not None):
            mesh.descs.append( {2 : '2_3', 3 : '3_4'}[mesh_in.dim] )
            mesh.mat_ids.append( -nm.ones_like( region.all_vertices ) )
            mesh.conns.append( make_point_cells( region.all_vertices, mesh_in.dim ) )

        if localize:
            mesh.localize( region.all_vertices )
        mesh._set_shape_info()
        return mesh
    from_region = staticmethod( from_region )

    ##
    # c: 02.01.2008, r: 02.01.2008
    def from_region_and_field( region, field ):
        mesh, ed, fa = field.domain.mesh, field.domain.ed, field.domain.fa
        mesh = Mesh.from_region( region, mesh, ed, fa )
        mesh.name = mesh.name + '_field'

        nodes = region.get_field_nodes( field, merge = True )

        aux = field.get_extra_nodes_as_simplices( nodes )
        mesh.nod0 = field.aps.coors
        mesh.descs.append( aux[0] )
        mesh.mat_ids.append( aux[1] )
        mesh.conns.append( aux[2] )

        mesh.localize( nodes )
        mesh._set_shape_info()
        return mesh
    from_region_and_field = staticmethod( from_region_and_field )

    ##
    # c: 21.02.2007, r: 08.02.2008
    def from_data( name, coors, conns, mat_ids, descs, igs = None ):
        if igs is None:
            igs = range( len( conns ) )
        mesh = Mesh( name = name,
                     nod0 = coors,
                     conns = [conns[ig] for ig in igs],
                     mat_ids = [mat_ids[ig] for ig in igs],
                     descs = [descs[ig] for ig in igs] )
        mesh._set_shape_info()
        return mesh
    from_data = staticmethod( from_data )
        

    ##
    # 22.02.2005
    # 16.06.2005
    # 26.09.2006
    def __init__( self, name = 'mesh', **kwargs ):
        Struct.__init__( self, **kwargs )
        self.name = name
        self.io = None
        self.setup_done = 0

    ##
    # 04.08.2006, c
    # 29.09.2006
    def _set_shape_info( self ):
        self.n_nod = self.nod0.shape[0]
        self.dim = self.nod0.shape[1] - 1
        self.n_els = nm.array( [conn.shape[0] for conn in self.conns] )
        self.n_e_ps = nm.array( [conn.shape[1] for conn in self.conns] )
        self.el_offsets = nm.cumsum( nm.r_[0, self.n_els] )
        self.n_el = nm.sum( self.n_els )

    ##
    # c: 15.02.2008, r: 15.02.2008
    def _set_data( self, coors, conns, mat_ids, descs ):
        self.nod0 = coors
        self.conns = conns
        self.mat_ids = mat_ids
        self.descs = descs
        
    ##
    # c: 23.01.2006, r: 23.06.2008
    def write( self, filename = None, io = None,
               coors = None, igs = None, out = None, float_format = None,
               **kwargs ):
        """Write mesh + optional results in 'out'.

        'io' == 'auto' respects the extension of 'filename'
        'coors' can be used instead of mesh coordinates,
        providing 'igs' filters some groups only"""
        if filename is None:
            filename = self.name + '.mesh'

        if io is None:
            io = self.io
        else:
            if io == 'auto':
                io = MeshIO.any_from_filename( filename )

        if coors is None:
            coors = self.nod0

        if igs is None:
            igs = range( len( self.conns ) )

        aux_mesh = Mesh.from_data( self.name, coors,
                                 self.conns, self.mat_ids, self.descs, igs )
        io.set_float_format( float_format )
        io.write( filename, aux_mesh, out, **kwargs )

    ##
    # 23.05.2007, c
    def get_bounding_box( self ):
        return nm.array( [nm.amin( self.nod0[:,:-1], 0 ),
                          nm.amax( self.nod0[:,:-1], 0 )] )


    ##
    # c: 02.01.2008, r: 02.01.2008
    def localize( self, inod ):
        """Strips nodes not in inod and remaps connectivities.
        TODO: fix the case when remap[conn] contains -1..."""
        remap = nm.empty( (self.nod0.shape[0],), dtype = nm.int32 )
        remap.fill( -1 )
        remap[inod] = nm.arange( inod.shape[0], dtype = nm.int32 )

        self.nod0 = self.nod0[inod]
        conns = []
        for conn in self.conns:
            conns.append( remap[conn] )
        self.conns = conns


    ##
    # c: 18.01.2008, r: 18.01.2008
    def transform_coords( self, mtx_t, ref_coords = None ):
        """x = T * x."""
        if ref_coords is None:
            ref_coords = self.nod0[:,:-1]

        self.nod0[:,:-1] = nm.dot( ref_coords, mtx_t.T )
