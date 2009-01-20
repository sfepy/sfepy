from sfepy.base.base import *
from sfepy.base.ioutils import skip_read_line, read_token, read_array, read_list, pt
import sfepy.base.la as la
from sfepy.base.progressbar import MyBar
import os.path as op

supported_formats = {
    '.mesh' : 'medit',
    '.vtk'  : 'vtk',
    '.node' : 'tetgen',
    '.txt'  : 'comsol',
    '.h5'   : 'hdf5',
}

##
# c: 15.02.2008, r: 15.02.2008
def sort_by_mat_id( conns_in ):

    # Sort by mat_id within a group, preserve order.
    conns = []
    mat_ids = []
    for ig, conn in enumerate( conns_in ):
        ii = nm.argsort( conn[:,-1], kind = 'mergesort' )
        conn = conn[ii]

        conns.append( conn[:,:-1].copy() )
        mat_ids.append( conn[:,-1].copy() )
    return conns, mat_ids

##
# conns_in must be sorted by mat_id within a group!
# c: 16.06.2005, r: 15.02.2008
def split_by_mat_id( conns_in, mat_ids_in, descs_in ):

    conns = []
    mat_ids = []
    descs = []

    for ig, conn in enumerate( conns_in ):
        one = nm.array( [-1], nm.int32 )
        ii = la.diff( nm.concatenate( (one, mat_ids_in[ig], one) ) ).nonzero()[0]
        n_gr = len( ii ) - 1;
#        print ii, n_gr
        for igr in range( 0, n_gr ):
            conns.append( conn[ii[igr]:ii[igr+1],:].copy() )
            mat_ids.append( mat_ids_in[ig][ii[igr]:ii[igr+1]] )
            descs.append( descs_in[ig] )

    return (conns, mat_ids, descs)


##
# 12.10.2005, c
def write_bb( fd, array, dtype ):

    fd.write( '3 %d %d %d\n' % (array.shape[1], array.shape[0], dtype) )
    format = ' '.join( ['%.5e'] * array.shape[1] + ['\n'] )

    for row in array:
        fd.write( format % tuple( row ) )

##
# c: 03.10.2005, r: 08.02.2008
def join_conn_groups( conns, descs, mat_ids, concat = False ):
    """Join groups of the same element type."""

    el = dict_from_keys_init( descs, list )
    for ig, desc in enumerate( descs ):
        el[desc].append( ig )
    groups = [ii for ii in el.values() if ii]
##     print el, groups

    descs_out, conns_out, mat_ids_out = [], [], []
    for group in groups:
        n_ep = conns[group[0]].shape[1]

        conn = nm.zeros( (0, n_ep), nm.int32 )
        mat_id = nm.zeros( (0,), nm.int32 )
        for ig in group:
            conn = nm.concatenate( (conn, conns[ig]) )
            mat_id = nm.concatenate( (mat_id, mat_ids[ig]) )

        if concat:
            conn = nm.concatenate( (conn, mat_id[:,nm.newaxis]), 1 )
        else:
            mat_ids_out.append( mat_id )
        conns_out.append( conn )
        descs_out.append( descs[group[0]] )

    if concat:
        return conns_out, descs_out
    else:
        return conns_out, descs_out, mat_ids_out

##
# c: 05.02.2008
class MeshIO( Struct ):
    """
    The abstract class for importing and exporting meshes.

    Read the docstring of the Mesh() class. Basically all you need to do is to
    implement the read() method:

    def read(self, mesh, **kwargs):
        nodes = ...
        conns = ...
        mat_ids = ...
        descs = ...
        mesh._set_data(nodes, conns, mat_ids, descs)
        return mesh

    See the Mesh() class's docstring how the nodes, conns, mat_ids and descs
    should look like. You just need to read them from your specific format from
    disk.

    To write a mesh to disk, just implement the write() method and use the
    information from the mesh instance (e.g. nodes, conns, mat_ids and descs)
    to construct your specific format.

    """
    format = None

    ##
    # c: 05.02.2008, r: 05.02.2008
    def __init__( self, filename, **kwargs ):
        Struct.__init__( self, filename = filename, **kwargs )
        self.set_float_format()

    ##
    # c: 03.07.2008, r: 03.07.2008
    def read_dimension( self, ret_fd = False ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

    ##
    # c: 22.07.2008, r: 22.07.2008
    def read_boundin_box( self, ret_fd = False, ret_dim = False ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

    ##
    # c: 05.02.2008, r: 26.03.2008
    def read( self, mesh, *args, **kwargs ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

    ##
    # c: 05.02.2008, r: 26.03.2008
    def write( self, filename, mesh, *args, **kwargs ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

    def set_float_format( self, format = None ):
        self.float_format = get_default( format, '%e' )

    def get_vector_format( self, dim ):
        return ' '.join( [self.float_format] * dim )
            
##
# c: 05.02.2008
class MeditMeshIO( MeshIO ):
    format = 'medit'

    ##
    # c: 03.07.2008, r: 10.07.2008
    def read_dimension( self, ret_fd = False ):
        fd = open( self.filename, 'r' )
        while 1:
            try:
                line = fd.readline()
            except:
                output( "reading " + fd.name + " failed!" )
                raise
            if len( line ) == 1: continue
            if line[0] == '#': continue
            aux = line.split()
            if aux[0] == 'Dimension':
                if len( aux ) == 2:
                    dim = int( aux[1] )
                else:
                    dim = int( fd.readline() )
                break

        if ret_fd:
            return dim, fd
        else:
            fd.close()
            return dim

    ##
    # c: 22.07.2008
    def read_boundin_box( self, ret_fd = False, ret_dim = False ):
        fd = open( self.filename, 'r' )

        while 1:
            try:
                line = fd.readline()
            except:
                output( "reading " + fd.name + " failed!" )
                raise
            if len( line ) == 0: break
            if len( line ) == 1: continue
            if line[0] == '#': continue
            aux = line.split()
            if (aux[0] == 'Dimension'):
                if len( aux ) == 2:
                    dim = int( aux[1] )
                else:
                    dim = int( fd.readline() )            
            elif (aux[0] == 'Vertices'):
                num = int( read_token( fd ) )
                nod = read_array( fd, num, dim + 1, nm.float64 )
                break

        bbox = [ [nod[0][0]]*2, [nod[0][1]]*2, [nod[0][2]]*2 ]

        for inod in nod[1:]:
            for idim in range( dim ):
                if inod[idim] < bbox[idim][0]:
                    bbox[idim][0] = inod[idim]
                if inod[idim] > bbox[idim][1]:
                    bbox[idim][1] = inod[idim]


        if ret_dim:
            if ret_fd:
                return bbox, dim, fd
            else:
                fd.close()
                return bbox, dim
        else:
            if ret_fd:
                return bbox, fd
            else:
                fd.close()
                return bbox

    def read( self, mesh, **kwargs ):
        dim, fd  = self.read_dimension( ret_fd = True )

        conns_in = []
        descs = []
        while 1:
            try:
                line = fd.readline()
                if (len( line ) == 0): break
                if len( line ) == 1: continue
            except EOFError:
                break
            except:
                output( "reading " + fd.name + " failed!" )
                raise
            ls = line.strip()
            if (ls == 'Vertices'):
                num = int( read_token( fd ) )
                nod = read_array( fd, num, dim + 1, nm.float64 )
    ##                 print nod
            elif (ls == 'Tetrahedra'):
                num = int( read_token( fd ) )
                conns_in.append( read_array( fd, num, 5, nm.int32 ) )
                conns_in[-1][:,:-1] -= 1
                descs.append( '3_4' )
            elif (ls == 'Hexahedra'):
                num = int( read_token( fd ) )
                conns_in.append( read_array( fd, num, 9, nm.int32 ) )
                conns_in[-1][:,:-1] -= 1
                descs.append( '3_8' )
            elif (ls == 'Triangles'):
                num = int( read_token( fd ) )
                conns_in.append( read_array( fd, num, 4, nm.int32 ) )
                conns_in[-1][:,:-1] -= 1
                descs.append( '2_3' )
            elif (ls == 'Quadrilaterals'):
                num = int( read_token( fd ) )
                conns_in.append( read_array( fd, num, 5, nm.int32 ) )
                conns_in[-1][:,:-1] -= 1
                descs.append( '2_4' )
            elif ls == 'End':
                break
            elif line[0] == '#':
                continue
            else:
                msg = "corrupted file (line '%s')!" % line
                raise ValueError( msg )
        fd.close()

        conns_in, mat_ids = sort_by_mat_id( conns_in )

        # Detect wedges and pyramides -> separate groups.
        if ('3_8' in descs):
            ic = descs.index( '3_8' )

            conn_in = conns_in.pop( ic )
            flag = nm.zeros( (conn_in.shape[0],), nm.int32 )
            for ii, el in enumerate( conn_in ):
                if (el[4] == el[5]):
                    if (el[5] == el[6]):
                        flag[ii] = 2
                    else:
                        flag[ii] = 1

            conn = []
            desc = []

            ib = nm.where( flag == 0 )[0]
            if (len( ib ) > 0):
                conn.append( conn_in[ib] )
                desc.append( '3_8' )

            iw = nm.where( flag == 1 )[0]
            if (len( iw ) > 0):
                ar = nm.array( [0,1,2,3,4,6,8], nm.int32 )
                conn.append( la.rect( conn_in, iw, ar ) )
                desc.append( '3_6' )

            ip = nm.where( flag == 2 )[0]
            if (len( ip ) > 0):
                ar = nm.array( [0,1,2,3,4,8], nm.int32 )
                conn.append( la.rect( conn_in, ip, ar ) )
                desc.append( '3_5' )

##             print "brick split:", ic, ":", ib, iw, ip, desc

            conns_in[ic:ic] = conn
            del( descs[ic] )
            descs[ic:ic] = desc

        conns, mat_ids, descs = split_by_mat_id( conns_in, mat_ids, descs )
        mesh._set_data( nod, conns, mat_ids, descs )

        return mesh

    def write( self, filename, mesh, out = None ):
        fd = open( filename, 'w' )

        nod = mesh.nod0
        conns, desc = join_conn_groups( mesh.conns, mesh.descs,
                                      mesh.mat_ids, concat = True )

        n_nod, dim = nod.shape
        dim -= 1

        fd.write( "MeshVersionFormatted 1\nDimension %d\n" % dim )

        fd.write( "Vertices\n%d\n" % n_nod )
        format = self.get_vector_format( dim ) + ' %d\n'
        for ii in range( n_nod ):
            nn = nod[ii]
            fd.write( format % tuple( nn ) )

        for ig, conn in enumerate( conns ):
            if (desc[ig] == "1_2"):
                fd.write( "Edges\n%d\n" % conn.shape[0] )
                for ii in range( conn.shape[0] ):
                    nn = conn[ii] + 1
                    fd.write( "%d %d %d\n" \
                              % (nn[0], nn[1], nn[2] - 1) )
            elif (desc[ig] == "2_4"):
                fd.write( "Quadrilaterals\n%d\n" % conn.shape[0] )
                for ii in range( conn.shape[0] ):
                    nn = conn[ii] + 1
                    fd.write( "%d %d %d %d %d\n" \
                              % (nn[0], nn[1], nn[2], nn[3], nn[4] - 1) )
            elif (desc[ig] == "2_3"):
                fd.write( "Triangles\n%d\n" % conn.shape[0] )
                for ii in range( conn.shape[0] ):
                    nn = conn[ii] + 1
                    fd.write( "%d %d %d %d\n" % (nn[0], nn[1], nn[2], nn[3] - 1) )
            elif (desc[ig] == "3_4"):
                fd.write( "Tetrahedra\n%d\n" % conn.shape[0] )
                for ii in range( conn.shape[0] ):
                    nn = conn[ii] + 1
                    fd.write( "%d %d %d %d %d\n"
                              % (nn[0], nn[1], nn[2], nn[3], nn[4] - 1) )
            elif (desc[ig] == "3_8"):
                fd.write( "Hexahedra\n%d\n" % conn.shape[0] )
                for ii in range( conn.shape[0] ):
                    nn = conn[ii] + 1
                    fd.write( "%d %d %d %d %d %d %d %d %d\n"
                              % (nn[0], nn[1], nn[2], nn[3], nn[4], nn[5],
                                 nn[6], nn[7], nn[8] - 1) )
            else:
                print 'unknown element type!', desc[ig]
                raise ValueError

        fd.close()

        if out is not None:
            for key, val in out.iteritems():
                raise NotImplementedError


vtk_header = r"""# vtk DataFile Version 2.0
generated by %s
ASCII
DATASET UNSTRUCTURED_GRID
"""
vtk_cell_types = {'2_2' : 3, '2_4' : 9, '2_3' : 5,
                '3_2' : 3, '3_4' : 10, '3_8' : 12 }
vtk_dims = {3 : 2, 9 : 2, 5 : 2, 3 : 3, 10 : 3, 12 : 3}
vtk_inverse_cell_types = {(3, 2) : '2_2', (9, 2) : '2_4', (5, 2) : '2_3',
                       (3, 3) : '3_2', (10, 3) : '3_4', (12, 3) : '3_8' }

##
# c: 05.02.2008
class VTKMeshIO( MeshIO ):
    format = 'vtk'

    ##
    # c: 03.07.2008, r: 15.07.2008
    def read_dimension( self, ret_fd = False ):
        fd = open( self.filename, 'r' )
        while 1:
            try:
                line = fd.readline().split()
                if not line: continue
                if line[0] == 'CELL_TYPES':
                    cell_types = read_array( fd, 1, -1, nm.int32 )
                    dim = vtk_dims[cell_types[0,0]]
                    break
            except:
                output( "reading " + fd.name + " failed!" )
                raise

        if ret_fd:
            return dim, fd
        else:
            fd.close()
            return dim

    ##
    # c: 22.07.2008
    def read_boundin_box( self, ret_fd = False, ret_dim = False ):
        fd = open( self.filename, 'r' )
        while 1:
            try:
                line = fd.readline().split()
                if line[0] == 'POINTS':
                    nod = read_array( fd, 1, -1, nm.float64 )
                    dim = nod.shape[1]
                    break
            except:
                output( "reading " + fd.name + " failed!" )
                raise

        bbox = [ [nod[0][0]]*2, [nod[0][1]]*2, [nod[0][2]]*2 ]

        for inod in nod[1:]:
            for idim in range( dim ):
                if inod[idim] < bbox[idim][0]:
                    bbox[idim][0] = inod[idim]
                if inod[idim] > bbox[idim][1]:
                    bbox[idim][1] = inod[idim]

        if ret_dim:
            if ret_fd:
                return bbox, dim, fd
            else:
                fd.close()
                return bbox, dim
        else:
            if ret_fd:
                return bbox, fd
            else:
                fd.close()
                return bbox

    ##
    # c: 05.02.2008, r: 10.07.2008
    def read( self, mesh, **kwargs ):
        fd = open( self.filename, 'r' )
        mode = 'header'
        mode_status = 0
        nod = conns = desc = mat_id = None
        while 1:
            try:
                line = fd.readline()
                if len( line ) == 0: break
                elif len( line ) == 1: continue
                if line[0] == '#': continue
            except EOFError:
                break
            except:
                output( "reading " + fd.name + " failed!" )
                raise

            if mode == 'header':
                if mode_status == 0:
                    if line.strip() == 'ASCII':
                        mode_status = 1
                elif mode_status == 1:
                    if line.strip() == 'DATASET UNSTRUCTURED_GRID':
                        mode_status = 0
                        mode = 'points'

            elif mode == 'points':
                line = line.split()
                if line[0] == 'POINTS':
                    n_nod = int( line[1] )
                    nod = read_array( fd, n_nod, -1, nm.float64 )
                    mode = 'cells'

            elif mode == 'cells':
                line = line.split()
                if line[0] == 'CELLS':
                    n_el, n_val = map( int, line[1:3] )
                    raw_conn = read_list( fd, n_val, int )
                    mode = 'cell_types'

            elif mode == 'cell_types':
                line = line.split()
                if line[0] == 'CELL_TYPES':
                    assert_( int( line[1] ) == n_el )
                    cell_types = read_array( fd, n_el, -1, nm.int32 )
                    mode = 'mat_id'

            elif mode == 'mat_id':
                if mode_status == 0:
                    line = line.split()
                    if line[0] == 'CELL_DATA':
                        assert_( int( line[1] ) == n_el )
                        mode_status = 1
                elif mode_status == 1:
                    if line.strip() == 'SCALARS mat_id float 1':
                        mode_status = 2
                elif mode_status == 2:
                    if line.strip() == 'LOOKUP_TABLE default':
                        mat_id = read_list( fd, n_el, int )
                        mode_status = 0
                        mode = 'finished'
            elif mode == 'finished':
                break
        fd.close()
 
        if mat_id is None:
            mat_id = [[0]] * n_el

        dim = vtk_dims[cell_types[0,0]]
        if dim == 3:
            nod = nm.concatenate( (nod, nm.zeros( (n_nod,1), dtype = nm.int32 ) ),
                                  1 )
        else:
            nod[:,2] = 0.0
        nod = nm.ascontiguousarray( nod )

        dim = nod.shape[1] - 1
        cell_types = cell_types.squeeze()

        dconns = {}
        for iel, row in enumerate( raw_conn ):
            ct = vtk_inverse_cell_types[(cell_types[iel],dim)]
            dconns.setdefault( ct, [] ).append( row[1:] + mat_id[iel] )

        desc = []
        conns = []
        for key, conn in dconns.iteritems():
            desc.append( key )
            conns.append( nm.array( conn, dtype = nm.int32 ) )

        conns_in, mat_ids = sort_by_mat_id( conns )
        conns, mat_ids, descs = split_by_mat_id( conns_in, mat_ids, desc )
        mesh._set_data( nod, conns, mat_ids, descs )

        return mesh

    def write( self, filename, mesh, out = None ):

        fd = open( filename, 'w' )
        fd.write( vtk_header % op.basename( sys.argv[0] ) )

        n_nod, nc = mesh.nod0.shape
        dim = nc - 1
        sym = dim * (dim + 1) / 2

        fd.write( '\nPOINTS %d float\n' % n_nod )

        aux = mesh.nod0[:,:dim]
        if dim == 2:
            aux = nm.hstack( (aux, nm.zeros( (n_nod, 1), dtype = nm.float64 ) ) )

        format = self.get_vector_format( 3 ) + '\n'
        for row in aux:
            fd.write( format % tuple( row ) )

        n_el, n_els, n_e_ps = mesh.n_el, mesh.n_els, mesh.n_e_ps
        total_size = nm.dot( n_els, n_e_ps + 1 )
        fd.write( '\nCELLS %d %d\n' % (n_el, total_size) )

        ct = []
        for ig, conn in enumerate( mesh.conns ):
            nn = n_e_ps[ig] + 1
            ct += [vtk_cell_types[mesh.descs[ig]]] * n_els[ig]
            format = ' '.join( ['%d'] * nn + ['\n'] )

            for row in conn:
                fd.write( format % ((nn-1,) + tuple( row )) )

        fd.write( '\nCELL_TYPES %d\n' % n_el )
        fd.write( ''.join( ['%d\n' % ii for ii in ct] ) )

        if out is None: return

        point_keys = [key for key, val in out.iteritems() if val.mode == 'vertex']
        if len( point_keys ):
            fd.write( '\nPOINT_DATA %d\n' % n_nod );
        for key in point_keys:
            val = out[key]
            nr, nc = val.data.shape

            if nc == 1:
                fd.write( '\nSCALARS %s float %d\n' % (key, nc) )
                fd.write( 'LOOKUP_TABLE default\n' )

                format = self.float_format + '\n'
                for row in val.data:
                    fd.write( format % row )

            elif nc == dim:
                fd.write( '\nVECTORS %s float\n' % key )
                if dim == 2:
                    aux = nm.hstack( (val.data,
                                      nm.zeros( (nr, 1), dtype = nm.float64 ) ) )
                else:
                    aux = val.data

                format = self.get_vector_format( 3 ) + '\n'
                for row in aux:
                    fd.write( format % tuple( row ) )

            else:
                raise NotImplementedError, nc

        cell_keys = [key for key, val in out.iteritems() if val.mode == 'cell']
        if len( cell_keys ):
            fd.write( '\nCELL_DATA %d\n' % n_el );
        for key in cell_keys:
            val = out[key]
            ne, aux, nr, nc = val.data.shape

            if (nr == 1) and (nc == 1):
                fd.write( '\nSCALARS %s float %d\n' % (key, nc) )
                fd.write( 'LOOKUP_TABLE default\n' )
                format = self.float_format + '\n'
                for row in val.data.squeeze():
                    fd.write( format % row )

            elif (nr == dim) and (nc == 1):
                fd.write( '\nVECTORS %s float\n' % key )
                if dim == 2:
                    aux = nm.hstack( (val.data.squeeze(),
                                      nm.zeros( (ne, 1), dtype = nm.float64 ) ) )
                else:
                    aux = val.data

                format = self.get_vector_format( 3 ) + '\n'
                for row in aux:
                    fd.write( format % tuple( row.squeeze() ) )

            elif (((nr == sym) or (nr == (dim * dim))) and (nc == 1)) \
                     or ((nr == dim) and (nc == dim)):
                # Below not tested!!!
                fd.write( '\nTENSORS %s float\n' % key );
                data = val.data.squeeze()

                if dim == 3:
                    if nr == sym:
                        aux = data[:,[0,3,4,3,1,5,4,5,2]]
                    elif nr == (dim * dim):
                        aux = data[:,[0,3,4,6,1,5,7,8,2]]
                    else:
                        aux = data[:,[0,1,2,3,4,5,6,7,8]]

                else:
                    zz = nm.zeros( (data.shape[0], 1), dtype = nm.float64 );
                    if nr == sym:
                        aux = nm.c_[data[:,[0,2]], zz, data[:,[2,1]],
                                    zz, zz, zz, zz]
                    elif nr == (dim * dim):
                        aux = nm.c_[data[:,[0,2]], zz, data[:,[3,1]],
                                    zz, zz, zz, zz]
                    else:
                        aux = nm.c_[data[:,0,[0,1]], zz, data[:,1,[0,1]],
                                    zz, zz, zz, zz]

                format = self.get_vector_format( 3 )
                format = '\n'.join( [format] * 3 ) + '\n\n';
                for row in aux:
                    fd.write( format % tuple( row ) )

            else:
                raise NotImplementedError, (nr, nc)

        fd.close()

##
# c: 15.02.2008
class TetgenMeshIO( MeshIO ):
    format = "tetgen"

    ##
    # c: 15.02.2008, r: 15.02.2008
    def read( self, mesh, **kwargs ):
        import os
        fname = os.path.splitext(self.filename)[0]
        nodes=self.getnodes(fname+".node", MyBar("       nodes:"))
        elements, regions = self.getele(fname+".ele", MyBar("       elements:"))
        descs = []
        conns = []
        mat_ids = []
        nodes = nm.c_[(nm.array(nodes, dtype = nm.float64),
                       nm.zeros(len(nodes), dtype = nm.float64))].copy()
        elements = nm.array( elements, dtype = nm.int32 )-1
        for key, value in regions.iteritems():
            descs.append( "3_4" )
            mat_ids.append( nm.ones_like(value) * key )
            conns.append( elements[nm.array(value)-1].copy() )

        mesh._set_data( nodes, conns, mat_ids, descs )
        return mesh

    ##
    # c: 15.02.2008, r: 15.02.2008
    @staticmethod
    def getnodes(fnods, up, verbose=True):
        """
        Reads t.1.nodes, returns a list of nodes.

        Example:

        >>> self.getnodes("t.1.node", MyBar("nodes:"))
        [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (-4.0, 0.0, 0.0),
        (0.0, 0.0, 4.0), (0.0, -4.0, 0.0), (0.0, -0.0, -4.0), (-2.0, 0.0,
        -2.0), (-2.0, 2.0, 0.0), (0.0, 2.0, -2.0), (0.0, -2.0, -2.0), (2.0,
        0.0, -2.0), (2.0, 2.0, 0.0), ... ]

        """
        f=open(fnods)
        l=[int(x) for x in f.readline().split()]
        npoints,dim,nattrib,nbound=l
        assert_( dim==3 )
        if verbose: up.init(npoints)
        nodes=[]
        for line in f:
            if line[0]=="#": continue
            l=[float(x) for x in line.split()]
            assert_( int(l[0])==len(nodes)+1 )
            l = l[1:]
            nodes.append(tuple(l))
            if verbose: up.update(len(nodes))
        assert_( npoints==len(nodes) )
        return nodes

    ##
    # c: 15.02.2008, r: 15.02.2008
    @staticmethod
    def getele(fele, up, verbose=True):
        """
        Reads t.1.ele, returns a list of elements.

        Example:

        >>> elements, regions = self.getele("t.1.ele", MyBar("elements:"))
        >>> elements
        [(20, 154, 122, 258), (86, 186, 134, 238), (15, 309, 170, 310), (146,
        229, 145, 285), (206, 207, 125, 211), (99, 193, 39, 194), (185, 197,
        158, 225), (53, 76, 74, 6), (19, 138, 129, 313), (23, 60, 47, 96),
        (119, 321, 1, 329), (188, 296, 122, 322), (30, 255, 177, 256), ...]
        >>> regions
        {100: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
        19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
        55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 7, ...],
        ...}

        """
        f=file(fele)
        l=[int(x) for x in f.readline().split()]
        ntetra,nnod,nattrib=l
        #we have either linear or quadratic tetrahedra:
        assert_( nnod in [4,10] )
        linear= (nnod==4)
        if not linear:
            raise Exception("Only linear tetrahedra reader is implemented")
        if verbose: up.init(ntetra)
        if nattrib!=1:
            raise "tetgen didn't assign an entity number to each element \
(option -A)"
        els=[]
        regions={}
        for line in f:
            if line[0]=="#": continue
            l=[int(x) for x in line.split()]
            assert_( len(l)-2 == 4 )
            els.append((l[1],l[2],l[3],l[4]))
            regionnum=l[5]
            if regionnum==0:
                print "see %s, element # %d"%(fele,l[0])
                raise "there are elements not belonging to any physical entity"
            if regions.has_key(regionnum):
                regions[regionnum].append(l[0])
            else:
                regions[regionnum]=[l[0]]
            assert_( l[0]==len(els) )
            if verbose: up.update(l[0])
        return els, regions

    ##
    # c: 26.03.2008, r: 26.03.2008
    def write( self, filename, mesh, out = None ):
        raise NotImplementedError

    def read_dimension(self):
        # TetGen only supports 3D mesh
        return 3

    ##
    # c: 22.07.2008
    def read_boundin_box( self ):
        raise NotImplementedError


##
# c: 20.03.2008
class ComsolMeshIO( MeshIO ):
    format = 'comsol'

    ##
    # c: 20.03.2008, r: 20.03.2008
    def _read_commented_int( self ):
        return int( skip_read_line( self.fd ).split( '#' )[0] )

    ##
    # c: 20.03.2008, r: 20.03.2008
    def read( self, mesh, **kwargs ):

        self.fd = fd = open( self.filename, 'r' )
        mode = 'header'

        nod = conns = desc = None
        while 1:
            if mode == 'header':
                line = skip_read_line( fd )

                n_tags = self._read_commented_int()
                for ii in xrange( n_tags ):
                    skip_read_line( fd )
                n_types = self._read_commented_int()
                for ii in xrange( n_types ):
                    skip_read_line( fd )

                skip_read_line( fd )
                assert_( skip_read_line( fd ).split()[1] == 'Mesh' )
                skip_read_line( fd )
                dim = self._read_commented_int()
                assert_( (dim == 2) or (dim == 3) )
                n_nod = self._read_commented_int()
                i0 = self._read_commented_int()
                mode = 'points'

            elif mode == 'points':
                nod = read_array( fd, n_nod, dim, nm.float64 )
                mode = 'cells'

            elif mode == 'cells':

                n_types = self._read_commented_int()
                conns = []
                descs = []
                mat_ids = []
                for it in xrange( n_types ):
                    t_name = skip_read_line( fd ).split()[1]
                    n_ep = self._read_commented_int()
                    n_el = self._read_commented_int()
                    if dim == 2:
                        aux = read_array( fd, n_el, n_ep, nm.int32 )
                        if t_name == 'tri':
                            conns.append( aux )
                            descs.append( '2_3' )
                            is_conn = True
                        else:
                            is_conn = False
                    else:
                        raise NotImplementedError

                    # Skip parameters.
                    n_pv = self._read_commented_int()
                    n_par = self._read_commented_int()
                    for ii in xrange( n_par ):
                        skip_read_line( fd )

                    n_domain = self._read_commented_int()
                    assert_( n_domain == n_el )
                    if is_conn:
                        mat_id = read_array( fd, n_domain, 1, nm.int32 )
                        mat_ids.append( mat_id )
                    else:
                        for ii in xrange( n_domain ):
                            skip_read_line( fd )

                    # Skip up/down pairs.
                    n_ud = self._read_commented_int()
                    for ii in xrange( n_ud ):
                        skip_read_line( fd )
                break

        nod = nm.concatenate( (nod, nm.zeros( (n_nod,1), dtype = nm.int32 ) ),
                              1 )
        nod = nm.ascontiguousarray( nod )

        dim = nod.shape[1] - 1

        conns2 = []
        for ii, conn in enumerate( conns ):
            conns2.append( nm.c_[conn, mat_ids[ii]] )

        conns_in, mat_ids = sort_by_mat_id( conns2 )
        conns, mat_ids, descs = split_by_mat_id( conns_in, mat_ids, descs )
        mesh._set_data( nod, conns, mat_ids, descs )

#        mesh.write( 'pokus.mesh', io = 'auto' )

        self.fd = None
        return mesh

    ##
    # c: 20.03.2008, r: 20.03.2008
    def write( self, filename, mesh, out = None ):
        raise NotImplementedError

##
# c: 23.06.2008
class HDF5MeshIO( MeshIO ):
    format = "hdf5"

    import string
    _all = ''.join( map( chr, range( 256 ) ) )
    _letters = string.letters + string.digits + '_'
    _rubbish = ''.join( [ch for ch in set( _all ) - set( _letters )] )
    _tr = string.maketrans( _rubbish, '_' * len( _rubbish ) )

    def read( self, mesh, **kwargs ):
        fd = pt.openFile( self.filename, mode = "r" )

        mesh_group = fd.root.mesh

        mesh.name = mesh_group.name.read()
        nodes = mesh_group.nod0.read()

        n_gr =  mesh_group.n_gr.read()

        conns = []
        descs = []
        mat_ids = []
        for ig in xrange( n_gr ):
            gr_name = 'group%d' % ig
            group = mesh_group._f_getChild( gr_name )
            conns.append( group.conn.read() )
            mat_ids.append( group.mat_id.read() )
            descs.append( group.desc.read() )

        fd.close()
        mesh._set_data( nodes, conns, mat_ids, descs )

        return mesh

    def write( self, filename, mesh, out = None, ts = None ):
        from time import asctime

        if pt is None:
            output( 'pytables not imported!' )
            raise ValueError

        step = get_default_attr( ts, 'step', 0 )
        if step == 0:
            # A new file.
            fd = pt.openFile( filename, mode = "w",
                              title = "SfePy output file" )

            mesh_group = fd.createGroup( '/', 'mesh', 'mesh' )

            fd.createArray( mesh_group, 'name', mesh.name, 'name' )
            fd.createArray( mesh_group, 'nod0', mesh.nod0, 'vertices' )
            fd.createArray( mesh_group, 'n_gr', len( mesh.conns ), 'n_gr' )
            for ig, conn in enumerate( mesh.conns ):
                conn_group = fd.createGroup( mesh_group, 'group%d' % ig,
                                            'connectivity group' )
                fd.createArray( conn_group, 'conn', conn, 'connectivity' )
                fd.createArray( conn_group, 'mat_id', mesh.mat_ids[ig], 'material id' )
                fd.createArray( conn_group, 'desc', mesh.descs[ig], 'element Type' )

            if ts is not None:
                ts_group = fd.createGroup( '/', 'ts', 'time stepper' )
                fd.createArray( ts_group, 't0', ts.t0, 'initial time' )
                fd.createArray( ts_group, 't1', ts.t1, 'final time'  )
                fd.createArray( ts_group, 'dt', ts.dt, 'time step' )
                fd.createArray( ts_group, 'n_step', ts.n_step, 'n_step' )

            tstat_group = fd.createGroup( '/', 'tstat', 'global time statistics' )
            fd.createArray( tstat_group, 'created', asctime(),
                            'file creation time' )
            fd.createArray( tstat_group, 'finished', '.' * 24,
                            'file closing time' )

            fd.createArray( fd.root, 'last_step', nm.array( [0], dtype = nm.int32 ),
                            'last saved step' )

            fd.close()

        if out is not None:
            if ts is None:
                step, time, nt  = 0, 0.0, 0.0
            else:
                step, time, nt = ts.step, ts.time, ts.nt

            # Existing file.
            fd = pt.openFile( filename, mode = "r+" )

            step_group = fd.createGroup( '/', 'step%d' % step, 'time step data' )
            name_dict = {}
            for key, val in out.iteritems():
    #            print key
                if val.dofs is None:
                    dofs = (-1,)
                else:
                    dofs = val.dofs

                group_name = '_' + key.translate( self._tr )
                data_group = fd.createGroup( step_group, group_name, '%s data' % key )
                fd.createArray( data_group, 'data', val.data, 'data' )
                fd.createArray( data_group, 'mode', val.mode, 'mode' )
                fd.createArray( data_group, 'dofs', dofs, 'dofs' )
                fd.createArray( data_group, 'name', val.name, 'object name' )
                fd.createArray( data_group, 'var_name',
                                val.var_name, 'object parent name' )
                fd.createArray( data_group, 'dname', key, 'data name' )
                name_dict[key] = group_name

            step_group._v_attrs.name_dict = name_dict
            fd.root.last_step[0] = step

            fd.removeNode( fd.root.tstat.finished )
            fd.createArray( fd.root.tstat, 'finished', asctime(),
                            'file closing time' )
            fd.close()

    def read_time_stepper( self, filename = None ):
        filename = get_default( filename, self.filename )
        fd = pt.openFile( filename, mode = "r" )

        ts_group = fd.root.ts
        out =  (ts_group.t0.read(), ts_group.t1.read(),
                ts_group.dt.read(), ts_group.n_step.read())
        fd.close()
        return out

    def _get_step_group( self, step, filename = None ):
        filename = get_default( filename, self.filename )
        fd = pt.openFile( filename, mode = "r" )

        gr_name = 'step%d' % step
        try:
            step_group = fd.getNode( fd.root, gr_name )
        except:
            output( 'step %d data not found - premature end of file?' % step )
            fd.close()
            return None, None

        return fd, step_group

    def read_data( self, step, filename = None ):
        fd, step_group = self._get_step_group( step, filename = filename )
        if fd is None: return None

        out = {}
        for data_group in step_group:
            key = data_group.dname.read()
            out[key] = Struct( name = data_group.name.read(),
                               mode = data_group.mode.read(),
                               data = data_group.data.read(),
                               dofs = tuple( data_group.dofs.read() ) )
            if out[key].dofs == (-1,):
                out[key].dofs = None

        fd.close()

        return out

    def read_data_header( self, dname, step = 0, filename = None ):
        fd, step_group = _get_step_group( step, filename = filename )
        if fd is None: return None

        groups = step_group._v_groups
        for name, data_group in groups.iteritems():
            key = data_group.dname.read()
            if key == dname:
                mode = data_group.mode.read()
                fd.close()
                return mode, name

        fd.close()
        raise KeyError, 'non-existent data: %s' % dname

    def read_time_history( self, node_name, indx, filename = None ):
        filename = get_default( filename, self.filename )
        fd = pt.openFile( filename, mode = "r" )

        th = dict_from_keys_init( indx, list )
        for step in xrange( fd.root.last_step[0] + 1 ):
            gr_name = 'step%d' % step

            step_group = fd.getNode( fd.root, gr_name )
            data = step_group._f_getChild( node_name ).data

            for ii in indx:
                th[ii].append( nm.array( data[ii] ) )

        fd.close()

        for key, val in th.iteritems():
            aux = nm.array( val )
            if aux.ndim == 4: # cell data.
                aux = aux[:,0,:,0]
            th[key] = aux

        return th

    def read_variables_time_history( self, var_names, ts, filename = None ):
        filename = get_default( filename, self.filename )
        fd = pt.openFile( filename, mode = "r" )

        assert_( (fd.root.last_step[0] + 1) == ts.n_step )

        ths = dict_from_keys_init( var_names, list )

        arr = nm.asarray
        for step in xrange( ts.n_step ):
            gr_name = 'step%d' % step
            step_group = fd.getNode( fd.root, gr_name )
            name_dict = step_group._v_attrs.name_dict
            for var_name in var_names:
                data = step_group._f_getChild( name_dict[var_name] ).data
                ths[var_name].append( arr( data.read() ) )

        fd.close()

        return ths

##
# c: 05.02.2008, r: 05.02.2008
var_dict = vars().items()
io_table = {}

for key, var in var_dict:
    try:
        if is_derived_class( var, MeshIO ):
            io_table[var.format] = var
    except TypeError:
        pass
del var_dict

##
# c: 05.02.2008, r: 05.02.2008
def any_from_filename( filename ):
    aux, ext = op.splitext( filename )
    format = supported_formats[ext]
    try:
        return io_table[format]( filename )
    except KeyError:
        output( 'unsupported mesh file suffix: %s' % ext )
        raise

insert_static_method( MeshIO, any_from_filename )
del any_from_filename
