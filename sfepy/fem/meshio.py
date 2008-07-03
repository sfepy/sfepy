from sfepy.base.base import *
from sfepy.base.ioutils import skipReadLine, readToken, readArray, readList, pt
import sfepy.base.la as la
from sfepy.base.progressbar import MyBar
import os.path as op

supportedFormats = {
    '.mesh' : 'medit',
    '.vtk'  : 'vtk',
    '.node' : 'tetgen',
    '.txt'  : 'comsol',
    '.h5'   : 'hdf5',
}

##
# c: 15.02.2008, r: 15.02.2008
def sortByMatID( connsIn ):

    # Sort by matId within a group, preserve order.
    conns = []
    matIds = []
    for ig, conn in enumerate( connsIn ):
        ii = nm.argsort( conn[:,-1], kind = 'mergesort' )
        conn = conn[ii]

        conns.append( conn[:,:-1].copy() )
        matIds.append( conn[:,-1].copy() )
    return conns, matIds

##
# connsIn must be sorted by matId within a group!
# c: 16.06.2005, r: 15.02.2008
def splitByMatId( connsIn, matIdsIn, descsIn ):

    conns = []
    matIds = []
    descs = []

    for ig, conn in enumerate( connsIn ):
        one = nm.array( [-1], nm.int32 )
        ii = la.diff( nm.concatenate( (one, matIdsIn[ig], one) ) ).nonzero()[0]
        nGr = len( ii ) - 1;
#        print ii, nGr
        for igr in range( 0, nGr ):
            conns.append( conn[ii[igr]:ii[igr+1],:].copy() )
            matIds.append( matIdsIn[ig][ii[igr]:ii[igr+1]] )
            descs.append( descsIn[ig] )

    return (conns, matIds, descs)


##
# 12.10.2005, c
def writeBB( fd, array, dtype ):

    fd.write( '3 %d %d %d\n' % (array.shape[1], array.shape[0], dtype) )
    format = ' '.join( ['%.5e'] * array.shape[1] + ['\n'] )

    for row in array:
        fd.write( format % tuple( row ) )

##
# c: 03.10.2005, r: 08.02.2008
def joinConnGroups( conns, descs, matIds, concat = False ):
    """Join groups of the same element type."""

    el = dictFromKeysInit( descs, list )
    for ig, desc in enumerate( descs ):
        el[desc].append( ig )
    groups = [ii for ii in el.values() if ii]
##     print el, groups

    descsOut, connsOut, matIdsOut = [], [], []
    for group in groups:
        nEP = conns[group[0]].shape[1]

        conn = nm.zeros( (0, nEP), nm.int32 )
        matId = nm.zeros( (0,), nm.int32 )
        for ig in group:
            conn = nm.concatenate( (conn, conns[ig]) )
            matId = nm.concatenate( (matId, matIds[ig]) )

        if concat:
            conn = nm.concatenate( (conn, matId[:,nm.newaxis]), 1 )
        else:
            matIdsOut.append( matId )
        connsOut.append( conn )
        descsOut.append( descs[group[0]] )

    if concat:
        return connsOut, descsOut
    else:
        return connsOut, descsOut, matIdsOut

##
# c: 05.02.2008
class MeshIO( Struct ):
    format = None

    ##
    # c: 05.02.2008, r: 05.02.2008
    def __init__( self, fileName, **kwargs ):
        Struct.__init__( self, fileName = fileName, **kwargs )

    ##
    # c: 03.07.2008, r: 03.07.2008
    def read_dimension( self, ret_fd = False ):
        print 'called an abstract MeshIO instance!'
        raise ValueError
        
    ##
    # c: 05.02.2008, r: 26.03.2008
    def read( self, mesh, *args, **kwargs ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

    ##
    # c: 05.02.2008, r: 26.03.2008
    def write( self, fileName, mesh, *args, **kwargs ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

##
# c: 05.02.2008
class MeditMeshIO( MeshIO ):
    format = 'medit'

    ##
    # c: 03.07.2008, r: 03.07.2008
    def read_dimension( self, ret_fd = False ):
        fd = open( self.fileName, 'r' )
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
                dim = int( aux[1] )
                break

        if ret_fd:
            return dim, fd
        else:
            fd.close()
            return dim

    ##
    # c: 17.02.2004, r: 03.07.2008
    def read( self, mesh, **kwargs ):
        dim, fd  = self.read_dimension( ret_fd = True )

        connsIn = []
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
            if (line[:-1] == 'Vertices'):
                num = int( readToken( fd ) )
                nod = readArray( fd, num, dim + 1, nm.float64 )
    ##                 print nod
            elif (line[:-1] == 'Tetrahedra'):
                num = int( readToken( fd ) )
                connsIn.append( readArray( fd, num, 5, nm.int32 ) )
                connsIn[-1][:,:-1] -= 1
                descs.append( '3_4' )
            elif (line[:-1] == 'Hexahedra'):
                num = int( readToken( fd ) )
                connsIn.append( readArray( fd, num, 9, nm.int32 ) )
                connsIn[-1][:,:-1] -= 1
                descs.append( '3_8' )
            elif (line[:-1] == 'Triangles'):
                num = int( readToken( fd ) )
                connsIn.append( readArray( fd, num, 4, nm.int32 ) )
                connsIn[-1][:,:-1] -= 1
                descs.append( '2_3' )
            elif (line[:-1] == 'Quadrilaterals'):
                num = int( readToken( fd ) )
                connsIn.append( readArray( fd, num, 5, nm.int32 ) )
                connsIn[-1][:,:-1] -= 1
                descs.append( '2_4' )
            elif line[0] == '#':
                continue
            else:
                output( "corrupted file (line '%s')!" % line )
                raise ValueError
        fd.close()

        connsIn, matIds = sortByMatID( connsIn )

        # Detect wedges and pyramides -> separate groups.
        if ('3_8' in descs):
            ic = descs.index( '3_8' )

            connIn = connsIn.pop( ic )
            flag = nm.zeros( (connIn.shape[0],), nm.int32 )
            for ii, el in enumerate( connIn ):
                if (el[4] == el[5]):
                    if (el[5] == el[6]):
                        flag[ii] = 2
                    else:
                        flag[ii] = 1

            conn = []
            desc = []

            ib = nm.where( flag == 0 )[0]
            if (len( ib ) > 0):
                conn.append( connIn[ib] )
                desc.append( '3_8' )

            iw = nm.where( flag == 1 )[0]
            if (len( iw ) > 0):
                ar = nm.array( [0,1,2,3,4,6,8], nm.int32 )
                conn.append( la.rect( connIn, iw, ar ) )
                desc.append( '3_6' )

            ip = nm.where( flag == 2 )[0]
            if (len( ip ) > 0):
                ar = nm.array( [0,1,2,3,4,8], nm.int32 )
                conn.append( la.rect( connIn, ip, ar ) )
                desc.append( '3_5' )

##             print "brick split:", ic, ":", ib, iw, ip, desc

            connsIn[ic:ic] = conn
            del( descs[ic] )
            descs[ic:ic] = desc

        conns, matIds, descs = splitByMatId( connsIn, matIds, descs )
        mesh._setData( nod, conns, matIds, descs )

        return mesh

    ##
    # c: 19.01.2005, r: 08.02.2008
    def write( self, fileName, mesh, out = None ):
        fd = open( fileName, 'w' )

        nod = mesh.nod0
        conns, desc = joinConnGroups( mesh.conns, mesh.descs,
                                      mesh.matIds, concat = True )

        nNod, dim = nod.shape
        dim -= 1

        fd.write( "MeshVersionFormatted 1\nDimension %d\n" % dim )

        fd.write( "Vertices\n%d\n" % nNod )
        if (dim == 2):
            for ii in range( nNod ):
                nn = nod[ii]
                fd.write( "%.8e %.8e %d\n" % (nn[0], nn[1], nn[2]) )
        else:
            for ii in range( nNod ):
                nn = nod[ii]
                fd.write( "%.8e %.8e %.8e %d\n" % (nn[0], nn[1], nn[2], nn[3]) )

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


vtkHeader = r"""# vtk DataFile Version 2.0
generated by %s
ASCII
DATASET UNSTRUCTURED_GRID
"""
vtkCellTypes = {'2_2' : 3, '2_4' : 9, '2_3' : 5,
                '3_2' : 3, '3_4' : 10, '3_8' : 12 }
vtkInverseCellTypes = {(3, 2) : '2_2', (9, 2) : '2_4', (5, 2) : '2_3',
                       (3, 3) : '3_2', (10, 3) : '3_4', (12, 3) : '3_8' }

##
# c: 05.02.2008
class VTKMeshIO( MeshIO ):
    format = 'vtk'

    ##
    # c: 03.07.2008, r: 03.07.2008
    def read_dimension( self, ret_fd = False ):
        fd = open( self.fileName, 'r' )
        while 1:
            try:
                line = fd.readline().split()
                if line[0] == 'POINTS':
                    nod = readArray( fd, 1, -1, nm.float64 )
                    dim = nod.shape[1]
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
    # c: 05.02.2008, r: 26.06.2008
    def read( self, mesh, **kwargs ):
        fd = open( self.fileName, 'r' )
        mode = 'header'
        modeStatus = 0
        nod = conns = desc = matId = None
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
                if modeStatus == 0:
                    if line.strip() == 'ASCII':
                        modeStatus = 1
                elif modeStatus == 1:
                    if line.strip() == 'DATASET UNSTRUCTURED_GRID':
                        modeStatus = 0
                        mode = 'points'

            elif mode == 'points':
                line = line.split()
                if line[0] == 'POINTS':
                    nNod = int( line[1] )
                    nod = readArray( fd, nNod, -1, nm.float64 )
                    mode = 'cells'

            elif mode == 'cells':
                line = line.split()
                if line[0] == 'CELLS':
                    nEl, nVal = map( int, line[1:3] )
                    rawConn = readList( fd, nVal, int )
                    mode = 'cell_types'

            elif mode == 'cell_types':
                line = line.split()
                if line[0] == 'CELL_TYPES':
                    assert int( line[1] ) == nEl
                    cellTypes = readArray( fd, nEl, -1, nm.int32 )
                    mode = 'matId'

            elif mode == 'matId':
                if modeStatus == 0:
                    line = line.split()
                    if line[0] == 'CELL_DATA':
                        assert int( line[1] ) == nEl
                        modeStatus = 1
                elif modeStatus == 1:
                    if line.strip() == 'SCALARS matId float 1':
                        modeStatus = 2
                elif modeStatus == 2:
                    if line.strip() == 'LOOKUP_TABLE default':
                        matId = readList( fd, nEl, int )
                        modeStatus = 0
                        mode = 'finished'
            elif mode == 'finished':
                break
        fd.close()

        if matId is None:
            matId = [[0]] * nEl

        nod = nm.concatenate( (nod, nm.zeros( (nNod,1), dtype = nm.int32 ) ),
                              1 )
        nod = nm.ascontiguousarray( nod )

        dim = nod.shape[1] - 1
        cellTypes = cellTypes.squeeze()

        dconns = {}
        for iel, row in enumerate( rawConn ):
            ct = vtkInverseCellTypes[(cellTypes[iel],dim)]
            dconns.setdefault( ct, [] ).append( row[1:] + matId[iel] )

        desc = []
        conns = []
        for key, conn in dconns.iteritems():
            desc.append( key )
            conns.append( nm.array( conn, dtype = nm.int32 ) )

        connsIn, matIds = sortByMatID( conns )
        conns, matIds, descs = splitByMatId( connsIn, matIds, desc )
        mesh._setData( nod, conns, matIds, descs )

        return mesh

    ##
    # c: 15.12.2005, r: 07.05.2008
    def write( self, fileName, mesh, out = None ):

        fd = open( fileName, 'w' )
        fd.write( vtkHeader % op.basename( sys.argv[0] ) )

        nNod, nc = mesh.nod0.shape
        dim = nc - 1
        sym = dim * (dim + 1) / 2

        fd.write( '\nPOINTS %d float\n' % nNod )

        aux = mesh.nod0[:,:dim]
        if dim == 2:
            aux = nm.hstack( (aux, nm.zeros( (nNod, 1), dtype = nm.float64 ) ) )
        for row in aux:
            fd.write( '%e %e %e\n' % tuple( row ) )

        nEl, nEls, nEPs = mesh.nEl, mesh.nEls, mesh.nEPs
        totalSize = nm.dot( nEls, nEPs + 1 )
        fd.write( '\nCELLS %d %d\n' % (nEl, totalSize) )

        ct = []
        for ig, conn in enumerate( mesh.conns ):
            nn = nEPs[ig] + 1
            ct += [vtkCellTypes[mesh.descs[ig]]] * nEls[ig]
            format = ' '.join( ['%d'] * nn + ['\n'] )

            for row in conn:
                fd.write( format % ((nn-1,) + tuple( row )) )

        fd.write( '\nCELL_TYPES %d\n' % nEl )
        fd.write( ''.join( ['%d\n' % ii for ii in ct] ) )

        if out is None: return

        pointKeys = [key for key, val in out.iteritems() if val.mode == 'vertex']
        if len( pointKeys ):
            fd.write( '\nPOINT_DATA %d\n' % nNod );
        for key in pointKeys:
            val = out[key]
            nr, nc = val.data.shape
            if nc == 1:
                fd.write( '\nSCALARS %s float %d\n' % (key, nc) )
                fd.write( 'LOOKUP_TABLE default\n' )
                for row in val.data:
                    fd.write( '%e\n' % row )
            elif nc == dim:
                fd.write( '\nVECTORS %s float\n' % key )
                if dim == 2:
                    aux = nm.hstack( (val.data,
                                      nm.zeros( (nr, 1), dtype = nm.float64 ) ) )
                else:
                    aux = val.data
                for row in aux:
                    fd.write( '%e %e %e\n' % tuple( row ) )
            else:
                raise NotImplementedError, nc

        cellKeys = [key for key, val in out.iteritems() if val.mode == 'cell']
        if len( cellKeys ):
            fd.write( '\nCELL_DATA %d\n' % nEl );
        for key in cellKeys:
            val = out[key]
            ne, aux, nr, nc = val.data.shape
            if (nr == 1) and (nc == 1):
                fd.write( '\nSCALARS %s float %d\n' % (key, nc) )
                fd.write( 'LOOKUP_TABLE default\n' )
                for row in val.data.squeeze():
                    fd.write( '%e\n' % row )
            elif (nr == dim) and (nc == 1):
                fd.write( '\nVECTORS %s float\n' % key )
                if dim == 2:
                    aux = nm.hstack( (val.data.squeeze(),
                                      nm.zeros( (ne, 1), dtype = nm.float64 ) ) )
                else:
                    aux = val.data
                for row in aux:
                    fd.write( '%e %e %e\n' % tuple( row.squeeze() ) )
            elif (((nr == sym) or (nr == (dim * dim))) and (nc == 1)) \
                     or ((nr == dim) and (nc == dim)):
                # Below not tested!!!
                fd.write( '\nTENSORS %s float\n' % key );
                form = '%e %e %e\n%e %e %e\n%e %e %e\n\n';
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
                for row in aux:
                    fd.write( form % tuple( row ) )

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
        fname = os.path.splitext(self.fileName)[0]
        nodes=self.getnodes(fname+".node", MyBar("       nodes:"))
        elements, regions = self.getele(fname+".ele", MyBar("       elements:"))
        descs = []
        conns = []
        matIds = []
        nodes = nm.c_[(nm.array(nodes, dtype = nm.float64),
                       nm.zeros(len(nodes), dtype = nm.float64))].copy()
        elements = nm.array( elements, dtype = nm.int32 )-1
        for key, value in regions.iteritems():
            descs.append( "3_4" )
            matIds.append( nm.ones_like(value) * key )
            conns.append( elements[nm.array(value)-1].copy() )

        mesh._setData( nodes, conns, matIds, descs )
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
        assert dim==3
        if verbose: up.init(npoints)
        nodes=[]
        for line in f:
            if line[0]=="#": continue
            l=[float(x) for x in line.split()]
            assert int(l[0])==len(nodes)+1
            l = l[1:]
            nodes.append(tuple(l))
            if verbose: up.update(len(nodes))
        assert npoints==len(nodes)
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
        assert nnod in [4,10]
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
            assert len(l)-2 == 4
            els.append((l[1],l[2],l[3],l[4]))
            regionnum=l[5]
            if regionnum==0:
                print "see %s, element # %d"%(fele,l[0])
                raise "there are elements not belonging to any physical entity"
            if regions.has_key(regionnum):
                regions[regionnum].append(l[0])
            else:
                regions[regionnum]=[l[0]]
            assert l[0]==len(els)
            if verbose: up.update(l[0])
        return els, regions

    ##
    # c: 26.03.2008, r: 26.03.2008
    def write( self, fileName, mesh, out = None ):
        raise NotImplementedError

##
# c: 20.03.2008
class ComsolMeshIO( MeshIO ):
    format = 'comsol'

    ##
    # c: 20.03.2008, r: 20.03.2008
    def _readCommentedInt( self ):
        return int( skipReadLine( self.fd ).split( '#' )[0] )

    ##
    # c: 20.03.2008, r: 20.03.2008
    def read( self, mesh, **kwargs ):

        self.fd = fd = open( self.fileName, 'r' )
        mode = 'header'

        nod = conns = desc = None
        while 1:
            if mode == 'header':
                line = skipReadLine( fd )

                nTags = self._readCommentedInt()
                for ii in xrange( nTags ):
                    skipReadLine( fd )
                nTypes = self._readCommentedInt()
                for ii in xrange( nTypes ):
                    skipReadLine( fd )

                skipReadLine( fd )
                assert skipReadLine( fd ).split()[1] == 'Mesh'
                skipReadLine( fd )
                dim = self._readCommentedInt()
                assert (dim == 2) or (dim == 3)
                nNod = self._readCommentedInt()
                i0 = self._readCommentedInt()
                mode = 'points'

            elif mode == 'points':
                nod = readArray( fd, nNod, dim, nm.float64 )
                mode = 'cells'

            elif mode == 'cells':

                nTypes = self._readCommentedInt()
                conns = []
                descs = []
                matIds = []
                for it in xrange( nTypes ):
                    tName = skipReadLine( fd ).split()[1]
                    nEP = self._readCommentedInt()
                    nEl = self._readCommentedInt()
                    if dim == 2:
                        aux = readArray( fd, nEl, nEP, nm.int32 )
                        if tName == 'tri':
                            conns.append( aux )
                            descs.append( '2_3' )
                            isConn = True
                        else:
                            isConn = False
                    else:
                        raise NotImplementedError

                    # Skip parameters.
                    nPV = self._readCommentedInt()
                    nPar = self._readCommentedInt()
                    for ii in xrange( nPar ):
                        skipReadLine( fd )

                    nDomain = self._readCommentedInt()
                    assert nDomain == nEl
                    if isConn:
                        matId = readArray( fd, nDomain, 1, nm.int32 )
                        matIds.append( matId )
                    else:
                        for ii in xrange( nDomain ):
                            skipReadLine( fd )

                    # Skip up/down pairs.
                    nUD = self._readCommentedInt()
                    for ii in xrange( nUD ):
                        skipReadLine( fd )
                break

        nod = nm.concatenate( (nod, nm.zeros( (nNod,1), dtype = nm.int32 ) ),
                              1 )
        nod = nm.ascontiguousarray( nod )

        dim = nod.shape[1] - 1

        conns2 = []
        for ii, conn in enumerate( conns ):
            conns2.append( nm.c_[conn, matIds[ii]] )

        connsIn, matIds = sortByMatID( conns2 )
        conns, matIds, descs = splitByMatId( connsIn, matIds, descs )
        mesh._setData( nod, conns, matIds, descs )

#        mesh.write( 'pokus.mesh', io = 'auto' )

        self.fd = None
        return mesh

    ##
    # c: 20.03.2008, r: 20.03.2008
    def write( self, fileName, mesh, out = None ):
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

    ##
    # c: 26.09.2006, r: 23.06.2008
    def read( self, mesh, **kwargs ):
        fd = pt.openFile( self.fileName, mode = "r" )

        meshGroup = fd.root.mesh

        mesh.name = meshGroup.name.read()
        nodes = meshGroup.nod0.read()

        nGr =  meshGroup.nGr.read()

        conns = []
        descs = []
        matIds = []
        for ig in xrange( nGr ):
            grName = 'group%d' % ig
            group = meshGroup._f_getChild( grName )
            conns.append( group.conn.read() )
            matIds.append( group.matId.read() )
            descs.append( group.desc.read() )

        fd.close()
        mesh._setData( nodes, conns, matIds, descs )

        return mesh

    ##
    # c: 26.09.2006, r: 23.06.2008
    def write( self, fileName, mesh, out = None, ts = None ):
        from time import asctime

        if pt is None:
            output( 'pytables not imported!' )
            raise ValueError

        step = getDefaultAttr( ts, 'step', 0 )
        if step == 0:
            # A new file.
            fd = pt.openFile( fileName, mode = "w",
                              title = "SfePy output file" )

            meshGroup = fd.createGroup( '/', 'mesh', 'mesh' )

            fd.createArray( meshGroup, 'name', mesh.name, 'name' )
            fd.createArray( meshGroup, 'nod0', mesh.nod0, 'vertices' )
            fd.createArray( meshGroup, 'nGr', len( mesh.conns ), 'nGr' )
            for ig, conn in enumerate( mesh.conns ):
                connGroup = fd.createGroup( meshGroup, 'group%d' % ig,
                                            'connectivity group' )
                fd.createArray( connGroup, 'conn', conn, 'connectivity' )
                fd.createArray( connGroup, 'matId', mesh.matIds[ig], 'material id' )
                fd.createArray( connGroup, 'desc', mesh.descs[ig], 'element Type' )

            if ts is not None:
                tsGroup = fd.createGroup( '/', 'ts', 'time stepper' )
                fd.createArray( tsGroup, 't0', ts.t0, 'initial time' )
                fd.createArray( tsGroup, 't1', ts.t1, 'final time'  )
                fd.createArray( tsGroup, 'dt', ts.dt, 'time step' )
                fd.createArray( tsGroup, 'nStep', ts.nStep, 'nStep' )

            tstatGroup = fd.createGroup( '/', 'tstat', 'global time statistics' )
            fd.createArray( tstatGroup, 'created', asctime(),
                            'file creation time' )
            fd.createArray( tstatGroup, 'finished', '.' * 24,
                            'file closing time' )

            fd.createArray( fd.root, 'lastStep', nm.array( [0], dtype = nm.int32 ),
                            'last saved step' )

            fd.close()

        if out is not None:
            if ts is None:
                step, time, nt  = 0, 0.0, 0.0
            else:
                step, time, nt = ts.step, ts.time, ts.nt

            # Existing file.
            fd = pt.openFile( fileName, mode = "r+" )

            stepGroup = fd.createGroup( '/', 'step%d' % step, 'time step data' )
            nameDict = {}
            for key, val in out.iteritems():
    #            print key
                if val.dofs is None:
                    dofs = (-1,)
                else:
                    dofs = val.dofs

                groupName = '_' + key.translate( self._tr )
                dataGroup = fd.createGroup( stepGroup, groupName, '%s data' % key )
                fd.createArray( dataGroup, 'data', val.data, 'data' )
                fd.createArray( dataGroup, 'mode', val.mode, 'mode' )
                fd.createArray( dataGroup, 'dofs', dofs, 'dofs' )
                fd.createArray( dataGroup, 'name', val.name, 'object name' )
                fd.createArray( dataGroup, 'varName',
                                val.varName, 'object parent name' )
                fd.createArray( dataGroup, 'dname', key, 'data name' )
                nameDict[key] = groupName

            stepGroup._v_attrs.nameDict = nameDict
            fd.root.lastStep[0] = step

            fd.removeNode( fd.root.tstat.finished )
            fd.createArray( fd.root.tstat, 'finished', asctime(),
                            'file closing time' )
            fd.close()

    ##
    # c: 26.09.2006, r: 23.06.2008
    def readTimeStepper( self, fileName = None ):
        fileName = getDefault( fileName, self.fileName )
        fd = pt.openFile( fileName, mode = "r" )

        tsGroup = fd.root.ts
        out =  (tsGroup.t0.read(), tsGroup.t1.read(),
                tsGroup.dt.read(), tsGroup.nStep.read())
        fd.close()
        return out

    ##
    # c: 26.09.2006, r: 23.06.2008
    def _getStepGroup( self, step, fileName = None ):
        fileName = getDefault( fileName, self.fileName )
        fd = pt.openFile( fileName, mode = "r" )

        grName = 'step%d' % step
        try:
            stepGroup = fd.getNode( fd.root, grName )
        except:
            output( 'step %d data not found - premature end of file?' % step )
            fd.close()
            return None, None

        return fd, stepGroup

    ##
    # c: 26.09.2006, r: 23.06.2008
    def readData( self, step, fileName = None ):
        fd, stepGroup = self._getStepGroup( step, fileName = fileName )
        if fd is None: return None

        out = {}
        for dataGroup in stepGroup:
            key = dataGroup.dname.read()
            out[key] = Struct( name = dataGroup.name.read(),
                               mode = dataGroup.mode.read(),
                               data = dataGroup.data.read(),
                               dofs = tuple( dataGroup.dofs.read() ) )
            if out[key].dofs == (-1,):
                out[key].dofs = None

        fd.close()

        return out

    ##
    # c: 26.09.2006, r: 23.06.2008
    def readDataHeader( self, dname, step = 0, fileName = None ):
        fd, stepGroup = _getStepGroup( step, fileName = fileName )
        if fd is None: return None

        groups = stepGroup._v_groups
        for name, dataGroup in groups.iteritems():
            key = dataGroup.dname.read()
            if key == dname:
                mode = dataGroup.mode.read()
                fd.close()
                return mode, name

        fd.close()
        raise KeyError, 'non-existent data: %s' % dname

    ##
    # c: 27.09.2006, r: 23.06.2008
    def readTimeHistory( self, nodeName, indx, fileName = None ):
        fileName = getDefault( fileName, self.fileName )
        fd = pt.openFile( fileName, mode = "r" )

        th = dictFromKeysInit( indx, list )
        for step in xrange( fd.root.lastStep[0] + 1 ):
            grName = 'step%d' % step

            stepGroup = fd.getNode( fd.root, grName )
            data = stepGroup._f_getChild( nodeName ).data

            for ii in indx:
                th[ii].append( nm.array( data[ii] ) )

        fd.close()

        for key, val in th.iteritems():
            aux = nm.array( val )
            if aux.ndim == 4: # cell data.
                aux = aux[:,0,:,0]
            th[key] = aux

        return th

    ##
    # c: 14.06.2007, r: 23.06.2008
    def readVariablesTimeHistory( self, varNames, ts, fileName = None ):
        fileName = getDefault( fileName, self.fileName )
        fd = pt.openFile( fileName, mode = "r" )

        assert (fd.root.lastStep[0] + 1) == ts.nStep

        ths = dictFromKeysInit( varNames, list )

        arr = nm.asarray
        for step in xrange( ts.nStep ):
            grName = 'step%d' % step
            stepGroup = fd.getNode( fd.root, grName )

            nameDict = stepGroup._v_attrs.nameDict
            for varName in varNames:
                data = stepGroup._f_getChild( nameDict[varName] ).data
                ths[varName].append( arr( data ) )

        fd.close()

        return ths

##
# c: 05.02.2008, r: 05.02.2008
varDict = vars().items()
ioTable = {}

for key, var in varDict:
    try:
        if isDerivedClass( var, MeshIO ):
            ioTable[var.format] = var
    except TypeError:
        pass
del varDict

##
# c: 05.02.2008, r: 05.02.2008
def anyFromFileName( fileName ):
    aux, ext = op.splitext( fileName )
    format = supportedFormats[ext]
    try:
        return ioTable[format]( fileName )
    except KeyError:
        output( 'unsupported mesh file suffix: %s' % ext )
        raise

insertStaticMethod( MeshIO, anyFromFileName )
del anyFromFileName
