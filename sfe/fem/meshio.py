from sfe.base.base import *
from sfe.base.ioutils import readToken, readArray, readList, vtkInverseCellTypes

import os.path as op

supportedFormats = {
    '.mesh' : 'medit',
    '.vtk'  : 'vtk',
}

##
# c: 05.02.2008
class MeshIO( Struct ):
    format = None
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def __init__( self, fileName, **kwargs ):
        Struct.__init__( self, fileName = fileName, **kwargs )

    ##
    # c: 05.02.2008, r: 05.02.2008
    def read( self, **kwargs ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

    ##
    # c: 05.02.2008, r: 05.02.2008
    def write( self, **kwargs ):
        print 'called an abstract MeshIO instance!'
        raise ValueError

##
# c: 05.02.2008
class MeditMeshIO( MeshIO ):
    format = 'medit'
    
    ##
    # c: 17.02.2004, r: 05.02.2008
    def read( self, **kwargs ):
        fd = open( self.fileName, 'r' )
        while 1:
            try:
                line = fd.readline()
            except:
                output( 'cannot read mesh header!' )
                raise
            if len( line ) == 1: continue
            if line[0] == '#': continue
            aux = line.split()
            if aux[0] == 'Dimension':
                dim = int( aux[1] )
                break

        conns = []
        desc = []
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
                conns.append( readArray( fd, num, 5, nm.int32 ) )
                conns[-1][:,:-1] -= 1
                desc.append( '3_4' )
            elif (line[:-1] == 'Hexahedra'):
                num = int( readToken( fd ) )
                conns.append( readArray( fd, num, 9, nm.int32 ) )
                conns[-1][:,:-1] -= 1
                desc.append( '3_8' )
            elif (line[:-1] == 'Triangles'):
                num = int( readToken( fd ) )
                conns.append( readArray( fd, num, 4, nm.int32 ) )
                conns[-1][:,:-1] -= 1
                desc.append( '2_3' )
            elif (line[:-1] == 'Quadrilaterals'):
                num = int( readToken( fd ) )
                conns.append( readArray( fd, num, 5, nm.int32 ) )
                conns[-1][:,:-1] -= 1
                desc.append( '2_4' )
            elif line[0] == '#':
                continue
            else:
                output( "corrupted file (line '%s')!" % line )
                raise ValueError
        fd.close()
        return nod, conns, desc
        
##
# c: 05.02.2008
class VTKMeshIO( MeshIO ):
    format = 'vtk'
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def read( self, **kwargs ):
        fd = open( self.fileName, 'r' )
        mode = 'header'
        modeStatus = 0
        nod = conns = desc = None
        while 1:
            try:
                line = fd.readline()
##                print line
            except EOFError:
                print 'sdsd'
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
                        matId = readList( fd, nEl, nm.int32 )
                        modeStatus = 0
                        mode = 'finished'
            elif mode == 'finished':
                break
        fd.close()

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

        return nod, conns, desc

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
