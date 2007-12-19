import numpy as nm
from numpy import array
import re, sys
import os.path as op
from base import output, Struct, pause, dictFromKeysInit
from time import asctime
try:
    import tables as pt
except:
    pass

##
# 27.04.2006, c
def getTrunk( fileName ):
    return op.splitext( op.basename( fileName ) )[0]

##
# 03.02.2004, c
def readToken( file ):
    
    out = "";
    # Skip initial whitespace.
    while (1):
        ch = file.read( 1 );
        if (ch.isspace()): continue
        elif (len( ch ) == 0): return out
        else: break

    while (not(ch.isspace())):
        out = out + ch;
        ch = file.read( 1 );
        if (len( ch ) == 0): break

    return out
            
        
##
# 03.02.2004, c
def readTuple( file, nItem, nTuple ):

    out = ();
    for it in range( 0, nTuple ):
        token = ();
        for ii in range( 0, nItem ):
            token = token + (readToken( file ),);
#            print token[ii];
        if (len( token[ii] ) == 0):
            output( "Corrupted file (token %d)!" % ii )
            raise "ERR_CorruptedFile"
        out = out + (token,);

    return out;

##
# 12.02.2004, c
def readArray( file, nRow, nCol, dtype ):
    val = [];
    for ir in range( 0, nRow ):
        try:
            while 1:
                line = file.readline();
                if (line[0] == "#"):
                    continue
                else:
                    break
        except:
            output( "Array (%d, %d) reading failed!" % (nRow, nCol) );
            raise
        row = [float( ii ) for ii in line.split()]
#        print ir, row
        val.append( row );

    val = array( val, dtype );
    return val

##
# 17.02.2004, c
# 08.12.2004
# 17.06.2005
# 13.11.2007
def readMesh( file ):
    while 1:
        try:
            line = file.readline()
        except:
            output( 'cannot read mesh header!' )
            raise
        if len( line ) == 1: continue
        if line[0] == '#': continue
        aux = line.split()
        if aux[0] == 'Dimension':
            dim = int( aux[1] )
            break
        
    conns = [];
    desc = [];
    while 1:
        try:
            line = file.readline()
            if (len( line ) == 0): break
            if len( line ) == 1: continue
        except EOFError:
            break
        except:
            output( "Reading " + file.name + " failed!" )
            raise
        if (line[:-1] == 'Vertices'):
            num = int( readToken( file ) );
            nod = readArray( file, num, dim + 1, 'd' );
##                 print nod
        elif (line[:-1] == 'Tetrahedra'):
            num = int( readToken( file ) );
            conns.append( readArray( file, num, 5, 'i' ) );
            desc.append( '3_4' )
        elif (line[:-1] == 'Hexahedra'):
            num = int( readToken( file ) );
            conns.append( readArray( file, num, 9, 'i' ) );
            desc.append( '3_8' )
        elif (line[:-1] == 'Triangles'):
            num = int( readToken( file ) );
            conns.append( readArray( file, num, 4, 'i' ) );
            desc.append( '2_3' )
        elif (line[:-1] == 'Quadrilaterals'):
            num = int( readToken( file ) );
            conns.append( readArray( file, num, 5, 'i' ) );
            desc.append( '2_4' )
        elif line[0] == '#':
            continue
        else:
            output( "corrupted file (line '%s')!" % line )
            raise ValueError
            
##                 print conns
    return nod, conns, desc

##
# 19.01.2005, c
# 03.03.2005
# 17.06.2005
# 03.10.2005
# 04.10.2005
# 05.06.2007
def writeMesh( fd, nod, conns, desc ):
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

##
# 12.10.2005, c
def writeBB( fd, array, dtype ):

    fd.write( '3 %d %d %d\n' % (array.shape[1], array.shape[0], dtype) )
    format = ' '.join( ['%.5e'] * array.shape[1] + ['\n'] )

    for row in array:
        fd.write( format % tuple( row ) )

vtkHeader = r"""# vtk DataFile Version 2.0
generated by %s
ASCII
DATASET UNSTRUCTURED_GRID
"""

vtkCellTypes = {'2_2' : 3, '2_4' : 9, '2_3' : 5,
                '3_2' : 3, '3_4' : 10, '3_8' : 12 }

##
# 15.12.2005, c
# 02.08.2006
# 07.09.2006
# 21.09.2006
# 29.09.2006
# 17.07.2007
def writeVTK( fd, mesh, out = None ):

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
                else:
                    aux = data[:,[0,3,4,6,1,5,7,8,2]]
            else:
                zz = nm.zeros( (data.shape[0], 1), dtype = nm.float64 );
                if nr == sym:
                    aux = nm.c_[data[:,[0,2]], zz, data[:,[2,1]],
                                zz, zz, zz, zz]
                else:
                    aux = nm.c_[data[:,[0,2]], zz, data[:,[3,1]],
                                zz, zz, zz, zz]
            for row in aux:
                fd.write( form % tuple( row ) )
            
        else:
            raise NotImplementedError, (nr, nc)

import string
_all = ''.join( map( chr, range( 256 ) ) )
_letters = string.letters + string.digits + '_'
_rubbish = ''.join( [ch for ch in set( _all ) - set( _letters )] )
_tr = string.maketrans( _rubbish, '_' * len( _rubbish ) )

##
# 26.09.2006, c
# 27.09.2006
# 18.06.2007
def writeHDF5( fileName, mesh = None, ts = None, out = None ):

    if mesh is not None:
        # A new file.
        fd = pt.openFile( fileName, mode = "w",
                          title = "SFE output file" )
        
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
            if val.dofTypes is None:
                dofTypes = (-1,)
            else:
                dofTypes = val.dofTypes

            groupName = '_' + key.translate( _tr )
            dataGroup = fd.createGroup( stepGroup, groupName, '%s data' % key )
            fd.createArray( dataGroup, 'data', val.data, 'data' )
            fd.createArray( dataGroup, 'mode', val.mode, 'mode' )
            fd.createArray( dataGroup, 'dofTypes', dofTypes, 'dofTypes' )
            fd.createArray( dataGroup, 'name', val.name, 'object name' )
            fd.createArray( dataGroup, 'dname', key, 'data name' )
            nameDict[key] = groupName

        stepGroup._v_attrs.nameDict = nameDict
        fd.root.lastStep[0] = step
        
        fd.removeNode( fd.root.tstat.finished )
        fd.createArray( fd.root.tstat, 'finished', asctime(),
                        'file closing time' )
        fd.close()


##
# 26.09.2006, c
def readMeshHDF5( fileName ):
    fd = pt.openFile( fileName, mode = "r" )

    meshGroup = fd.root.mesh

    name = meshGroup.name.read()
    nod0 = meshGroup.nod0.read()
    
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
    return name, nod0, conns, matIds, descs

##
# 26.09.2006, c
def readTimeStepperHDF5( fileName ):
    fd = pt.openFile( fileName, mode = "r" )

    tsGroup = fd.root.ts
    out =  (tsGroup.t0.read(), tsGroup.t1.read(),
            tsGroup.dt.read(), tsGroup.nStep.read())
    fd.close()
    return out 

##
# 26.09.2006, c
def _getStepGroup( fileName, step ):
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
# 26.09.2006, c
def readDataHDF5( fileName, step ):
    fd, stepGroup = _getStepGroup( fileName, step )
    if fd is None: return None
    
    out = {}
    for dataGroup in stepGroup:
        key = dataGroup.dname.read()
        out[key] = Struct( name = dataGroup.name.read(),
                           mode = dataGroup.mode.read(),
                           data = dataGroup.data.read(),
                           dofTypes = tuple( dataGroup.dofTypes.read() ) )
        if out[key].dofTypes == (-1,):
            out[key].dofTypes = None
            
    fd.close()

    return out

##
# 26.09.2006, c
def readDataHeaderHDF5( fileName, dname, step = 0 ):
    fd, stepGroup = _getStepGroup( fileName, step )
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
# 27.09.2006, c
def readTimeHistoryHDF5( fileName, nodeName, indx ):
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
# 14.06.2007, c
# 15.06.2007
# 18.06.2007
def readVariablesTimeHistoryHDF5( fileName, varNames, ts ):
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
# 27.09.2006, c
def writeDictHDF5( fileName, adict, level = 0, group = None, fd = None ):

    if level == 0:
        fd = pt.openFile( fileName, mode = "w",
                          title = "Recursive dict dump" )
        group = '/'

    for key, val in adict.iteritems():
#        print level * ' ', key, '->', group
        if isinstance( val, dict ):
            group2 = fd.createGroup( group, '_' + str( key ), '%s group' % key )
            writeDictHDF5( fileName, val, level + 1, group2, fd )
        else:
            fd.createArray( group, '_' + str( key ), val, '%s data' % key )
            
    if level == 0:
        fd.close()

##
# 09.07.2007, c
def readDictHDF5( fileName, level = 0, group = None, fd = None ):
    out = {}

    if level == 0:
        fd = pt.openFile( fileName, mode = "r" )
        group = fd.root

    for name, gr in group._v_groups.iteritems():
#        print level * ' ', gr, '->', group
        name = name.replace( '_', '' )
        out[name] = readDictHDF5( fileName, level + 1, gr, fd )

    for name, data in group._v_leaves.iteritems():
        name = name.replace( '_', '' )
        out[name] = data.read()

    if level == 0:
        fd.close()
        
    return out

##
# 02.07.2007, c
def writeSparseMatrixHDF5( fileName, mtx, name = 'a sparse matrix' ):
    """Assume CSR/CSC."""
    fd = pt.openFile( fileName, mode = "w", title = name )
    try:
        info = fd.createGroup( '/', 'info' )
        fd.createArray( info, 'dtype', mtx.dtype.str )
        fd.createArray( info, 'shape', mtx.shape )
        fd.createArray( info, 'format', mtx.format )

        data = fd.createGroup( '/', 'data' )
        fd.createArray( data, 'data', mtx.data )
        fd.createArray( data, 'indptr', mtx.indptr )
        fd.createArray( data, 'indices', mtx.indices )

    except:
        print 'matrix must be in SciPy sparse CSR/CSC format!'
        print mtx.__repr__()
        raise

    fd.close()

##
# 02.07.2007, c
# 08.10.2007
def readSparseMatrixHDF5( fileName, outputFormat = None ):
    import scipy.sparse as sp
    constructors = {'csr' : sp.csr_matrix, 'csc' : sp.csc_matrix}
    
    fd = pt.openFile( fileName, mode = "r" )
    info = fd.root.info
    data = fd.root.data

    format = info.format.read()
    if not isinstance( format, str ):
        format = format[0]

    dtype = info.dtype.read()
    if not isinstance( dtype, str ):
        dtype = dtype[0]

    if outputFormat is None:
        constructor = constructors[format]
    else:
        constructor = constructors[outputFormat]

    if format in ['csc', 'csr']:
        mtx = constructor( (data.data.read(),
                            data.indices.read(), data.indptr.read()),
                           dims = info.shape.read(), dtype = dtype )
    elif format == 'coo':
        mtx = constructor( (data.data.read(),
                            nm.c_[data.rows.read(), data.cols.read()].T),
                           dims = info.shape.read(), dtype = dtype )
    else:
        print format
        raise ValueError
    fd.close()

    if outputFormat in ['csc', 'csr']:
        mtx.sort_indices()
    
    return mtx
