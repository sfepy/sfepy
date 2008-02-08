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
# c: 05.02.2008, r: 05.02.2008
def readList( fd, nItem, dtype ):
    vals = []
    ii = 0
    while ii < nItem:
        line = [dtype( ic ) for ic in fd.readline().split()]
        vals.append( line )
        ii += len( line )
    if ii > nItem:
        output( 'corrupted row?', line, ii, nItem  )
        raise ValueError

    return vals

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
