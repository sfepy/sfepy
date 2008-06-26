import numpy as nm
from numpy import array
import re, sys
import os.path as op
from base import output, Struct, pause, dictFromKeysInit
try:
    import tables as pt
except:
    pt = None

##
# 27.04.2006, c
def getTrunk( fileName ):
    return op.splitext( op.basename( fileName ) )[0]

##
# c: 20.03.2008, r: 20.03.2008
def skipReadLine( fd ):
    while 1:
        try:
            line = fd.readline().strip()
        except EOFError:
            break
        except:
            output( "reading " + fd.name + " failed!" )
            raise
        if (len( line ) == 0) or (line[0] == '#'): continue
        return line

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
# c: 12.02.2004, r: 20.03.2008
def readArray( file, nRow, nCol, dtype ):
    """nCol is basically ignored, nCol == -1 -> intentionally ignored."""
    val = []
    for ir in range( 0, nRow ):
        try:
            while 1:
                line = file.readline().split()
                if (len( line ) == 0) or (line[0] == "#"):
                    continue
                else:
                    break
        except:
            output( "Array (%d, %d) reading failed!" % (nRow, nCol) )
            raise
        row = [float( ii ) for ii in line]
        val.append( row )

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
