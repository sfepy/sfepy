import numpy as nm
from numpy import array
import re, sys
import os.path as op
from base import output, Struct, pause, dict_from_keys_init
try:
    import tables as pt
except:
    pt = None

##
# 27.04.2006, c
def get_trunk( file_name ):
    return op.splitext( op.basename( file_name ) )[0]

##
# c: 20.03.2008, r: 20.03.2008
def skip_read_line( fd ):
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
def read_token( file ):
    
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
def read_tuple( file, n_item, n_tuple ):

    out = ();
    for it in range( 0, n_tuple ):
        token = ();
        for ii in range( 0, n_item ):
            token = token + (read_token( file ),);
#            print token[ii];
        if (len( token[ii] ) == 0):
            output( "Corrupted file (token %d)!" % ii )
            raise "ERR_CorruptedFile"
        out = out + (token,);

    return out;

##
# c: 12.02.2004, r: 20.03.2008
def read_array( file, n_row, n_col, dtype ):
    """n_col is basically ignored, n_col == -1 -> intentionally ignored."""
    val = []
    for ir in range( 0, n_row ):
        try:
            while 1:
                line = file.readline().split()
                if (len( line ) == 0) or (line[0] == "#"):
                    continue
                else:
                    break
        except:
            output( "Array (%d, %d) reading failed!" % (n_row, n_col) )
            raise
        row = [float( ii ) for ii in line]
        val.append( row )

    val = array( val, dtype );
    return val

##
# c: 05.02.2008, r: 05.02.2008
def read_list( fd, n_item, dtype ):
    vals = []
    ii = 0
    while ii < n_item:
        line = [dtype( ic ) for ic in fd.readline().split()]
        vals.append( line )
        ii += len( line )
    if ii > n_item:
        output( 'corrupted row?', line, ii, n_item  )
        raise ValueError

    return vals

##
# 27.09.2006, c
def write_dict_hdf5( file_name, adict, level = 0, group = None, fd = None ):

    if level == 0:
        fd = pt.open_file( file_name, mode = "w",
                          title = "Recursive dict dump" )
        group = '/'

    for key, val in adict.iteritems():
#        print level * ' ', key, '->', group
        if isinstance( val, dict ):
            group2 = fd.create_group( group, '_' + str( key ), '%s group' % key )
            write_dict_hdf5( file_name, val, level + 1, group2, fd )
        else:
            fd.create_array( group, '_' + str( key ), val, '%s data' % key )
            
    if level == 0:
        fd.close()

##
# 09.07.2007, c
def read_dict_hdf5( file_name, level = 0, group = None, fd = None ):
    out = {}

    if level == 0:
        fd = pt.open_file( file_name, mode = "r" )
        group = fd.root

    for name, gr in group._v_groups.iteritems():
#        print level * ' ', gr, '->', group
        name = name.replace( '_', '' )
        out[name] = read_dict_hdf5( file_name, level + 1, gr, fd )

    for name, data in group._v_leaves.iteritems():
        name = name.replace( '_', '' )
        out[name] = data.read()

    if level == 0:
        fd.close()
        
    return out

##
# 02.07.2007, c
def write_sparse_matrix_hdf5( file_name, mtx, name = 'a sparse matrix' ):
    """Assume CSR/CSC."""
    fd = pt.open_file( file_name, mode = "w", title = name )
    try:
        info = fd.create_group( '/', 'info' )
        fd.create_array( info, 'dtype', mtx.dtype.str )
        fd.create_array( info, 'shape', mtx.shape )
        fd.create_array( info, 'format', mtx.format )

        data = fd.create_group( '/', 'data' )
        fd.create_array( data, 'data', mtx.data )
        fd.create_array( data, 'indptr', mtx.indptr )
        fd.create_array( data, 'indices', mtx.indices )

    except:
        print 'matrix must be in SciPy sparse CSR/CSC format!'
        print mtx.__repr__()
        raise

    fd.close()

##
# 02.07.2007, c
# 08.10.2007
def read_sparse_matrix_hdf5( file_name, output_format = None ):
    import scipy.sparse as sp
    constructors = {'csr' : sp.csr_matrix, 'csc' : sp.csc_matrix}
    
    fd = pt.open_file( file_name, mode = "r" )
    info = fd.root.info
    data = fd.root.data

    format = info.format.read()
    if not isinstance( format, str ):
        format = format[0]

    dtype = info.dtype.read()
    if not isinstance( dtype, str ):
        dtype = dtype[0]

    if output_format is None:
        constructor = constructors[format]
    else:
        constructor = constructors[output_format]

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

    if output_format in ['csc', 'csr']:
        mtx.sort_indices()
    
    return mtx
