import numpy as nm
from numpy import array
import re, sys, os
import fnmatch
import os.path as op
from base import output, Struct, pause, dict_from_keys_init
try:
    import tables as pt
except:
    pt = None

class InDir(Struct):
    """
    Store the directory name a file is in, and prepend this name to other
    files.
    
    Examples
    --------

    >>> indir = InDir('output/file1')
    >>> print indir('file2')
    """
    def __init__(self, filename):
        self.dir = op.split(op.join(os.getcwd(), filename))[0]

    def __call__(self, filename):
         return op.join(self.dir, filename)

def ensure_path(filename):
    """
    Check if path to `filename` exists and if not, create the necessary
    intermediate directories.
    """
    dirname = os.path.dirname(filename)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def locate_files(pattern, root_dir=os.curdir):
    """
    Locate all files matching fiven filename pattern in and below
    supplied root directory.
    """
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(root_dir)):
        for filename in fnmatch.filter(filenames, pattern):
            yield os.path.join(dirpath, filename)

##
# 27.04.2006, c
def get_trunk( filename ):
    return op.splitext( op.basename( filename ) )[0]

def get_print_info(n_step, fill=None):
    """
    Returns the max. number of digits in range(n_step) and the corresponding
    format string.
    
    Examples:
    
    >>> get_print_info(11)
    (2, '%2d')
    >>> get_print_info(8)
    (1, '%1d')
    >>> get_print_info(100)
    (2, '%2d')
    >>> get_print_info(101)
    (3, '%3d')
    >>> get_print_info(101, fill='0')
    (3, '%03d')
    """    
    if n_step > 1:
        n_digit = int(nm.log10(n_step - 1) + 1)
        if fill is None:
            format = '%%%dd' % n_digit
        else:
            format = '%%%s%dd' % (fill, n_digit)
    else:
        n_digit, format = 0, None
    return n_digit, format

def skip_read_line(fd, no_eof=False):
    """
    Read the first non-empty line (if any) from the given file
    object. Return an empty string at EOF, if `no_eof` is False. If it
    is True, raise the EOFError instead.
    """
    ls = ''
    while 1:
        try:
            line = fd.readline()

        except EOFError:
            break

        if not line:
            if no_eof:
                raise EOFError

            else:
                break

        ls = line.strip()
        if ls and (ls[0] != '#'):
            break

    return ls

def read_token(fd):
    """
    Read a single token (sequence of non-whitespace characters) from the
    given file object.

    Notes
    -----
    Consumes the first whitespace character after the token.
    """
    out = "";
    # Skip initial whitespace.

    while (1):
        ch = fd.read(1)
        if (ch.isspace()): continue
        elif (len(ch) == 0): return out
        else: break

    while (not(ch.isspace())):
        out = out + ch;
        ch = fd.read(1)
        if (len(ch) == 0): break

    return out

def read_array(fd, n_row, n_col, dtype):
    """
    Read a NumPy array of shape `(n_row, n_col)` from the given file
    object and cast it to type `dtype`.
    """
    count = n_row * n_col
    val = nm.fromfile(fd, sep=' ', count=count)

    if val.shape[0] < count:
        raise ValueError("(%d, %d) array reading failed!" % (n_row, n_col))

    val = nm.asarray(val, dtype=dtype)
    val.shape = (n_row, n_col)

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

def write_dict_hdf5( filename, adict, level = 0, group = None, fd = None ):

    if level == 0:
        fd = pt.openFile( filename, mode = "w",
                          title = "Recursive dict dump" )
        group = '/'

    for key, val in adict.iteritems():
#        print level * ' ', key, '->', group
        if isinstance( val, dict ):
            group2 = fd.createGroup( group, '_' + str( key ), '%s group' % key )
            write_dict_hdf5( filename, val, level + 1, group2, fd )
        else:
            fd.createArray( group, '_' + str( key ), val, '%s data' % key )
            
    if level == 0:
        fd.close()

def read_dict_hdf5( filename, level = 0, group = None, fd = None ):
    out = {}

    if level == 0:
        fd = pt.openFile( filename, mode = "r" )
        group = fd.root

    for name, gr in group._v_groups.iteritems():
#        print level * ' ', gr, '->', group
        name = name.replace( '_', '', 1 )
        out[name] = read_dict_hdf5( filename, level + 1, gr, fd )

    for name, data in group._v_leaves.iteritems():
        name = name.replace( '_', '', 1 )
        out[name] = data.read()

    if level == 0:
        fd.close()
        
    return out

##
# 02.07.2007, c
def write_sparse_matrix_hdf5( filename, mtx, name = 'a sparse matrix' ):
    """Assume CSR/CSC."""
    fd = pt.openFile( filename, mode = "w", title = name )
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
def read_sparse_matrix_hdf5( filename, output_format = None ):
    import scipy.sparse as sp
    constructors = {'csr' : sp.csr_matrix, 'csc' : sp.csc_matrix}
    
    fd = pt.openFile( filename, mode = "r" )
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
        mtx = constructor((data.data.read(),
                           data.indices.read(), data.indptr.read()),
                          shape=info.shape.read(), dtype=dtype)
    elif format == 'coo':
        mtx = constructor((data.data.read(),
                           nm.c_[data.rows.read(), data.cols.read()].T),
                          shape=info.shape.read(), dtype=dtype)
    else:
        print format
        raise ValueError
    fd.close()

    if output_format in ['csc', 'csr']:
        mtx.sort_indices()
    
    return mtx
