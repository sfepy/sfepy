from __future__ import absolute_import
from __future__ import print_function
import numpy as nm
import sys
import os
import os.path as op
import fnmatch
import shutil
import glob
from .base import output, ordered_iteritems, Struct, basestr
import six
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
    if dirname:
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        if not os.path.isdir(dirname):
            raise IOError('cannot ensure path for "%s"!' % filename)

def locate_files(pattern, root_dir=os.curdir, **kwargs):
    """
    Locate all files matching fiven filename pattern in and below
    supplied root directory.

    The `**kwargs` arguments are passed to ``os.walk()``.
    """
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(root_dir),
                                                **kwargs):
        for filename in fnmatch.filter(filenames, pattern):
            yield os.path.join(dirpath, filename)

def remove_files(root_dir, **kwargs):
    """
    Remove all files and directories in supplied root directory.

    The `**kwargs` arguments are passed to ``os.walk()``.
    """
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(root_dir),
                                                **kwargs):
        for filename in filenames:
            os.remove(os.path.join(root_dir, filename))

        for dirname in dirnames:
            shutil.rmtree(os.path.join(root_dir, dirname))

def remove_files_patterns(root_dir, patterns, ignores=None,
                          verbose=False):
    """
    Remove files with names satisfying the given glob patterns in a supplied
    root directory. Files with patterns in `ignores` are omitted.
    """
    from itertools import chain

    if ignores is None: ignores = []
    for _f in chain(*[glob.glob(os.path.join(root_dir, pattern))
                      for pattern in patterns]):
        can_remove = True
        for ignore in ignores:
            if fnmatch.fnmatch(_f, os.path.join(root_dir, ignore)):
                can_remove = False
                break

        if can_remove:
            output('removing "%s"' % _f, verbose=verbose)
            os.remove(_f)

def save_options(filename, options_groups, save_command_line=True):
    """
    Save groups of options/parameters into a file.

    Each option group has to be a sequence with two items: the group name and
    the options in ``{key : value}`` form.
    """
    with open(filename, 'w') as fd:
        if save_command_line:
            fd.write('command line\n')
            fd.write('------------\n\n')
            fd.write(' '.join(sys.argv) + '\n')

        for options_group in options_groups:
            name, options = options_group
            fd.write('\n%s\n' % name)
            fd.write(('-' * len(name)) + '\n\n')
            for key, val in ordered_iteritems(options):
                fd.write('%s: %s\n' % (key, val))

def enc(string, encoding='utf-8'):
    """
    Encode given string or bytes using the specified encoding.
    """
    val = string.encode(encoding)
    return val

def dec(val, encoding='utf-8'):
    """
    Decode given bytes using the specified encoding.
    """
    if isinstance(val, bytes):
        return val.decode(encoding)

    else:
        return val

##
# 27.04.2006, c
def get_trunk(filename):
    return op.splitext(op.basename(filename))[0]

def edit_filename(filename, prefix='', suffix='', new_ext=None):
    """
    Edit a file name by add a prefix, inserting a suffix in front of a file
    name extension or replacing the extension.

    Parameters
    ----------
    filename : str
        The file name.
    prefix : str
        The prefix to be added.
    suffix : str
        The suffix to be inserted.
    new_ext : str, optional
        If not None, it replaces the original file name extension.

    Returns
    -------
    new_filename : str
        The new file name.
    """
    path, filename = os.path.split(filename)
    base, ext = os.path.splitext(filename)

    if new_ext is None:
        new_filename = prefix + base + suffix + ext

    else:
        new_filename = prefix + base + suffix + new_ext

    return os.path.join(path, new_filename)

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

def look_ahead_line(fd):
    """
    Read and return a line from the given file object. Saves the current
    position in the file before the reading occurs and then, after the reading,
    restores the saved (original) position.
    """
    lastpos = fd.tell()
    line = fd.readline()
    fd.seek(lastpos)

    return line

def read_token(fd):
    """
    Read a single token (sequence of non-whitespace characters) from the
    given file object.

    Notes
    -----
    Consumes the first whitespace character after the token.
    """
    out = ''
    # Skip initial whitespace.

    while 1:
        ch = fd.read(1)
        if ch.isspace(): continue
        elif len(ch) == 0: return out
        else: break

    while not ch.isspace():
        out = out + ch
        ch = fd.read(1)
        if len(ch) == 0: break

    return out

def read_array(fd, n_row, n_col, dtype):
    """
    Read a NumPy array of shape `(n_row, n_col)` from the given file
    object and cast it to type `dtype`.
    If `n_col` is None, determine the number of columns automatically.
    """
    if n_col is None:
        idx = fd.tell()
        row = fd.readline().split()
        fd.seek(idx)
        n_col = len(row)

    count = n_row * n_col
    val = nm.fromfile(fd, sep=' ', count=count)

    if val.shape[0] < count:
        raise ValueError('(%d, %d) array reading failed!' % (n_row, n_col))

    val = nm.asarray(val, dtype=dtype)
    val.shape = (n_row, n_col)

    return val

##
# c: 05.02.2008, r: 05.02.2008
def read_list(fd, n_item, dtype):
    vals = []
    ii = 0
    while ii < n_item:
        line = [dtype(ic) for ic in fd.readline().split()]
        vals.append(line)
        ii += len(line)
    if ii > n_item:
        output('corrupted row?', line, ii, n_item)
        raise ValueError

    return vals

def write_dict_hdf5(filename, adict, level=0, group=None, fd=None):

    if level == 0:
        fd = pt.open_file(filename, mode='w', title='Recursive dict dump')
        group = '/'

    for key, val in six.iteritems(adict):
        if isinstance(val, dict):
            group2 = fd.create_group(group, '_' + str(key), '%s group' % key)
            write_dict_hdf5(filename, val, level + 1, group2, fd)
        else:
            if not isinstance(val, basestr):
                fd.create_array(group, '_' + str(key), val, '%s data' % key)

            else:
                fd.create_array(group, '_' + str(key), enc(val), '%s data' % key)

    if level == 0:
        fd.close()

def read_dict_hdf5(filename, level=0, group=None, fd=None):
    out = {}

    if level == 0:
        fd = pt.open_file(filename, mode='r')
        group = fd.root

    for name, gr in six.iteritems(group._v_groups):
        name = name.replace('_', '', 1)
        out[name] = read_dict_hdf5(filename, level + 1, gr, fd)

    for name, data in six.iteritems(group._v_leaves):
        name = name.replace('_', '', 1)
        val = data.read()
        if isinstance(val, bytes):
            val = dec(val)

        out[name] = val

    if level == 0:
        fd.close()

    return out

##
# 02.07.2007, c
def write_sparse_matrix_hdf5(filename, mtx, name='a sparse matrix'):
    """Assume CSR/CSC."""
    fd = pt.open_file(filename, mode='w', title=name)
    try:
        info = fd.create_group('/', 'info')
        fd.create_array(info, 'dtype', enc(mtx.dtype.str))
        fd.create_array(info, 'shape', mtx.shape)
        fd.create_array(info, 'format', enc(mtx.format))

        data = fd.create_group('/', 'data')
        fd.create_array(data, 'data', mtx.data)
        fd.create_array(data, 'indptr', mtx.indptr)
        fd.create_array(data, 'indices', mtx.indices)

    except:
        print('matrix must be in SciPy sparse CSR/CSC format!')
        print(mtx.__repr__())
        raise

    fd.close()

##
# 02.07.2007, c
# 08.10.2007
def read_sparse_matrix_hdf5(filename, output_format=None):
    import scipy.sparse as sp
    constructors = {'csr' : sp.csr_matrix, 'csc' : sp.csc_matrix}

    fd = pt.open_file(filename, mode='r')
    info = fd.root.info
    data = fd.root.data

    format = dec(info.format.read())
    dtype = dec(info.dtype.read())

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
        print(format)
        raise ValueError
    fd.close()

    if output_format in ['csc', 'csr']:
        mtx.sort_indices()

    return mtx
