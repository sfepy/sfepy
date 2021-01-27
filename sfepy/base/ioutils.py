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
import pickle
import warnings
import scipy.sparse as sp

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

def save_options(filename, options_groups, save_command_line=True,
                 quote_command_line=False):
    """
    Save groups of options/parameters into a file.

    Each option group has to be a sequence with two items: the group name and
    the options in ``{key : value}`` form.
    """
    with open(filename, 'w') as fd:
        if save_command_line:
            fd.write('command line\n')
            fd.write('------------\n\n')
            if quote_command_line:
                fd.write(' '.join('"%s"' % ii for ii in sys.argv) + '\n')

            else:
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
    if sys.version_info > (3, 0):
        string = string.encode(encoding)
    return string

def dec(val, encoding='utf-8'):
    """
    Decode given bytes using the specified encoding.
    """
    if isinstance(val, bytes) and sys.version_info > (3, 0):
        val = val.decode(encoding)
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
                fd.create_array(group, '_' + str(key), enc(val),
                                '%s data' % key)

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
    with pt.open_file(filename, mode='w', title=name) as fd:
        write_sparse_matrix_to_hdf5(fd, fd.root, mtx)

def write_sparse_matrix_to_hdf5(fd, group, mtx):
    """
    Write sparse matrix to given data group of hdf5 file

    Parameters
    ----------
    group: tables.group.group
        The hdf5 file group the matrix will be read from.
    mtx: scipy.sparse.base.spmatrix
        The writed matrix
    """
    try:
        info = fd.create_group(group, 'info')
        fd.create_array(info, 'dtype', enc(mtx.dtype.str))
        fd.create_array(info, 'shape', mtx.shape)
        fd.create_array(info, 'format', enc(mtx.format))

        data = fd.create_group(group, 'data')
        fd.create_array(data, 'data', mtx.data)
        fd.create_array(data, 'indptr', mtx.indptr)
        fd.create_array(data, 'indices', mtx.indices)

    except:
        print('matrix must be in SciPy sparse CSR/CSC format!')
        print(mtx.__repr__())
        raise

##
# 02.07.2007, c
# 08.10.2007
def read_sparse_matrix_hdf5(filename, output_format=None):

    with pt.open_file(filename, mode='r') as fd:
        out = read_sparse_matrix_from_hdf5(fd, fd.root, output_format)
    return out

def read_sparse_matrix_from_hdf5(fd, group, output_format=None):
    """
    Read sparse matrix from given data group of hdf5 file

    Parameters
    ----------
    fd: tables.File
        The hdf5 file handle the matrix will be read from.
    group: tables.group.group
        The hdf5 file group of the file the matrix will be read from.
    output_format: {'csr', 'csc', None}, optional
        The resulting matrix will be in CSR or CSC format
        if this parameter is not None (which is default),
        otherwise it will be in the format the matrix was
        stored.

    Returns
    -------
    scipy.sparse.base.spmatrix
        Readed matrix
    """
    info = group.info
    data = group.data

    constructors = {'csr' : sp.csr_matrix, 'csc' : sp.csc_matrix}

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

    if output_format in ['csc', 'csr']:
        mtx.sort_indices()

    return mtx

def path_of_hdf5_group(group):
    return group._v_pathname

class HDF5BaseData(object):
    """
    When storing values to HDF5, special classes can be used that wrap the
    stored data and modify the way the storing is done. This class is the base
    of those.
    """

    def unpack_data(self):
        """
        One can request unpacking of the wrappers during saving.

        Returns
        -------
        object
            The original object, if possible, or self.
        """
        return self

class DataMarker(HDF5BaseData):
    """
    The Base class for classes for marking data to be handled in a special way
    during saving to a HDF5 file by write_to_hdf5(). The usage is simple: just
    "decorate" the desired data element, e.g.::

        data = [data1, Cached(data2)]
        write_to_hdf5(... , ... , data)
    """

    def __init__(self, data):
        self.data = data

    def unpack_data(self):
        return self.data

class Cached(DataMarker):
    """
    The wrapper class that marks data, that should be checked
    during saving, whether it has been stored to the hdf5 file already
    and if so, a softlink to the already created instance is created
    instead of saving.
    """

class Uncached(DataMarker):
    """
    The wrapper class that marks data, that should be always stored
    to the hdf5 file, even if the object has been already stored at a
    different path in the file and so it would have been stored
    by a softlink otherwise (IGDomain, Mesh and sparse matrices behave
    so).
    """

class HDF5Data(HDF5BaseData):
    """
    Some data written to the HDF5 file can have a custom format.
    Descendants of this class should have the method `.write_data()`
    or redefine the `.write()` method.
    """

    def write(self, fd, group, name, cache=None):
        """
        Write a data structure to the HDF5 file.

        Create the following structure in the HDF5 file:
        {type: self.get_type(), anything writed by self.write_data()}

        Parameters
        ----------
        fd: tables.File
            The hdf5 file handle the data should be writed in.
        group: tables.group.Group
            The group the data will be stored to
        name: str
            Name of node that will be appended to group and will contain
            the data
        cache: dict or None, optional
            Store for already cached objects with structs id(obj) : /path/to
            Can be used for not storing the one object twice.
        """
        dgroup = fd.create_group(group, name)
        fd.create_array(dgroup, 'type', nm.array(self.get_type()))
        self.write_data(fd, dgroup, cache)
        return dgroup

    def write_data(self, fd, group):
        """
        Write data to the HDF5 file. Redefine this function in sub-classes.

        Parameters
        ----------
        fd: tables.File
            The hdf5 file handle the data should be writed to.
        group: tables.group.Group
            The group the data should be stored to.
        """
        raise Exception('Unimplemented')

class SoftLink(HDF5Data):
    """
    This object is written to the HDF5 file as a softlink to the given path.
    """
    def __init__(self, destination):
        """
        Parameters
        ----------
        destination: str
            The link destination.
        """
        self.destination = destination

    def write(self, fd, group, name, cache=None):
        """Create the softlink to the destination."""
        return fd.create_soft_link(group, name, self.destination)

class DataSoftLink(HDF5Data):
    """
    This object is written to the HDF5 file as a softlink to the given path.
    The destination of the softlink should contain only data, so
    the structure {type: type, data: softlink_to(destination)}
    is created in the place where the softlink is written.
    """

    def __init__(self, type, destination, cache=None):
        """
        Parameters
        ----------
        type: str
            Type of the object that should be written.
            See `read_from_hdf5()` for the possibilities.
        destination: str
            The Destination of the soft link where the data are stored.
        cache: object, True or None, default None
            - If None, do nothing special.
            - If True or an object cache with destination as key, storing
              such a softlink to several places does not lead to several copies
              of the object.
            - If an object is given, just identify the softlink with the
              object.
        """
        self.type = type
        self.destination = destination
        self.cache = cache

    def unpack_data(self):
        if self.cache is not None and self.cache is not True:
            return self.cache
        return self

    def get_type(self):
        return self.type

    def write_data(self, fd, group, cache=None):
        """Create the softlink to the destination and handle the caching."""

        if self.cache:
            keys = self.destination
            if self.cache is not True:
                if cache is not None:
                    cache_id = id(self.cache)
                    if cache_id in cache:
                        keys = (keys, cache[cache_id])
                    else:
                        cache[cache_id] = path_of_hdf5_group(group)

            fd.create_array(group, 'cache', nm.array(keys))
        return fd.create_soft_link(group, "data", self.destination)

def write_to_hdf5(fd, group, name, data, cache=None,
                  unpack_markers=False):
    """
    Save custom data to a HDF5 file group to be restored by read_from_hdf5().

    Allows saving lists, dicts, numpy arrays, scalars, sparse matrices,
    meshes and iga domains and all pickleable objects.

    Parameters
    ----------
    fd: tables.File
        The hdf5 file handle the data should be written in.
    group: tables.group.Group
        The group the data will be stored to.
    name: str
        The name of the node that will be appended to the group and will
        contain the data.
    data: object
        Data to be stored in the HDF5 file.
    cache: dict or None
        The cache where the paths to stored objects (currently meshes and iga
        domains) are stored, so subsequent attempts to store such objects create
        only softlinks to the initially stored object.
        The id() of objects serve as the keys into the cache.
        Mark the object with Cached() or Uncached() for (no) softlinking.
    unpack_markers:
        If True, the input data is modified so that Cached and Uncached markers
        are removed from all sub-elements of the data.

    Returns
    -------
    tables.group.Group
        The HDF5 group the data was stored to.
    """

    #imports must be done here due to circular references
    from sfepy.discrete.iga.domain import IGDomain
    from sfepy.discrete.fem.mesh import Mesh
    from sfepy.discrete.fem.meshio import HDF5MeshIO

    if cache is None: cache = {}

    def _write_to_hdf5(group, name, data):

        def save_value(type, data):
            dgroup = fd.create_group(group, name)
            fd.create_array(dgroup, 'type', nm.array(type))
            fd.create_array(dgroup, 'data', data)
            return dgroup

        def save_dict(type, data):
            dgroup = fd.create_group(group, name)
            fd.create_array(dgroup, 'type', nm.array(type))
            data_group = fd.create_group(dgroup, 'data')
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for d in data:
                    _write_to_hdf5(data_group, d, data[d])
            if unpack_markers:
                for d in data:
                    if isinstance(data[d], HDF5BaseData):
                        data[d] = data[d].unpack_data()
            return dgroup

        def save_list(type, data):
            dgroup = fd.create_group(group, name)
            fd.create_array(dgroup, 'type', nm.array(type))
            fd.create_array(dgroup, 'len', len(data))
            data_group = fd.create_group(dgroup, 'data')

            #suppress warning that nodes with numeric ids isn't accessible as
            #python attribute
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                for i, d in enumerate(data):
                    _write_to_hdf5(data_group, str(i), d)
            if unpack_markers:
                for i,d in enumerate(data):
                    if isinstance(d, DataMarker):
                        data[i] = d.data
            return dgroup

        def save_by_function(type, fn=None, data_arg=None, cached=False):
            """
            If type is None (cached must be true), try to cache data_arg,
            otherwise save just save the data_arg. Otherwise save the data
            using the given function or just save the raw data to an array.
            """
            if cached:
                obj = data if type is not None else data_arg
                if id(obj) in cache:
                    return fd.create_soft_link(group, name,cache[id(obj)])
                if type is None:
                    #store the argument into the group
                    out = _write_to_hdf5(group, name, data_arg)
                    cache[id(data_arg)] = path_of_hdf5_group(out)
                    return out

            dgroup = fd.create_group(group, name)
            cache[id(data)] = path_of_hdf5_group(dgroup)
            fd.create_array(dgroup, 'type', nm.array(type))
            if fn is None:
                fd.create_array(dgroup, 'data', data)
            else:
                data_group = fd.create_group(dgroup, 'data')
                if data_arg is None:
                    fn(fd, data_group)
                else:
                    fn(fd, data_group, data_arg)
            return dgroup

        def save_type(type):
            dgroup = fd.create_group(group, name)
            fd.create_array(dgroup, 'type', nm.array(type))
            return dgroup

        for t in (True, False, None):
            if data is t:
                return save_type(str(t))

        if isinstance(data, dict):
            return save_dict('dict', data)

        if isinstance(data, list):
            return save_list('list', data)

        if isinstance(data, tuple):
            return save_list('tuple', data)

        if isinstance(data, (int, float, complex, nm.ndarray)):
            return save_value('raw', data)

        if isinstance(data, str):
            return save_value('str', nm.array(enc(data)))

        if isinstance(data, (sp.csr_matrix, sp.csc_matrix)):
            return save_by_function('sparse_matrix',
                                    write_sparse_matrix_to_hdf5,
                                    data, cached=True)

        if isinstance(data, IGDomain):
            return save_by_function('IGDomain', data.write_domain_to_hdf5,
                                    cached=True)

        if isinstance(data, Mesh):
            return save_by_function('Mesh', HDF5MeshIO.write_mesh_to_hdf5,
                                    data_arg=data, cached=True)

        if isinstance(data, Uncached):
            return write_to_hdf5(fd, group, name, data.data)

        if isinstance(data, Cached):
            return save_by_function(None, None, data.data, cached=True)

        if isinstance(data, HDF5Data):
            return data.write(fd, group, name, cache)

        if isinstance(data, Struct):
            return save_dict('Struct', data.to_dict())

        return save_value('pickle', nm.array(pickle.dumps(data)))

    return _write_to_hdf5(group, name, data)

def read_from_hdf5(fd, group, cache=None):
    """
    Read custom data from a HDF5 file group saved by write_to_hdf5().

    The data are stored in a general (possibly nested) structure:
    {
      'type' : string type identificator
      'data' : stored data
      'cache': string, optional - another posible location of object
    }

    Parameters
    ----------
    fd: tables.File
        The hdf5 file handle the data should be restored from.
    group: tables.group.Group
        The group in the hdf5 file the data will be restored from.
    cache: dict or None
        Some objects (e.g. Mesh instances) can be stored on more places in the
        HDF5 file tree using softlinks, so when the data are restored, the
        restored objects are stored and searched in cache so that they are
        created only once. The keys to cache are the (real) paths of the
        created objects.
        Moreover, if some stored object has a 'cache' key (see e.g.
        DataSoftLink class), and the object with a given 'path' has been
        already created, it is returned instead of creating a new object.
        Otherwise, the newly created object is associated both with its real
        path and with the cache key path.

        The caching is not active for scalar data types.

    Returns
    -------
    data : object
        The restored custom data.
    """
    if cache is None:
       cache = {}

    from sfepy.discrete.iga.domain import IGDomain
    from sfepy.discrete.fem.meshio import HDF5MeshIO
    types = {b'True' : True, b'False' : False, b'None' : None}

    def _read_from_hdf5(group):
        while isinstance(group, pt.link.SoftLink):
            group = group()
        path = path_of_hdf5_group(group)
        if path in cache:
            return cache[path]

        cache_ids = tuple()
        if 'cache' in group:
            cache_ids = group.cache.read().reshape(-1)
            cache_ids = [dec(i) for i in cache_ids]
            #for DataSoftlink, there can be 2 paths, one for the softlink target
            #and one for the object identified with the softlink
            for cache_id in cache_ids:
                if cache_id in cache:
                    out = cache[cache_id]
                    break
            else:
                out = _read_value_from_hdf5(group)
        else:
            out = _read_value_from_hdf5(group)
        if not isinstance(out, (int, float, str, bool, None.__class__)):
           cache[path] = out
           for cache_id in cache_ids:
               cache[cache_id] = out
        return out

    def _read_value_from_hdf5(group):

        type = group.type.read().item()
        if type in types:
            return types[type]

        data = group.data
        while isinstance(data, pt.link.SoftLink):
            data = data()

        def load_list():
            out = [None] * group.len.read()
            for i in data:
                out[int(i._v_name)] = _read_from_hdf5(i)
            return out

        def load_dict():
            out = {}
            for i in fd.iter_nodes(data):
                out[i._v_name] = _read_from_hdf5(i)
            return out

        def load_object(fn=None):
            if fn:
                out = fn(fd, data)
            else:
                out = data.read()
            return out

        if type == b'raw':
            return data.read()

        if type == b'object':
            return load_object()

        if type == b'str':
            return dec(data.read().item())

        if type == b'pickle':
            return pickle.loads(data.read().item())

        if type == b'dict':
            return load_dict()

        if type == b'Struct':
            return Struct(**load_dict())

        if type == b'list':
            return load_list()

        if type == b'tuple':
            return tuple(load_list())

        if type == b'sparse_matrix':
            return load_object(read_sparse_matrix_from_hdf5)

        if type == b'IGDomain':
            return load_object(IGDomain.read_domain_from_hdf5)

        if type == b'Mesh':
            return load_object(HDF5MeshIO.read_mesh_from_hdf5)

        raise Exception('Unknown h5 group type {}'.format(type.decode('utf8')))

    return _read_from_hdf5(group)

class HDF5ContextManager:
    def __init__(self, filename, *args, **kwargs):
        self.filename = filename
        self.file = None
        self.args = args
        self.kwargs = kwargs

    def __enter__(self):
        if isinstance(self.filename, pt.File):
           return self.filename
        else:
           self.file = pt.open_file(self.filename, *self.args, **self.kwargs)
           return self.file

    def __exit__(self, type, value, traceback):
        if self.file:
           self.file.close()
           self.file = None

def get_or_create_hdf5_group(fd, path, from_group=None):
    if from_group is None:
       from_group = fd.root
    if path == '':
       return from_group
    if path[0] == '/':
       path = path[1:]
    for name in path.split('/'):
       if name in from_group:
           from_group = getattr(from_group, name)
       else:
           from_group = pt.create_group(from_group, name)
    return from_group
