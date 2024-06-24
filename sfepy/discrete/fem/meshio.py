import sys
from copy import copy
from io import StringIO
from packaging import version
import numpy as nm

from sfepy.base.base import (complex_types, dict_from_keys_init,
                             assert_, is_derived_class, ordered_iteritems,
                             insert_static_method, output, get_default,
                             get_default_attr, Struct, basestr)
from sfepy.base.ioutils import (skip_read_line, look_ahead_line, read_token,
                                read_array, pt, enc, dec,
                                edit_filename,
                                read_from_hdf5, write_to_hdf5,
                                HDF5ContextManager, get_or_create_hdf5_group)

import os.path as op
import six
from six.moves import range
import meshio as meshiolib
try:
    from meshio import CellBlock as meshio_Cells  # for meshio >= 4.0.3
except:
    from meshio import Cells as meshio_Cells  # for 4.0.3 > meshio > 4.0.0

_supported_formats = {
    # format name: IO class, suffix, modes[, variants]
    # modes: r = read, w = write, c = test cell groups, v = test vertex groups
    'abaqus': ('meshio', None, 'cv'),
    'avs-ucd': ('meshio', None, ''),
    'cgns': ('meshio', None, ''),
    'dolfin-xml': ('meshio', None, ''),
    'exodus': ('meshio', None, 'v'),
    'flac3d': ('meshio', None, ''),
    'gmsh': ('meshio', None, 'cv', ['gmsh4-binary', 'gmsh4-ascii',
                                    'gmsh2-binary', 'gmsh2-ascii']),
    'h5m': ('meshio', None, ''),
    'mdpa': ('meshio', None, ''),
    'med': ('meshio', None, 'cv'),
    'medit': ('meshio', None, 'cv'),
    'nastran': ('meshio', None, 'cv'),
    'netgen': ('meshio', None, ''),
    'obj': ('meshio', None, ''),
    'permas': ('meshio', None, ''),
    'su2': ('meshio', None, ''),
    'tetgen': ('meshio', None, ''),
    'ugrid': ('meshio', None, ''),
    'vtk': ('meshio', None, 'cv', ['vtk-binary', 'vtk-ascii']),
    'vtu': ('meshio', None, 'cv'),
    'wtk': ('meshio', None, ''),
    'xdmf': ('meshio', None, 'cv'),
    # The '*' prevents overriding meshio in ext2io dict, see
    # any_from_filename().
    'gmsh-dg': ('gmshio', '.msh', '*rw'),
    'ansys': ('ansys_cdb', '.cdb', 'r'),
    'hdf5': ('hdf5', '.h5', 'rwcv'),
    'hdf5-xdmf': ('hdf5-xdmf', '.h5x', 'rwcv'),
    'xyz': ('xyz', '.xyz', 'rw'),
    'comsol': ('comsol', '.txt', 'r'),
    'hmascii': ('hmascii', '.hmascii', 'r'),
    'gambit': ('gambit', '.neu', 'r'),
    'mesh3d': ('mesh3d', '.mesh3d', 'r'),
}


def update_supported_formats(formats):
    from meshio._helpers import reader_map, _writer_map
    f2e = {}

    try:
        from meshio._helpers import extension_to_filetype
        for k, v in extension_to_filetype.items():
            f2e.setdefault(v, []).append(k)
    except:
        from meshio._helpers import extension_to_filetypes
        for k, v in extension_to_filetypes.items():
            f2e.setdefault(v[-1], []).append(k)

    out = {}
    for format, info in formats.items():
        io, ext, _flag = info[:3]
        variants = info[3] if len(info) >= 4 else []
        for f in [format] + variants:
            if io == 'meshio':
                flag = _flag[:]
                if ext is None:
                    if format in f2e:
                        ext = f2e[format]
                    else:
                        continue

                if format in _writer_map:
                    flag = 'w' + flag
                if format in reader_map:
                    flag = 'r' + flag
                if not f == format:
                    flag = '*' + flag  # label format variants
            else:
                flag = _flag

            if not isinstance(ext, list):
                ext = [ext]

            out[f] = (io, ext, flag)

    return out


supported_formats = update_supported_formats(_supported_formats)
del _supported_formats


def output_mesh_formats(mode='r'):
    for key, vals in ordered_iteritems(supported_formats):
        if mode in vals[2]:
            output('%s (%s)' % (key,
                   vals[1] if len(vals[1]) > 1 else vals[1][0]))


def split_conns_mat_ids(conns_in):
    """
    Split connectivities (columns except the last ones in `conns_in`) from cell
    groups (the last columns of `conns_in`).
    """
    conns, mat_ids = [], []
    for conn in conns_in:
        conn = nm.asarray(conn, dtype=nm.int32)
        conns.append(conn[:, :-1])
        mat_ids.append(conn[:, -1])

    return conns, mat_ids


def convert_complex_output(out_in):
    """
    Convert complex values in the output dictionary `out_in` to pairs of
    real and imaginary parts.
    """
    out = {}
    for key, val in six.iteritems(out_in):

        if val.data.dtype in complex_types:
            rval = copy(val)
            rval.data = val.data.real
            out['real.%s' % key] = rval

            ival = copy(val)
            ival.data = val.data.imag
            out['imag.%s' % key] = ival

        else:
            out[key] = val

    return out


def check_format_suffix(file_format, suffix):
    """
    Check compatibility of a mesh file format and a mesh file suffix.
    """
    if file_format is None:
        return

    try:
        suffixes = supported_formats[file_format][1]

    except KeyError:
        raise ValueError('unknown file format! (%s)'
                         % file_format)

    if suffix is None:
        return

    suffix = suffix if suffix.startswith('.') else '.' + suffix
    if suffix not in suffixes:
        raise ValueError('"%s" format is not compatible with "%s" suffix!'
                         % (file_format, suffix))


class MeshIO(Struct):
    """
    The abstract class for importing and exporting meshes.

    Read the docstring of the Mesh() class. Basically all you need to do is to
    implement the read() method::

        def read(self, mesh, **kwargs):
            nodes = ...
            ngroups = ...
            conns = ...
            mat_ids = ...
            descs = ...
            mesh._set_io_data(nodes, ngroups, conns, mat_ids, descs)
            return mesh

    See the Mesh class' docstring how the nodes, ngroups, conns, mat_ids and
    descs should look like. You just need to read them from your specific
    format from disk.

    To write a mesh to disk, just implement the write() method and use the
    information from the mesh instance (e.g. nodes, conns, mat_ids and descs)
    to construct your specific format.

    Optionally, subclasses can implement read_data() to read also computation
    results. This concerns mainly the subclasses with implemented write()
    supporting the 'out' kwarg.

    The default implementation od read_last_step() just returns 0. It should be
    reimplemented in subclasses capable of storing several steps.
    """
    format = None
    call_msg = 'called an abstract MeshIO instance!'

    def __init__(self, filename, **kwargs):
        Struct.__init__(self, filename=filename, **kwargs)
        self.set_float_format()

    def get_filename_trunk(self):
        if isinstance(self.filename, basestr):
            trunk = op.splitext(self.filename)[0]
        else:
            trunk = 'from_descriptor'

        return trunk

    def read_last_step(self):
        """The default implementation: just return 0 as the last step."""
        return 0

    def read_times(self, filename=None):
        """
        Read true time step data from individual time steps.

        Returns
        -------
        steps : array
            The time steps.
        times : array
            The times of the time steps.
        nts : array
            The normalized times of the time steps, in [0, 1].

        Notes
        -----
        The default implementation returns empty arrays.
        """
        aux = nm.array([0.0], dtype=nm.float64)
        return aux.astype(nm.int32), aux, aux

    def read(self, mesh, omit_facets=False, **kwargs):
        raise ValueError(MeshIO.call_msg)

    def write(self, filename, mesh, **kwargs):
        raise ValueError(MeshIO.call_msg)

    def read_data(self, step, filename=None, cache=None):
        raise ValueError(MeshIO.call_msg)

    def set_float_format(self, format=None):
        self.float_format = get_default(format, '%e')

    def get_vector_format(self, dim):
        return ' '.join([self.float_format] * dim)


class UserMeshIO(MeshIO):
    """
    Special MeshIO subclass that enables reading and writing a mesh using a
    user-supplied function.
    """
    format = 'function'

    def __init__(self, filename, **kwargs):
        assert_(hasattr(filename, '__call__'))
        self.function = filename

        MeshIO.__init__(self, filename='function:%s' % self.function.__name__,
                        **kwargs)

    def get_filename_trunk(self):
        return self.filename

    def read(self, mesh, *args, **kwargs):
        aux = self.function(mesh, mode='read')
        if aux is not None:
            mesh = aux

        self.filename = mesh.name

        return mesh

    def write(self, filename, mesh, *args, **kwargs):
        self.function(mesh, mode='write')


def _suppress_meshio_warnings(f):
    def __suppress_meshio_warnings(*args, **kwargs):
        old_stderr = sys.stderr
        sys.stderr = mystderr = StringIO()
        try:
            out = f(*args, **kwargs)
            to_stderr = [k for k in mystderr.getvalue().split('\n')
                         if k.startswith('Error')]
            sys.stderr = old_stderr
            if len(to_stderr) > 0:
                output('\n'.join(to_stderr))
        except:
            to_stderr = mystderr.getvalue()
            if len(to_stderr) > 0:
                output(to_stderr)
            sys.stderr = old_stderr
            raise

        return out

    return __suppress_meshio_warnings


def _decorate_all(module, decorator):
    import types
    for name in dir(module):
        obj = getattr(module, name)
        if isinstance(obj, types.FunctionType):
            setattr(module, name, decorator(obj))


_decorate_all(meshiolib, _suppress_meshio_warnings)
del _decorate_all, _suppress_meshio_warnings


class MeshioLibIO(MeshIO):
    format = 'meshio'

    cell_types = {
        ('hexahedron', 3): '3_8',
        ('wedge', 3): '3_6',
        ('tetra', 3): '3_4',
        ('triangle', 3): '2_3',
        ('triangle', 2): '2_3',
        ('quad', 3): '2_4',
        ('quad', 2): '2_4',
        ('line', 3): '1_2',
        ('line', 2): '1_2',
        ('line', 1): '1_2',
    }

    def __init__(self, filename, file_format=None, **kwargs):
        MeshIO.__init__(self, filename=filename, **kwargs)
        import pathlib

        if file_format is None:
            try:
                from meshio._helpers import _filetype_from_path
                file_format = _filetype_from_path(pathlib.Path(filename))
            except:
                from meshio._helpers import _filetypes_from_path
                file_format = _filetypes_from_path(pathlib.Path(filename))[-1]

        self.file_format = file_format
        # If True, ignore axes with all-equal coordinates.
        self.squeeze = self.file_format in ('vtk', 'vtu')

    def read_bounding_box(self, ret_dim=False):
        m = meshiolib.read(self.filename, file_format=self.file_format)
        dim = self._get_dimension(m.points, squeeze=self.squeeze)

        points = m.points[:, :dim]
        bbox = nm.vstack([nm.amin(points, 0), nm.amax(points, 0)])

        if ret_dim:
            return bbox, dim
        else:
            return bbox

    @staticmethod
    def _get_dimension(points, squeeze=True):
        if squeeze: # Ignore axes with all-equal coordinates.
            dim = nm.sum(nm.max(points, axis=0)
                         - nm.min(points, axis=0) > 1e-15)

        else:
            dim = points.shape[1]

        return dim

    def read_dimension(self):
        m = meshiolib.read(self.filename, file_format=self.file_format)
        dim = self._get_dimension(m.points, squeeze=self.squeeze)

        return dim

    def read(self, mesh, omit_facets=False, **kwargs):
        m = meshiolib.read(self.filename, file_format=self.file_format)
        dim = self._get_dimension(m.points, squeeze=self.squeeze)
        ngkey = None
        for k in m.point_data.keys():
            if k == 'node_groups' or k.endswith(':ref'):
                ngkey = k
                break

        if ngkey is not None:
            ngroups = nm.asarray(m.point_data[ngkey]).flatten()
        elif hasattr(m, 'point_sets') and len(m.point_sets) > 0:
            ngroups = nm.zeros((len(m.points),), dtype=nm.int32)
            keys = list(m.point_sets.keys())
            keys.sort()
            try:
                ngrps = [int(ii) for ii in keys]
            except:
                ngrps = nm.arange(len(keys)) + 1

            for ik, k in enumerate(keys):
                ngroups[m.point_sets[k]] = ngrps[ik]
        else:
            ngroups = None

        cells, cgroups, cell_types = [], [], []

        # meshio.__version__ > 3.3.2
        cgkey = None
        for k in list(m.cell_data.keys()):
            if k == 'mat_id' or k.endswith(':ref'):
                cgkey = k
                break

        if cgkey is not None:
            cgdata = m.cell_data[cgkey]
        elif len(m.cell_sets) > 0:
            cgdata = []
            keys = list(m.cell_sets.keys())
            keys.sort()
            try:
                cgrps = [int(ii) for ii in keys]
            except:
                cgrps = nm.arange(len(keys)) + 1

            for ic, c in enumerate(m.cells):
                cgdata0 = nm.zeros((len(c.data),), dtype=nm.int32)
                for ik, k in enumerate(keys):
                    cgdata0[m.cell_sets[k][ic]] = cgrps[ik]
                cgdata.append(cgdata0)
        else:
            cgdata = None

        for ic, c in enumerate(m.cells):
            if (c.type, dim) not in self.cell_types:
                output('warning: unknown cell type %s with dimension %d'
                       % (c.type, dim))

                continue

            if (not omit_facets) or (c.dim == dim):
                cells.append(c.data)
                cell_types.append(self.cell_types[(c.type, dim)])

            if cgdata is not None:
                cgroups.append(nm.asarray(cgdata[ic]).flatten())
            else:
                cgroups.append(nm.ones((len(c.data),), dtype=nm.int32))

        mesh._set_io_data(m.points[:,:dim], ngroups,
                          cells, cgroups, cell_types)

        output('number of vertices: %d' % m.points.shape[0])
        output('number of cells:')
        for ii, k in enumerate(cell_types):
            output('  %s: %d' % (k, cells[ii].shape[0]))

        return mesh

    def write(self, filename, mesh, out=None, ts=None, **kwargs):
        (coors, cells,
         point_data,
         point_sets,
         cell_data,
         cell_sets) = self._create_out_data(mesh, out,
                                            format=self.file_format)

        if version.parse(meshiolib.__version__) >= version.parse('4.0.3') and\
           ('-ascii' in self.file_format or '-binary' in self.file_format):
            self.file_format, ab_str = self.file_format.split('-')
            kwargs['binary'] = True if 'binary' in ab_str else False

        args0 = {'point_data': point_data, 'point_sets': point_sets,
                 'cell_data': cell_data, 'cell_sets': cell_sets,
                 'file_format': self.file_format}
        args = args0.copy()
        args.update(kwargs)
        try:
            meshiolib.write_points_cells(filename, coors, cells,
                                         **args)
        except TypeError:
            meshiolib.write_points_cells(filename, coors, cells,
                                         **args0)

    def _create_out_data(self, mesh, out, format=None):
        inv_cell_types = {v: k for k, v in self.cell_types.items()}
        coors, ngroups, conns, _, descs = mesh._get_io_data()
        out = {} if out is None else out
        point_data = {k: v.data for k, v in out.items() if v.mode == 'vertex'}
        cell_data_keys = [k for k, v in out.items() if v.mode == 'cell']
        if self.file_format in ['vtk', 'vtk-ascii', 'vtk-binary', 'vtu']:
            ngkey = 'node_groups'
            cgkey = 'mat_id'
            if coors.shape[1] != 3:
                nnd, ndim = coors.shape
                ndim = 3 - ndim
                coors = nm.hstack([coors,
                                   nm.zeros((nnd, ndim), dtype=nm.float64)])
            ngroups = nm.asarray(ngroups, dtype=nm.int64)
        else:
            ngkey = '%s:ref' % self.file_format
            cgkey = '%s:ref' % self.file_format
        point_data[ngkey] = ngroups
        point_sets = {str(k): nm.where(ngroups == k)[0]
                      for k in nm.unique(ngroups)}
        cell_groups = [mesh.get_cmesh(desc).cell_groups for desc in descs]
        cgrps = nm.unique(nm.hstack(cell_groups))
        # meshio.__version__ > 3.3.2
        cells = []
        cgroups = []
        cell_data = {k: [] for k in cell_data_keys}
        cell_sets = {str(k): [] for k in cgrps}
        for ii, desc in enumerate(descs):
            cmesh = mesh.get_cmesh(desc)
            cells.append(meshio_Cells(inv_cell_types[desc][0], conns[ii]))
            cidxs = nm.where(cmesh.cell_types == cmesh.key_to_index[desc])
            cidxs = cidxs[0].astype(nm.uint32)

            cgroups.append(cmesh.cell_groups[cidxs])
            for k in cell_data_keys:
                cell_data[k].append(out[k].data[cidxs, 0, :, 0])

            for k in cgrps:
                idxs = nm.where(cmesh.cell_groups[cidxs] == k)[0]
                cell_sets[str(k)].append(cidxs[idxs])

        if self.file_format in ['vtk', 'vtk-ascii', 'vtk-binary', 'vtu']:
            point_sets, cell_sets = None, None
            cgroups = [nm.asarray(k, dtype=nm.int64) for k in cgroups]

        cell_data[cgkey] = cgroups

        return coors, cells, point_data, point_sets, cell_data, cell_sets

    def read_data(self, step, filename=None, cache=None):
        """
        Renames cell resp. vertex data with name "*:ref"
        to mat_id resp. node_groups

        Parameters
        ----------
        step: has no effect
        filename : string, optional
            The file name to use instead of self.filename.
        cache: has no effect


        Returns
        -------
        out : dictionary
            Data loaded from file, keys are names. values are Structs with
            name repeated, mode ('vertex' or 'cell') and the data itself.
        """

        filename = get_default(filename, self.filename)
        m = meshiolib.read(filename, file_format=self.file_format)
        dim = self._get_dimension(m.points, squeeze=self.squeeze)

        def _fix_shape(data):
            if data.ndim == 2:
                data = data[:, :dim]

            elif data.ndim == 3:
                data = data[:, None, ...]

            return data

        out = {}
        for key, data in m.point_data.items():
            aux = _fix_shape(data).astype(nm.float64)
            if key.endswith(':ref'):
                key = 'node_groups'
            out[key] = Struct(name=key, mode='vertex', data=aux)

        for key, data in m.cell_data.items():
            aux = _fix_shape(data[0]).astype(nm.float64)
            if key.endswith(':ref'):
                key = 'mat_id'
            out[key] = Struct(name=key, mode='cell', data=aux)

        return out


class ComsolMeshIO(MeshIO):
    format = 'comsol'

    def _read_commented_int(self):
        return int(skip_read_line(self.fd).split('#')[0])

    def _skip_comment(self):
        read_token(self.fd)
        self.fd.readline()

    def read(self, mesh, **kwargs):

        self.fd = fd = open(self.filename, 'r')
        mode = 'header'

        coors = conns = None
        while 1:
            if mode == 'header':
                line = skip_read_line(fd)

                n_tags = self._read_commented_int()
                for ii in range(n_tags):
                    skip_read_line(fd)
                n_types = self._read_commented_int()
                for ii in range(n_types):
                    skip_read_line(fd)

                skip_read_line(fd)
                assert_(skip_read_line(fd).split()[1] == 'Mesh')
                skip_read_line(fd)
                dim = self._read_commented_int()
                assert_((dim == 2) or (dim == 3))
                n_nod = self._read_commented_int()
                i0 = self._read_commented_int()
                mode = 'points'

            elif mode == 'points':
                self._skip_comment()
                coors = read_array(fd, n_nod, dim, nm.float64)
                mode = 'cells'

            elif mode == 'cells':

                n_types = self._read_commented_int()
                conns = []
                descs = []
                mat_ids = []
                for it in range(n_types):
                    t_name = skip_read_line(fd).split()[1]
                    n_ep = self._read_commented_int()
                    n_el = self._read_commented_int()

                    self._skip_comment()
                    aux = read_array(fd, n_el, n_ep, nm.int32)
                    if t_name == 'tri':
                        conns.append(aux)
                        descs.append('2_3')
                        is_conn = True
                    elif t_name == 'quad':
                        # Rearrange element node order to match SfePy.
                        aux = aux[:, (0, 1, 3, 2)]
                        conns.append(aux)
                        descs.append('2_4')
                        is_conn = True
                    elif t_name == 'hex':
                        # Rearrange element node order to match SfePy.
                        aux = aux[:, (0, 1, 3, 2, 4, 5, 7, 6)]
                        conns.append(aux)
                        descs.append('3_8')
                        is_conn = True
                    elif t_name == 'tet':
                        conns.append(aux)
                        descs.append('3_4')
                        is_conn = True
                    else:
                        is_conn = False

                    # Skip parameters.
                    n_pv = self._read_commented_int()
                    n_par = self._read_commented_int()
                    for ii in range(n_par):
                        skip_read_line(fd)

                    n_domain = self._read_commented_int()
                    assert_(n_domain == n_el)
                    if is_conn:
                        self._skip_comment()
                        mat_id = read_array(fd, n_domain, 1, nm.int32)
                        mat_ids.append(mat_id.squeeze())
                    else:
                        for ii in range(n_domain):
                            skip_read_line(fd)

                    # Skip up/down pairs.
                    n_ud = self._read_commented_int()
                    for ii in range(n_ud):
                        skip_read_line(fd)
                break

        fd.close()
        self.fd = None

        mesh._set_io_data(coors, None, conns, mat_ids, descs)

        return mesh

    def write(self, filename, mesh, out=None, **kwargs):

        def write_elements(fd, ig, conn, mat_ids, type_name,
                           npe, format, norder, nm_params):
            fd.write("# Type #%d\n\n" % ig)
            fd.write("%s # type name\n\n\n" % type_name)
            fd.write("%d # number of nodes per element\n" % npe)
            fd.write("%d # number of elements\n" % conn.shape[0])
            fd.write("# Elements\n")
            for ii in range(conn.shape[0]):
                nn = conn[ii]  # Zero based
                fd.write(format % tuple(nn[norder]))
            fd.write("\n%d # number of parameter values per element\n"
                     % nm_params)
            # Top level always 0?
            fd.write("0 # number of parameters\n")
            fd.write("# Parameters\n\n")
            fd.write("%d # number of domains\n"
                     % sum([mi.shape[0] for mi in mat_ids]))
            fd.write("# Domains\n")
            for mi in mat_ids:
                # Domains in comsol have to be > 0
                if (mi <= 0).any():
                    mi += mi.min() + 1
                for dom in mi:
                    fd.write("%d\n" % abs(dom))
            fd.write("\n0 # number of up/down pairs\n")
            fd.write("# Up/down\n")

        fd = open(filename, 'w')

        coors, ngroups, conns, mat_ids, desc = mesh._get_io_data()

        n_nod, dim = coors.shape

        # Header
        fd.write("# Created by SfePy\n\n\n")
        fd.write("# Major & minor version\n")
        fd.write("0 1\n")
        fd.write("1 # number of tags\n")
        fd.write("# Tags\n")
        fd.write("2 m1\n")
        fd.write("1 # number of types\n")
        fd.write("# Types\n")
        fd.write("3 obj\n\n")

        # Record
        fd.write("# --------- Object 0 ----------\n\n")
        fd.write("0 0 1\n") # version unused serializable
        fd.write("4 Mesh # class\n")
        fd.write("1 # version\n")
        fd.write("%d # sdim\n" % dim)
        fd.write("%d # number of mesh points\n" % n_nod)
        fd.write("0 # lowest mesh point index\n\n")  # Always zero in SfePy

        fd.write("# Mesh point coordinates\n")

        format = self.get_vector_format(dim) + '\n'
        for ii in range(n_nod):
            nn = tuple(coors[ii])
            fd.write(format % tuple(nn))

        fd.write("\n%d # number of element types\n\n\n" % len(conns))

        for ig, conn in enumerate(conns):
            if (desc[ig] == "2_4"):
                write_elements(fd, ig, conn, mat_ids,
                               "4 quad", 4, "%d %d %d %d\n", [0, 1, 3, 2], 8)

            elif (desc[ig] == "2_3"):
                # TODO: Verify number of parameters for tri element
                write_elements(fd, ig, conn, mat_ids,
                               "3 tri", 3, "%d %d %d\n", [0, 1, 2], 4)

            elif (desc[ig] == "3_4"):
                # TODO: Verify number of parameters for tet element
                write_elements(fd, ig, conn, mat_ids,
                               "3 tet", 4, "%d %d %d %d\n", [0, 1, 2, 3], 16)

            elif (desc[ig] == "3_8"):
                write_elements(fd, ig, conn, mat_ids,
                               "3 hex", 8, "%d %d %d %d %d %d %d %d\n",
                               [0, 1, 3, 2, 4, 5, 7, 6], 24)

            else:
                raise ValueError('unknown element type! (%s)' % desc[ig])

        fd.close()

        if out is not None:
            for key, val in six.iteritems(out):
                raise NotImplementedError


class HDF5MeshIO(MeshIO):
    format = "hdf5"

    import string
    _all = ''.join(map(chr, list(range(256))))
    _letters = string.ascii_letters + string.digits + '_'
    _rubbish = ''.join([ch for ch in set(_all) - set(_letters)])
    if sys.version_info[0] >= 3:
        _tr = str.maketrans(_rubbish, '_' * len(_rubbish))
    else:
        _tr = string.maketrans(_rubbish, '_' * len(_rubbish))

    @staticmethod
    def read_mesh_from_hdf5(filename, group=None, mesh=None):
        """
        Read the mesh from a HDF5 file.

        filename: str or tables.File
            The HDF5 file to read the mesh from.
        group: tables.group.Group or str, optional
            The HDF5 file group to read the mesh from.
            If None, the root group is used.
        mesh: sfepy.dicrete.fem.Mesh or None
            If None, the new mesh is created and returned, otherwise
            content of this argument is replaced by the read mesh.

        Returns
        -------
        sfepy.dicrete.fem.Mesh
            readed mesh
        """
        with HDF5ContextManager(filename, mode='r') as fd:
            if group is None:
                group = fd.root
            elif not isinstance(group, pt.group.Group):
                group = fd.get_node(group)

            set_shape_info = mesh is None
            if mesh is None:
                from .mesh import Mesh
                mesh = Mesh('mesh')

            mesh.name = dec(group.name.read())
            coors = group.coors.read()
            ngroups = group.ngroups.read()

            n_gr = group.n_gr.read()

            conns = []
            descs = []
            mat_ids = []
            for ig in range(n_gr):
                gr_name = 'group%d' % ig
                conn_group = group._f_get_child(gr_name)
                conns.append(conn_group.conn.read())
                mat_ids.append(conn_group.mat_id.read())
                descs.append(dec(conn_group.desc.read()))

            nodal_bcs = {}
            try:
                node_sets_groups = group.node_sets

            except:
                pass

            else:
                for group in node_sets_groups:
                    key = dec(group.key.read())
                    nods = group.nods.read()
                    nodal_bcs[key] = nods

            mesh._set_io_data(coors, ngroups, conns, mat_ids, descs,
                              nodal_bcs=nodal_bcs)

            if set_shape_info:
                mesh._set_shape_info()
        return mesh

    @staticmethod
    def write_mesh_to_hdf5(filename, group, mesh, force_3d=False):
        """
        Write mesh to a hdf5 file.

        filename: str or tables.File
            The HDF5 file to write the mesh to.
        group: tables.group.Group or None or str
            The HDF5 file group to write the mesh to.
            If None, the root group is used.
            The group can be given as a path from root, e.g. /path/to/mesh
        mesh: sfepy.dicrete.fem.Mesh
            The mesh to write.
        """
        with HDF5ContextManager(filename, mode='w') as fd:
            if group is None:
                group = fd.root
            elif not isinstance(group, pt.group.Group):
                group = get_or_create_hdf5_group(fd, group)

            coors, ngroups, conns, mat_ids, descs = mesh._get_io_data()
            if force_3d and coors.shape[1] == 2:
                coors = nm.hstack([coors, nm.zeros((len(coors), 1))])
            fd.create_array(group, 'name', enc(mesh.name), 'name')
            fd.create_array(group, 'coors', coors, 'coors')
            fd.create_array(group, 'ngroups', ngroups, 'ngroups')
            fd.create_array(group, 'n_gr', len(conns), 'n_gr')
            for ig, conn in enumerate(conns):
                conn_group = fd.create_group(group, 'group%d' % ig,
                                             'connectivity group')
                fd.create_array(conn_group, 'conn', conn, 'connectivity')
                fd.create_array(conn_group, 'mat_id', mat_ids[ig],
                                'material id')
                fd.create_array(conn_group, 'desc', enc(descs[ig]),
                                'element Type')

            node_sets_groups = fd.create_group(group, 'node_sets',
                                               'node sets groups')
            ii = 0
            for key, nods in six.iteritems(mesh.nodal_bcs):
                group = fd.create_group(node_sets_groups, 'group%d' % ii,
                                        'node sets group')
                fd.create_array(group, 'key', enc(key), 'key')
                fd.create_array(group, 'nods', nods, 'nods')
                ii += 1

    def read_dimension(self, ret_fd=False):
        fd = pt.open_file(self.filename, mode="r")

        dim = fd.root.mesh.coors.shape[1]

        if ret_fd:
            return dim, fd

        else:
            fd.close()
            return dim

    def read_bounding_box(self, ret_fd=False, ret_dim=False):
        fd = pt.open_file(self.filename, mode="r")

        mesh_group = fd.root.mesh

        coors = mesh_group.coors.read()
        bbox = nm.vstack((nm.amin(coors, 0),
                          nm.amax(coors, 0)))

        if ret_dim:
            dim = coors.shape[1]

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

    def read(self, mesh=None, **kwargs):
        return self.read_mesh_from_hdf5(self.filename, '/mesh', mesh=mesh)

    @staticmethod
    def write_xdmf_file(filename, **kwargs):
        def get_path(node):
            _, fname = op.split(filename)
            return '%s:%s' % (fname, node._v_pathname)

        def get_data_dim(shape):
            if len(shape) == 4:
                return shape[2]
            if len(shape) == 2:
                return shape[1]
            else:
                return 1

        def data_item(data):
            dtype = data.dtype
            if nm.issubdtype(dtype, nm.integer):
                data_type = 'Int'
            elif nm.issubdtype(dtype, nm.floating):
                data_type = 'Float'
            else:
                raise ValueError('wrong data type! (%s)' % dtype)

            dim = get_data_dim(data.shape)
            sh = (data.shape[0], dim)
            ditem = et.Element('DataItem',
                               attrib={'DataType': data_type,
                                       'Dimensions': '%d %d' % sh,
                                       'Format': 'HDF',
                                       'Precision': str(dtype.itemsize)})
            ditem.text = get_path(data)

            return ditem

        def attr_item(data, center=None, name=None):
            if isinstance(data, pt.Group):
                if center is None and 'mode' in data:
                    mode = data.mode.read().decode('ascii')
                    if mode == 'custom':
                        return None
                    center = {'vertex': 'Node',
                              'cell': 'Cell'}[mode]
                if name is None and 'dname' in data:
                    name = data.dname.read().decode('ascii').lstrip('_')
                data = data.data

            dim = get_data_dim(data.shape)
            atype = {1: 'Scalar', 3: 'Vector', 6: 'Tensor6', 9: 'Tensor'}[dim]

            aitem = et.Element('Attribute',
                               attrib={'AttributeType': atype,
                                       'Center': center,
                                       'Name': name})
            aitem.append(data_item(data))

            return aitem

        import xml.etree.ElementTree as et
        from xml.dom import minidom

        topology_table = {
            '2_2': 'Line',
            '3_2': 'Line',
            '2_3': 'Triangle',
            '2_4': 'Quadrilateral',
            '3_4': 'Tetrahedron',
            '3_8': 'Hexahedron',
        }

        et_root = et.Element('Xdmf', attrib={'Version': '3.0'})
        if 'extra_data' in kwargs:
            for k, v in kwargs['extra_data'].items():
                d = et.SubElement(et_root, 'Information', attrib={'Name': k})
                d.text = str(v)

        et_domain = et.SubElement(et_root, 'Domain')

        with HDF5ContextManager(filename, mode='r') as fd:
            root = fd.root
            mesh = root.mesh
            name = mesh.name.read().decode('ascii')
            et_mesh = []
            geom = et.Element('Geometry', attrib={'GeometryType': 'XYZ'})
            geom.append(data_item(mesh.coors))
            et_mesh.append(geom)
            et_mesh.append(attr_item(mesh.ngroups, 'Node', 'node_groups'))

            n_gr = mesh.n_gr.read()
            for ig in range(n_gr):
                gr_name = 'group%d' % ig
                conn_group = mesh._f_get_child(gr_name)
                nc, nnd = conn_group.conn.shape
                ttype = topology_table[dec(conn_group.desc.read())]
                et_conn = et.Element('Topology',
                                     attrib={'NumberOfElements': str(nc),
                                             'TopologyType': ttype})
                et_conn.append(data_item(conn_group.conn))
                et_mesh.append(et_conn)
                et_mesh.append(attr_item(conn_group.mat_id, 'Cell', 'mat_id'))

            steps = [k for k in root if k._v_name.startswith('step')]
            et_ts = et.SubElement(et_domain, 'Grid',
                                  attrib={'Name': 'TimeSeries',
                                          'GridType': 'Collection',
                                          'CollectionType': 'Temporal'})

            for step in steps:
                istep = int(step._v_name[4:])
                et_grid = et.SubElement(et_ts, 'Grid',
                                        attrib={'Name': 'grid%d' % istep,
                                                'GridType': 'Uniform'})
                et.SubElement(et_grid, 'Time', attrib={'Value': '%d' % istep})
                et_grid.extend(et_mesh)

                for val in filter(lambda x: x._v_name.startswith('__'), step):
                    aitem = attr_item(val)
                    if aitem is not None:
                        et_grid.append(aitem)

        out = minidom.parseString(et.tostring(et_root)).toprettyxml(indent="  ")
        xdmf_filename = op.splitext(filename)[0] + '.xdmf'
        with open(xdmf_filename, 'w') as f:
            f.write(out[(out.find('\n') + 1):])

    def write(self, filename, mesh, out=None, ts=None, cache=None,
              xdmf=False, **kwargs):
        def expand_data_3d(data):
            expand_tab = {
                2: (3, [0, 1]),
                3: (6, [0, 1, 3]),
                4: (9, [0, 1, 3, 4]),
            }

            n, dim = data.shape
            if dim in expand_tab:
                dim3, order = expand_tab[dim]
                out = nm.zeros((n, dim3), dtype=data.dtype)
                out[:, order] = data
                return out
            else:
                return data

        from time import asctime

        if pt is None:
            raise ValueError('pytables not imported!')

        step = get_default_attr(ts, 'step', 0)
        if (step == 0) or not op.exists(filename):
            # A new file.
            with pt.open_file(filename, mode="w",
                              title="SfePy output file") as fd:
                mesh_group = fd.create_group('/', 'mesh', 'mesh')
                self.write_mesh_to_hdf5(fd, mesh_group, mesh, force_3d=xdmf)

                if ts is not None:
                    ts_group = fd.create_group('/', 'ts', 'time stepper')
                    fd.create_array(ts_group, 't0', ts.t0, 'initial time')
                    fd.create_array(ts_group, 't1', ts.t1, 'final time')
                    fd.create_array(ts_group, 'dt', ts.dt, 'time step')
                    fd.create_array(ts_group, 'n_step', ts.n_step, 'n_step')

                tstat_group = fd.create_group('/', 'tstat',
                                              'global time statistics')
                fd.create_array(tstat_group, 'created', enc(asctime()),
                                'file creation time')
                fd.create_array(tstat_group, 'finished', enc('.' * 24),
                                'file closing time')

                fd.create_array(fd.root, 'last_step',
                                nm.array([step], dtype=nm.int32),
                                'last saved step')

        if out is not None:
            if ts is None:
                step, time, nt = 0, 0.0, 0.0
            else:
                step, time, nt = ts.step, ts.time, ts.nt

            # Existing file.
            fd = pt.open_file(filename, mode="r+")

            step_group_name = 'step%d' % step
            if step_group_name in fd.root:
                raise ValueError('step %d is already saved in "%s" file!'
                                 ' Possible help: remove the old file or'
                                 ' start saving from the initial time.'
                                 % (step, filename))
            step_group = fd.create_group('/', step_group_name, 'time step data')

            ts_group = fd.create_group(step_group, 'ts', 'time stepper')
            fd.create_array(ts_group, 'step', step, 'step')
            fd.create_array(ts_group, 't', time, 'time')
            fd.create_array(ts_group, 'nt', nt, 'normalized time')

            name_dict = {}
            for key, val in six.iteritems(out):
                if xdmf and mesh.coors.shape[1] == 2:
                    data = expand_data_3d(val.data)
                else:
                    data = val.data

                group_name = '__' + key.translate(self._tr)
                data_group = fd.create_group(step_group, group_name,
                                             '%s data' % key)
                fd.create_array(data_group, 'dname', enc(key), 'data name')
                fd.create_array(data_group, 'mode', enc(val.mode), 'mode')
                name = val.get('name', 'output_data')
                fd.create_array(data_group, 'name', enc(name), 'object name')
                if val.mode == 'custom':
                    write_to_hdf5(fd, data_group, 'data', data,
                                  cache=cache,
                                  unpack_markers=getattr(val, 'unpack_markers',
                                                         False))
                    continue

                shape = val.get('shape', data.shape)
                dofs = val.get('dofs', None)
                if dofs is None:
                    dofs = [''] * nm.squeeze(shape)[-1]
                var_name = val.get('var_name', '')

                fd.create_array(data_group, 'data', data, 'data')
                fd.create_array(data_group, 'dofs', [enc(ic) for ic in dofs],
                                'dofs')
                fd.create_array(data_group, 'shape', shape, 'shape')
                fd.create_array(data_group, 'var_name',
                                enc(var_name), 'object parent name')
                if val.mode == 'full':
                    fd.create_array(data_group, 'field_name',
                                    enc(val.field_name), 'field name')

                reg_name = val.get('region_name', '')
                fd.create_array(data_group, 'region_name',
                                enc(reg_name), 'region name')

                name_dict[key] = group_name

            step_group._v_attrs.name_dict = name_dict
            fd.root.last_step[0] = step

            fd.remove_node(fd.root.tstat.finished)
            fd.create_array(fd.root.tstat, 'finished', enc(asctime()),
                            'file closing time')
            fd.close()

        if xdmf:
            self.write_xdmf_file(filename, **kwargs)

    def read_last_step(self, filename=None):
        filename = get_default(filename, self.filename)
        fd = pt.open_file(filename, mode="r")
        last_step = fd.root.last_step[0]
        fd.close()
        return last_step

    def read_time_stepper(self, filename=None):
        filename = get_default(filename, self.filename)
        fd = pt.open_file(filename, mode="r")

        try:
            ts_group = fd.root.ts
            out = (ts_group.t0.read(), ts_group.t1.read(),
                   ts_group.dt.read(), ts_group.n_step.read())

        except:
            raise ValueError('no time stepper found!')

        finally:
            fd.close()

        return out

    def _get_step_group_names(self, fd):
        return sorted([name for name in fd.root._v_groups.keys()
                       if name.startswith('step')],
                      key=lambda name: int(name[4:]))

    def read_times(self, filename=None):
        """
        Read true time step data from individual time steps.

        Returns
        -------
        steps : array
            The time steps.
        times : array
            The times of the time steps.
        nts : array
            The normalized times of the time steps, in [0, 1].
        """
        filename = get_default(filename, self.filename)
        fd = pt.open_file(filename, mode='r')

        steps = []
        times = []
        nts = []
        for gr_name in self._get_step_group_names(fd):
            ts_group = fd.get_node(fd.root, gr_name + '/ts')

            steps.append(ts_group.step.read())
            times.append(ts_group.t.read())
            nts.append(ts_group.nt.read())
        fd.close()

        steps = nm.asarray(steps, dtype=nm.int32)
        times = nm.asarray(times, dtype=nm.float64)
        nts = nm.asarray(nts, dtype=nm.float64)

        return steps, times, nts

    def _get_step_group(self, step, filename=None):
        filename = get_default(filename, self.filename)
        fd = pt.open_file(filename, mode="r")

        if step is None:
            step = int(self._get_step_group_names(fd)[0][4:])

        gr_name = 'step%d' % step
        try:
            step_group = fd.get_node(fd.root, gr_name)
        except:
            output('step %d data not found - premature end of file?' % step)
            fd.close()
            return None, None

        return fd, step_group

    def read_data(self, step, filename=None, cache=None):
        fd, step_group = self._get_step_group(step, filename=filename)
        if fd is None:
            return None

        out = {}
        for data_group in step_group:
            try:
                key = dec(data_group.dname.read())

            except pt.exceptions.NoSuchNodeError:
                continue

            mode = dec(data_group.mode.read())
            if mode == 'custom':
                out[key] = read_from_hdf5(fd, data_group.data, cache=cache)
                continue

            name = dec(data_group.name.read())
            data = data_group.data.read()
            dofs = tuple([dec(ic) for ic in data_group.dofs.read()])
            try:
                shape = tuple(int(ii) for ii in data_group.shape.read())

            except pt.exceptions.NoSuchNodeError:
                shape = data.shape

            if mode == 'full':
                field_name = dec(data_group.field_name.read())

            else:
                field_name = None

            out[key] = Struct(name=name, mode=mode, data=data,
                              dofs=dofs, shape=shape, field_name=field_name)

            if out[key].dofs == (-1,):
                out[key].dofs = None

        fd.close()

        return out

    def read_data_header(self, dname, step=None, filename=None):
        fd, step_group = self._get_step_group(step, filename=filename)
        if fd is None:
            return None

        groups = step_group._v_groups
        for name, data_group in six.iteritems(groups):
            try:
                key = dec(data_group.dname.read())

            except pt.exceptions.NoSuchNodeError:
                continue

            if key == dname:
                mode = dec(data_group.mode.read())
                fd.close()
                return mode, name

        fd.close()
        raise KeyError('non-existent data: %s' % dname)

    def read_time_history(self, node_name, indx, filename=None):
        filename = get_default(filename, self.filename)
        fd = pt.open_file(filename, mode="r")

        th = dict_from_keys_init(indx, list)
        for gr_name in self._get_step_group_names(fd):
            step_group = fd.get_node(fd.root, gr_name)
            data = step_group._f_get_child(node_name).data

            for ii in indx:
                th[ii].append(nm.array(data[ii]))

        fd.close()

        for key, val in six.iteritems(th):
            aux = nm.array(val)
            if aux.ndim == 4: # cell data.
                aux = aux[:,0,:,0]
            th[key] = aux

        return th

    def read_variables_time_history(self, var_names, ts, filename=None):
        filename = get_default(filename, self.filename)
        fd = pt.open_file(filename, mode="r")

        assert_((fd.root.last_step[0] + 1) == ts.n_step)

        ths = dict_from_keys_init(var_names, list)

        arr = nm.asarray
        for step in range(ts.n_step):
            gr_name = 'step%d' % step
            step_group = fd.get_node(fd.root, gr_name)
            name_dict = step_group._v_attrs.name_dict
            for var_name in var_names:
                data = step_group._f_get_child(name_dict[var_name]).data
                ths[var_name].append(arr(data.read()))

        fd.close()

        return ths


class HDF5XdmfMeshIO(HDF5MeshIO):
    format = "hdf5-xdmf"

    def write(self, filename, mesh, out=None, ts=None, cache=None, **kwargs):
        HDF5MeshIO.write(self, filename, mesh, out=out, ts=ts, cache=cache,
                         xdmf=True, **kwargs)


class Mesh3DMeshIO(MeshIO):
    format = "mesh3d"

    def read(self, mesh, **kwargs):
        # read the whole file:
        with open(self.filename) as f:
            vertices = self._read_section(f, integer=False)
            tetras = self._read_section(f)
            hexes = self._read_section(f)
            prisms = self._read_section(f)
            tris = self._read_section(f)
            quads = self._read_section(f)

        # substract 1 from all elements, because we count from 0:
        conns = []
        mat_ids = []
        descs = []
        if len(tetras) > 0:
            conns.append(tetras - 1)
            mat_ids.append([0]*len(tetras))
            descs.append("3_4")
        if len(hexes) > 0:
            conns.append(hexes - 1)
            mat_ids.append([0]*len(hexes))
            descs.append("3_8")
        mesh._set_io_data(vertices, None, conns, mat_ids, descs)
        return mesh

    def read_dimension(self):
        return 3

    def _read_line(self, f):
        """
        Reads one non empty line (if it's a comment, it skips it).
        """
        l = f.readline().strip()
        while l == "" or l[0] == "#":  # comment or an empty line
            l = f.readline().strip()
        return l

    def _read_section(self, f, integer=True):
        """
        Reads one section from the mesh3d file.

        integer ... if True, all numbers are passed to int(), otherwise to
            float(), before returning

        Some examples how a section can look like:

        2
        1 2 5 4 7 8 11 10
        2 3 6 5 8 9 12 11

        or

        5
        1 2 3 4     1
        1 2 6 5     1
        2 3 7 6     1
        3 4 8 7     1
        4 1 5 8     1

        or

        0

        """
        if integer:
            dtype=int
        else:
            dtype=float
        l = self._read_line(f)
        N = int(l)
        rows = []
        for i in range(N):
            l = self._read_line(f)
            row = nm.fromstring(l, sep=" ", dtype=dtype)
            rows.append(row)
        return nm.array(rows)


def mesh_from_groups(mesh, ids, coors, ngroups,
                     tris, mat_tris, quads, mat_quads,
                     tetras, mat_tetras, hexas, mat_hexas, remap=None):
    ids = nm.asarray(ids, dtype=nm.int32)
    coors = nm.asarray(coors, dtype=nm.float64)

    if remap is None:
        n_nod = coors.shape[0]
        remap = nm.zeros((ids.max()+1,), dtype=nm.int32)
        remap[ids] = nm.arange(n_nod, dtype=nm.int32)

    tris = remap[nm.array(tris, dtype=nm.int32)]
    quads = remap[nm.array(quads, dtype=nm.int32)]
    tetras = remap[nm.array(tetras, dtype=nm.int32)]
    hexas = remap[nm.array(hexas, dtype=nm.int32)]

    conns = [tris, quads, tetras, hexas]
    mat_ids = [nm.array(ar, dtype=nm.int32)
               for ar in [mat_tris, mat_quads, mat_tetras, mat_hexas]]
    descs = ['2_3', '2_4', '3_4', '3_8']

    # Remove empty groups.
    conns, mat_ids, descs = zip(*[(conns[ig], mat_ids[ig], descs[ig])
                                  for ig in range(4)
                                  if conns[ig].shape[0] > 0])

    mesh._set_io_data(coors, ngroups, conns, mat_ids, descs)
    return mesh


class HypermeshAsciiMeshIO(MeshIO):
    format = 'hmascii'

    def read(self, mesh, **kwargs):
        fd = open(self.filename, 'r')

        ids = []
        coors = []
        tetras = []
        mat_tetras = []
        hexas = []
        mat_hexas = []
        quads = []
        mat_quads = []
        trias = []
        mat_trias = []
        mat_id = 0

        for line in fd:
            if line and (line[0] == '*'):
                if line[1:10] == 'component':
                    line = line.strip()[11:-1].split(',')
                    mat_id = int(line[0])
                if line[1:5] == 'node':
                    line = line.strip()[6:-1].split(',')
                    ids.append(int(line[0]))
                    coors.append([float(coor) for coor in line[1:4]])

                elif line[1:7] == 'tetra4':
                    line = line.strip()[8:-1].split(',')
                    mat_tetras.append(mat_id)
                    tetras.append([int(ic) for ic in line[2:6]])

                elif line[1:6] == 'hexa8':
                    line = line.strip()[7:-1].split(',')
                    mat_hexas.append(mat_id)
                    hexas.append([int(ic) for ic in line[2:10]])

                elif line[1:6] == 'quad4':
                    line = line.strip()[7:-1].split(',')
                    mat_quads.append(mat_id)
                    quads.append([int(ic) for ic in line[2:6]])

                elif line[1:6] == 'tria3':
                    line = line.strip()[7:-1].split(',')
                    mat_trias.append(mat_id)
                    trias.append([int(ic) for ic in line[2:5]])
        fd.close()

        mesh = mesh_from_groups(mesh, ids, coors, None,
                                trias, mat_trias, quads, mat_quads,
                                tetras, mat_tetras, hexas, mat_hexas)

        return mesh

    def read_dimension(self):
        return 3

    def write(self, filename, mesh, out=None, **kwargs):
        raise NotImplementedError


class NEUMeshIO(MeshIO):
    format = 'gambit'

    def read_dimension(self, ret_fd=False):

        fd = open(self.filename, 'r')

        row = fd.readline().split()
        while 1:
            if not row:
                break
            if len(row) == 0:
                continue

            if (row[0] == 'NUMNP'):
                row = fd.readline().split()
                n_nod, n_el, dim = row[0], row[1], int(row[4])
                break

        if ret_fd:
            return dim, fd
        else:
            fd.close()
            return dim

    def read(self, mesh, **kwargs):

        el = {'3_8': [], '3_4': [], '2_4': [], '2_3': []}
        nod = []

        conns_in = []
        descs = []

        group_ids = []
        group_n_els = []
        groups = []
        nodal_bcs = {}

        fd = open(self.filename, 'r')

        row = fd.readline()
        while 1:
            if not row:
                break
            row = row.split()

            if len(row) == 0:
                row = fd.readline()
                continue

            if (row[0] == 'NUMNP'):
                row = fd.readline().split()
                n_nod, n_el, dim = int(row[0]), int(row[1]), int(row[4])

            elif (row[0] == 'NODAL'):
                row = fd.readline().split()
                while not(row[0] == 'ENDOFSECTION'):
                    nod.append(row[1:])
                    row = fd.readline().split()

            elif (row[0] == 'ELEMENTS/CELLS'):
                row = fd.readline().split()
                while not(row[0] == 'ENDOFSECTION'):
                    elid = [row[0]]
                    gtype = int(row[1])
                    if gtype == 6:
                        el['3_4'].append(row[3:]+elid)
                    elif gtype == 4:
                        rr = row[3:]
                        if (len(rr) < 8):
                            rr.extend(fd.readline().split())
                        el['3_8'].append(rr+elid)
                    elif gtype == 3:
                        el['2_3'].append(row[3:]+elid)
                    elif gtype == 2:
                        el['2_4'].append(row[3:]+elid)
                    row = fd.readline().split()

            elif (row[0] == 'GROUP:'):
                group_ids.append(row[1])
                g_n_el = int(row[3])
                group_n_els.append(g_n_el)
                name = fd.readline().strip()

                els = []
                row = fd.readline().split()
                row = fd.readline().split()
                while not(row[0] == 'ENDOFSECTION'):
                    els.extend(row)
                    row = fd.readline().split()
                if g_n_el != len(els):
                    msg = 'wrong number of group elements! (%d == %d)'\
                          % (n_el, len(els))
                    raise ValueError(msg)
                groups.append(els)

            elif (row[0] == 'BOUNDARY'):
                row = fd.readline().split()
                key = row[0]
                num = int(row[2])
                inod = read_array(fd, num, None, nm.int32) - 1
                nodal_bcs[key] = inod.squeeze()

                row = fd.readline().split()
                assert_(row[0] == 'ENDOFSECTION')

            row = fd.readline()

        fd.close()

        if int(n_el) != sum(group_n_els):
            print('wrong total number of group elements! (%d == %d)'\
                  % (int(n_el), len(group_n_els)))

        mat_ids = nm.zeros(n_el, dtype=nm.int32)
        for ii, els in enumerate(groups):
            els = nm.array(els, dtype=nm.int32)
            mat_ids[els - 1] = group_ids[ii]

        for elem in el.keys():
            if len(el[elem]) > 0:
                els = nm.array(el[elem], dtype=nm.int32)
                els[:, :-1] -= 1
                els[:, -1] = mat_ids[els[:, -1]-1]

                if elem == '3_8':
                    els = els[:, [0, 1, 3, 2, 4, 5, 7, 6, 8]]

                conns_in.append(els)
                descs.append(elem)

        nod = nm.array(nod, nm.float64)

        conns, mat_ids = split_conns_mat_ids(conns_in)
        mesh._set_io_data(nod, None, conns, mat_ids, descs,
                          nodal_bcs=nodal_bcs)

        return mesh

    def write(self, filename, mesh, out=None, **kwargs):
        raise NotImplementedError


class ANSYSCDBMeshIO(MeshIO):
    format = 'ansys_cdb'

    @staticmethod
    def guess(filename):
        fd = open(filename, 'r')

        for ii in range(1000):
            row = fd.readline()
            if not row: break
            if len(row) == 0: continue

            row = row.split(',')
            kw = row[0].lower()

            if (kw == 'nblock'):
                ok = True
                break

        else:
            ok = False

        fd.close()

        return ok

    @staticmethod
    def make_format(format, nchar=1000):
        idx = []
        dtype = []
        start = 0

        for iform in format:
            ret = iform.partition('i')
            if not ret[1]:
                ret = iform.partition('e')
            if not ret[1]:
                raise ValueError
            aux = ret[2].partition('.')
            step = int(aux[0])
            for j in range(int(ret[0])):
                if (start + step) > nchar:
                    break
                idx.append((start, start+step))
                start += step
                dtype.append(ret[1])

        return idx, dtype

    def write(self, filename, mesh, out=None, **kwargs):
        raise NotImplementedError

    def read_bounding_box(self):
        raise NotImplementedError

    def read_dimension(self, ret_fd=False):
        return 3

    def read(self, mesh, **kwargs):
        ids = []
        coors = []
        tetras = []
        hexas = []
        qtetras = []
        qhexas = []
        nodal_bcs = {}

        fd = open(self.filename, 'r')

        while True:
            row = fd.readline()
            if not row:
                break
            if len(row) == 0:
                continue

            row = row.split(',')
            kw = row[0].lower()

            if (kw == 'nblock'):
                # Solid keyword -> 3, otherwise 1 is the starting coors index.
                ic = 3 if len(row) == 3 else 1
                fmt = fd.readline()
                fmt = fmt.strip()[1:-1].split(',')
                row = look_ahead_line(fd)
                nchar = len(row)
                idx, dtype = self.make_format(fmt, nchar)
                ii0, ii1 = idx[0]
                while True:
                    row = fd.readline()
                    if ((row[0] == '!') or (row[:2] == '-1')
                        or len(row) != nchar):
                        break

                    line = [float(row[i0:i1]) for i0, i1 in idx[ic:]]

                    ids.append(int(row[ii0:ii1]))
                    coors.append(line)

            elif (kw == 'eblock'):
                if (len(row) <= 2) or row[2].strip().lower() != 'solid':
                    continue

                fmt = fd.readline()
                fmt = [fmt.strip()[1:-1]]
                row = look_ahead_line(fd)
                nchar = len(row)
                idx, dtype = self.make_format(fmt, nchar)

                imi0, imi1 = idx[0] # Material id.
                inn0, inn1 = idx[8] # Number of nodes in line.
                ien0, ien1 = idx[10] # Element number.
                ic0 = 11
                while True:
                    row = fd.readline()
                    if ((row[0] == '!') or (row[:2] == '-1')
                        or (len(row) != nchar)):
                        break

                    line = [int(row[imi0:imi1])]
                    n_nod = int(row[inn0:inn1])

                    line.extend(int(row[i0:i1])
                                for i0, i1 in idx[ic0:ic0 + n_nod])
                    if n_nod == 4:
                        tetras.append(line)

                    elif n_nod == 8:
                        hexas.append(line)

                    elif n_nod == 10:
                        row = fd.readline()
                        line.extend(int(row[i0:i1])
                                    for i0, i1 in idx[:2])
                        qtetras.append(line)

                    elif n_nod == 20:
                        row = fd.readline()
                        line.extend(int(row[i0:i1])
                                    for i0, i1 in idx[:12])
                        qhexas.append(line)

                    else:
                        raise ValueError('unsupported element type! (%d nodes)'
                                         % n_nod)

            elif kw == 'cmblock':
                if row[2].lower() != 'node':  # Only node sets support.
                    continue

                n_nod = int(row[3].split('!')[0])
                fd.readline() # Format line not needed.

                nods = read_array(fd, n_nod, 1, nm.int32)
                nodal_bcs[row[1].strip()] = nods.ravel()

        fd.close()

        coors = nm.array(coors, dtype=nm.float64)

        tetras = nm.array(tetras, dtype=nm.int32)
        if len(tetras):
            mat_ids_tetras = tetras[:, 0]
            tetras = tetras[:, 1:]

        else:
            tetras.shape = (0, 4)
            mat_ids_tetras = nm.array([])

        hexas = nm.array(hexas, dtype=nm.int32)
        if len(hexas):
            mat_ids_hexas = hexas[:, 0]
            hexas = hexas[:, 1:]

        else:
            hexas.shape = (0, 8)
            mat_ids_hexas = nm.array([])

        if len(qtetras):
            qtetras = nm.array(qtetras, dtype=nm.int32)
            tetras.shape = (max(0, tetras.shape[0]), 4)
            tetras = nm.r_[tetras, qtetras[:, 1:5]]
            mat_ids_tetras = nm.r_[mat_ids_tetras, qtetras[:, 0]]

        if len(qhexas):
            qhexas = nm.array(qhexas, dtype=nm.int32)
            hexas.shape = (max(0, hexas.shape[0]), 8)
            hexas = nm.r_[hexas, qhexas[:, 1:9]]
            mat_ids_hexas = nm.r_[mat_ids_hexas, qhexas[:, 0]]

        if len(qtetras) or len(qhexas):
            ii = nm.union1d(tetras.ravel(), hexas.ravel())
            n_nod = len(ii)

            remap = nm.zeros((ii.max()+1,), dtype=nm.int32)
            remap[ii] = nm.arange(n_nod, dtype=nm.int32)

            ic = nm.searchsorted(ids, ii)
            coors = coors[ic]

        else:
            n_nod = coors.shape[0]
            remap = nm.zeros((nm.array(ids).max() + 1,), dtype=nm.int32)
            remap[ids] = nm.arange(n_nod, dtype=nm.int32)

        # Convert tetras as degenerate hexas to true tetras.
        ii = nm.where((hexas[:, 2] == hexas[:, 3])
                      & (hexas[:, 4] == hexas[:, 5])
                      & (hexas[:, 4] == hexas[:, 6])
                      & (hexas[:, 4] == hexas[:, 7]))[0]

        if len(ii) == len(hexas):
            tetras = nm.r_[tetras, hexas[ii[:, None], [0, 1, 2, 4]]]
            mat_ids_tetras = nm.r_[mat_ids_tetras, mat_ids_hexas[ii]]

            hexas = nm.delete(hexas, ii, axis=0)
            mat_ids_hexas = nm.delete(mat_ids_hexas, ii)

        else:
            output('WARNING: mesh "%s" has both tetrahedra and hexahedra!'
                   % mesh.name)

        ngroups = nm.zeros(len(coors), dtype=nm.int32)

        mesh = mesh_from_groups(mesh, ids, coors, ngroups,
                                [], [], [], [],
                                tetras, mat_ids_tetras,
                                hexas, mat_ids_hexas, remap=remap)

        mesh.nodal_bcs = {}
        for key, nods in six.iteritems(nodal_bcs):
            nods = nods[nods < len(remap)]
            mesh.nodal_bcs[key] = remap[nods]

        return mesh


class GmshIO(MeshioLibIO):
    """
    Used to read and write data in .msh format when file_format gmsh-dg is
    specified. Tailored for use with Discontinous galerking methods, mesh and
    ElementNodeData with InterpolationScheme can be written and read.
    It however omits mat_ids and node_groups.

    For details on format see [1].

    For details on representing and visualization of DG FEM data using gmsh see [2].

    [1] http://gmsh.info/doc/texinfo/gmsh.html#File-formats

    [2] Remacle, J.-F., Chevaugeon, N., Marchandise, E., & Geuzaine, C. (2007).
    Efficient visualization of high-order finite elements. International Journal
    for Numerical Methods in Engineering, 69(4), 750-771.
    https://doi.org/10.1002/nme.1787

    """
    format = 'gmshio'

    load_slices = {"all": slice(0, None),
                   "first": slice(0, 1),
                   "last": slice(-1, None)}

    def __init__(self, filename, file_format=None, **kwargs):
        MeshioLibIO.__init__(self, filename=filename, file_format=None,
                             **kwargs)

    def _get_filename_format(self, filename):
        try:
            basename, step_num, extension = filename.split(".")
        except ValueError:
            raise ValueError("Filename of automatically loaded GMSH data must"
                             " be: <base name>.<step number>.msh, {} has to "
                             "correspond to that"
                             .format(filename))
        n_digits = len(step_num)
        return basename + ".{:0"+str(n_digits)+"d}." + extension

    def _get_filename_wildcard(self, filename):
        try:
            basename, step_num, extension = filename.split(".")
        except ValueError:
            raise ValueError("Filename of automatically loaded GMSH data must"
                             " be: <base name>.<step number>.msh, {} has to "
                             "correspond to that"
                             .format(filename))
        return basename + ".*[0-9]." + extension

    def read_data(self, step=None, filename=None, cache=None):
        """
        Reads file or files with basename filename or self.filename. Considers
        all files to contain data from time steps of solution of single transient
        problem i.e. all data have the same shape, mesh and same interpolation
        scheme in case of ElementNodeData. Does not read mulitple
        NodeData or ElementData. For stationary problems just reads one file
        with time 0.0 and time step 0.

        Providing filename allows reading multiple files of format
        `basename.*[0-9].msh`

        Parameters
        ----------
        step : String, int,  optional
            "all", "last", "first" or number of step to read:
            if "all" read all files with the basename and varying step,
            if "last" read only last step of all files with the filename,
            if "first" reads step=0,
            if None reads file with filename provided or specified in object.
        filename : string, optional
             Filename of the files to use, if None filename from object is used.
             Basename is extracted as `basename.*[0-9].msh`
        cache : has no effect

        Returns
        -------
        out : dictionary
            Keys represent name of data, values are Structs with attributes:

            data : list, array
                For ElementNodeData with shape (n_cell, n_cell_dof) contains
                for each time step.
                For other contains array of data from last time step.
            time : list
                Contains times.
            time_n : list
                Contains time step numbers.
            scheme : Struct
                Interpolation scheme used in data, only one interpolation
                scheme is allowed.
            scheme_name : str
                Name of the interpolation scheme, repeated fo convenience.
            mode : str
                 Represents of type of data. cell_nodes : for ElementNodeData;
                 vertex or cell : Note that for vertex and cell data reading
                 multiple time steps does not work yet.

        Notes
        -----
        The interpolation scheme `Struct` contains the following items:
            name : string
                Name of the scheme.
            F : array
                Coefficients matrix as defined in [1] and [2].
            P : array
                Exponents matrix as defined in [1] and [2].
        """
        filename = get_default(filename, self.filename)

        out = {}

        def append_data_structs(struct1, struct2):
            struct = struct1 + struct2
            if hasattr(struct, "data"):
                struct.data = struct1.data + struct2.data
            if hasattr(struct, "time"):
                struct.time = struct1.time + struct2.time
            if hasattr(struct, "time_n"):
                struct.time_n = struct1.time_n + struct2.time_n
            return struct

        if step in ["all", "last", "first"]:
            import glob
            from os.path import join as pjoin
            filename_wildcard = self._get_filename_wildcard(filename)
            filenames = glob.glob(filename_wildcard)[self.load_slices[step]]

            for filename in filenames:
                element_node_out = self._read_element_node_data(filename)
                for key, val in element_node_out.items():
                    out[key] = append_data_structs(
                        out.setdefault(key, Struct(data=[],
                                                   time=[],
                                                   time_n=[])), val)

                # read vertex or cell data
                vertex_cell_out = super(GmshIO, self).read_data(step, filename)
                out.update(vertex_cell_out)

        elif isinstance(step, int) and not op.exists(filename):
            filename_format = self._get_filename_format(filename)
            filename = filename_format.format(step)
            try:
                element_node_out = self._read_element_node_data(filename)
                out.update(element_node_out)

                # read vertex or cell data
                vertex_cell_out = super(GmshIO, self).read_data(step, filename)
                out.update(vertex_cell_out)
            except FileNotFoundError as e:
                raise FileNotFoundError(str(e) +
                                        " Maybe time step {} is not in output."
                                        .format(step))

        elif step is None or op.exists(filename):
            element_node_out = self._read_element_node_data(filename)
            out.update(element_node_out)

            # read vertex or cell data
            vertex_cell_out = super(GmshIO, self).read_data(step, filename)
            out.update(vertex_cell_out)

        else:
            raise ValueError("Unsupported vaule for step : {}".format(step))
        return out

    def _read_element_node_data(self, filename):
        try:
            fd = open(filename, "r")
        except FileNotFoundError:
            raise FileNotFoundError("[Errno 2] No such file or directory: {}."
                                    .format(filename))

        out = {}
        schemes = {}
        while 1:
            line = skip_read_line(fd).split()
            if not line:
                break

            ls = line[0]
            if ls == "$InterpolationScheme":
                scheme = Struct(name=None, desc=None, F=None, P=None)
                scheme.name = skip_read_line(fd).strip('"\'')
                n_int_tags = int(skip_read_line(fd))
                scheme.desc = int(skip_read_line(fd))
                n_matrices = int(skip_read_line(fd))
                f_shape = [int(i) for i in skip_read_line(fd).split(" ")]
                scheme.F = read_array(fd, f_shape[0], f_shape[1], nm.float64)
                p_shape = [int(i) for i in skip_read_line(fd).split(" ")]
                scheme.P = read_array(fd, p_shape[0], p_shape[1], nm.float64)
                schemes[scheme.name] = scheme
            elif ls == "$ElementNodeData":
                n_str_tags = int(skip_read_line(fd))
                data_name = skip_read_line(fd).strip('"\'')
                if n_str_tags == 2:
                    scheme_name = skip_read_line(fd).strip('"\'')
                n_float_tags =  int(skip_read_line(fd))
                time = float(skip_read_line(fd))
                n_int_tags = int(skip_read_line(fd))
                time_n = int(skip_read_line(fd))
                comp = int(skip_read_line(fd))
                n_el = int(skip_read_line(fd))

                n_el_nod = int(look_ahead_line(fd).split()[1])
                # read data including indexing
                data = read_array(fd, n_el, n_el_nod + 2, nm.float64)
                # strip indexing columns
                data = data[:, 2:]

                out[data_name] = Struct(name=data_name,
                                        data=[data],
                                        time=[time],
                                        time_n=[time_n],
                                        scheme_name=scheme_name,
                                        scheme=schemes.get(scheme_name),
                                        mode="cell_nodes")
            elif line[0] == '#' or ls[:4] == '$End':
                pass
        fd.close()

        # add schemes read later than data
        for key, val in out.items():
            if val.scheme is None:
                val.scheme = schemes[val.scheme_name]

        return out

    def _write_interpolation_scheme(self, fd, scheme):
        """
        Unpacks matrices from scheme struct and writes them in correct format
        for gmsh to read.

        Parameters
        ----------
        fd :
            File opened for writing.
        scheme : Struct
            Struct with interpolation scheme used in data, only one interpolation
            scheme is allowed,
            contains :
                name - name of the scheme,
                F - coeficients matrix,
                P - exponents matrix as defined in [1] and [2].
        """
        fd.write('$InterpolationScheme\n')
        fd.write('"{}"\n'.format(scheme.name))
        fd.write("1\n")  # one int tag
        fd.write("{}\n".format(scheme.desc[-1]))
        fd.write("2\n")  # number of matrices
        fd.write("{} {}\n".format(*scheme.F.shape))
        sF = "{} " * scheme.F.shape[1] + "\n"
        for row in scheme.F:
            fd.write(sF.format(*row))
        fd.write("{} {}\n".format(*scheme.P.shape))
        sP = "{} " * scheme.P.shape[1] + "\n"
        for row in scheme.P:
            fd.write(sP.format(*row))
        fd.write('$EndInterpolationScheme\n')

    def _write_elementnode_data(self, fd, out, ts):
        """
        Writes "cell_nodes" data in out as $ElementNodeData,
        including interpolation scheme.
        """
        for key, value in out.items():
            if not value.mode == "cell_nodes":
                continue
            if value.scheme is not None:
                self._write_interpolation_scheme(fd, value.scheme)
                scheme_name = value.scheme.name
            data = value.data
            n_el_nod = nm.shape(data)[1]
            fd.write("$ElementNodeData\n")
            fd.write("{}\n".format(1 if scheme_name is None else 2))
            fd.write('"{}"\n'.format(key))  # name
            if scheme_name is not None:
                fd.write('"{}"\n'.format(scheme_name))
            fd.write("1\n")  # number of real tags
            fd.write("{}\n".format(ts.time if ts is not None else 0.0))
            fd.write("3\n")  # number of integer tags
            fd.write("{}\n".format(ts.step if ts is not None else 0))
            fd.write("1\n")  # number of components
            fd.write("{}\n".format(data.shape[0]))
            s = "{} {}" + n_el_nod * " {}" + "\n"
            for i, el_node_vals in enumerate(data, 1):
                fd.write(s.format(i, n_el_nod, *el_node_vals))
            fd.write("$EndElementNodeData\n")

    def write(self, filename, mesh, out=None, ts=None, **kwargs):
        """
        Writes mesh and data, handles cell DOFs data from DGField
        as ElementNodeData.

        Omits gmsh:ref for cells and vertices i.e. mat_ids and
        node_groups to prevent cluttering the GMSH postprocessing.

        Parameters
        ----------
        filename : string
            Path to file.
        mesh : sfepy.discrete.fem.mesh.Mesh
            Computational mesh to write.
        out : dictionary
           Keys represent name of the data, values are Structs with attributes:

           data : array
             For ElementNodeData shape is (n_cell, n_cell_dof)
           mode : str
             Represents type of data, cell_nodes for ElementNodeData.

           For ElementNodeData:

           scheme : Struct
             Interpolation scheme used in data, only one interpolation
             scheme is allowed.
           scheme_name : str
             Name of the interpolation scheme, associated with data,
             repeated fo convenience.

        ts : sfepy.solvers.ts.TimeStepper instance, optional
            Provides data to write time step.

        Notes
        -----
        The interpolation scheme `Struct` contains the following items:
            name : string
                Name of the scheme.
            F : array
                Coefficients matrix as defined in [1] and [2].
            P : array
                Exponents matrix as defined in [1] and [2].
        """
        # fd.writelines(self.msh20header)
        # self._write_mesh(fd, mesh)
        (coors, cells,
         point_data,
         point_sets,
         cell_data,
         cell_sets) = self._create_out_data(mesh, out)

        # gmsh:ref creates clutter in GMSH, especially for transient problems
        point_data.pop("gmsh:ref", None)
        cell_data.pop("gmsh:ref", None)

        meshiolib.write_points_cells(filename, coors, cells,
                                     point_data=point_data,
                                     point_sets=point_sets,
                                     cell_data=cell_data,
                                     cell_sets=cell_sets,
                                     file_format=self.file_format,
                                     binary=False)

        if out:
            with open(filename, 'a') as fd:
                self._write_elementnode_data(fd, out, ts)
        return


class XYZMeshIO(MeshIO):
    """
    Trivial XYZ format working only with coordinates (in a .XYZ file) and the
    connectivity stored in another file with the same base name and .IEN
    suffix.
    """
    format = 'xyz'

    def _read_coors(self):
        coors = nm.loadtxt(self.filename, ndmin=2)
        if (coors[:, -1] == 0).all():
            coors = coors[:, :-1].copy()

        return coors

    def read_dimension(self, ret_fd=False):
        coors = self._read_coors()
        dim = coors.shape[1]

        if ret_fd:
            fd = open(self.filename, 'r')
            return dim, fd

        else:
            return dim

    def read_bounding_box(self, ret_fd=False, ret_dim=False):
        coors = self._read_coors()
        bbox = nm.vstack((nm.amin(coors, 0), nm.amax(coors, 0)))

        if ret_fd:
            fd = open(self.filename, 'r')
        if ret_dim:
            dim = coors.shape[1]
            if ret_fd:
                return bbox, dim, fd
            else:
                return bbox, dim
        else:
            if ret_fd:
                return bbox, fd
            else:
                return bbox

    def read(self, mesh, omit_facets=False, **kwargs):
        coors = self._read_coors()
        n_nod, dim = coors.shape

        conn_ext = '.IEN' if op.splitext(self.filename)[1].isupper() else '.ien'
        conn = nm.loadtxt(edit_filename(self.filename, new_ext=conn_ext),
                          dtype=nm.int32, ndmin=2) - 1
        desc = '%d_%d' % (dim, conn.shape[1])

        mesh._set_io_data(coors, nm.zeros(n_nod, dtype=nm.int32),
                          [conn], [nm.zeros(conn.shape[0], dtype=nm.int32)],
                          [desc])
        return mesh

    def write(self, filename, mesh, out=None, **kwargs):
        coors, ngroups, conns, mat_ids, desc = mesh._get_io_data()
        n_nod, dim = coors.shape

        zz = nm.zeros((n_nod, 3-dim))
        nm.savetxt(filename, nm.c_[coors, zz])

        conn_ext = '.IEN' if op.splitext(filename)[1].isupper() else '.ien'
        nm.savetxt(edit_filename(self.filename, new_ext=conn_ext),
                   conns[0] + 1, fmt='%d')

        if out is not None:
            raise NotImplementedError


var_dict = list(vars().items())
io_table = {}

for key, var in var_dict:
    try:
        if is_derived_class(var, MeshIO):
            io_table[var.format] = var
    except TypeError:
        pass
del var_dict


def any_from_filename(filename, prefix_dir=None, file_format=None, mode='r'):
    """
    Create a MeshIO instance according to the kind of `filename`.

    Parameters
    ----------
    filename : str, function or MeshIO subclass instance
        The name of the mesh file. It can be also a user-supplied function
        accepting two arguments: `mesh`, `mode`, where `mesh` is a Mesh
        instance and `mode` is one of 'read','write', or a MeshIO subclass
        instance.
    prefix_dir : str
        The directory name to prepend to `filename`.

    Returns
    -------
    io : MeshIO subclass instance
        The MeshIO subclass instance corresponding to the kind of `filename`.
    """
    if not isinstance(filename, basestr):
        if isinstance(filename, MeshIO):
            return filename

        else:
            return UserMeshIO(filename)

    if prefix_dir is not None:
        filename = op.normpath(op.join(prefix_dir, filename))

    kwargs = {}
    if file_format is not None:
        if file_format in supported_formats:
            io_class = supported_formats[file_format][0]
            kwargs['file_format'] = file_format
        else:
            raise ValueError('unknown mesh format! (%s)' % file_format)
    else:
        ext2io = {e: (v[0], k) for k, v in supported_formats.items()
                  for e in v[1] if '*' not in v[2]}
        ext = op.splitext(filename)[1].lower()
        if ext in ext2io:
            io_class = ext2io[ext][0]
            file_format = ext2io[ext][1]
            kwargs['file_format'] = file_format
        else:
            raise ValueError('unknown mesh format! (%s)' % ext)

    if mode == 'w' and 'w' not in supported_formats[file_format][2]:
        output('writable mesh formats:')
        output_mesh_formats('w')
        msg = 'write support not implemented for output mesh format "%s",' \
              ' see above!' % file_format
        raise ValueError(msg)

    return io_table[io_class](filename, **kwargs)


insert_static_method(MeshIO, any_from_filename)
del any_from_filename
