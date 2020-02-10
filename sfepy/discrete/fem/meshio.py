from __future__ import print_function
from __future__ import absolute_import
import sys
from copy import copy

import numpy as nm

from sfepy.base.base import (complex_types, dict_from_keys_init,
                             assert_, is_derived_class, ordered_iteritems,
                             insert_static_method, output, get_default,
                             get_default_attr, Struct, basestr)
from sfepy.base.ioutils import (skip_read_line, look_ahead_line, read_token,
                                read_array, read_list, pt, enc, dec,
                                edit_filename,
                                read_from_hdf5, write_to_hdf5,
                                HDF5ContextManager, get_or_create_hdf5_group)
import os.path as op
import six
from six.moves import range
import meshio as meshiolib

_supported_formats = {
    # format name: IO class, suffix, modes[, variants]
    # modes: r = read, w = write, c = test cell groups, v = test vertex groups
    'abaqus': ('meshio', None, 'cv'),
    'exodus': ('meshio', None, 'v'),
    # 'ansys': ('meshio', None, ''),
    'gmsh': ('meshio', None, 'cv', ['gmsh4-binary', 'gmsh4-ascii',
                                    'gmsh2-binary', 'gmsh2-ascii']),
    'medit': ('meshio', None, 'cv'),
    'nastran': ('meshio', None, 'cv'),
    'vtk': ('meshio', None, 'cv', ['vtk-binnary', 'vtk-ascii']),
    'vtu': ('meshio', None, 'cv'),
    'med': ('meshio', None, 'cv'),
    'xdmf': ('meshio', None, 'cv'),
    'tetgen': ('meshio', None, ''),
    'hdf5': ('hdf5', '.h5', 'rwcv'),
    'xyz': ('xyz', '.xyz', 'rw'),
    'comsol': ('comsol', '.txt', 'r'),
    'hmascii': ('hmascii', '.hmascii', 'r'),
    'gambit': ('gambit', '.neu', 'r'),
    'mesh3d': ('mesh3d', '.mesh3d', 'r'),
}

def update_supported_formats(formats):
    from meshio._helpers import reader_map, _writer_map,\
        _extension_to_filetype

    f2e = {}
    for k, v in _extension_to_filetype.items():
        f2e.setdefault(v, []).append(k)

    out = {}
    for format, info in formats.items():
        io, ext, _flag = info[:3]
        variants = info[3] if len(info) >= 4 else []
        for f in [format] + variants:
            if io is 'meshio':
                flag = _flag[:]
                if ext is None:
                    ext = f2e[format]
                if f in _writer_map:
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
            output('%s (%s)' % (key, vals[1]))


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

        if val.data.dtype in  complex_types:
            rval = copy(val)
            rval.data = val.data.real
            out['real.%s' % key] = rval

            ival = copy(val)
            ival.data = val.data.imag
            out['imag.%s' % key] = ival

        else:
            out[key] = val

    return out


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


class MeshioLibIO(MeshIO):
    format = 'meshio'

    cell_types = {
        ('hexahedron', 3): '3_8',
        ('tetra', 3): '3_4',
        ('triangle', 3): '2_3',
        ('triangle', 2): '2_3',
        ('quad', 3): '2_4',
        ('quad', 2): '2_4',
        ('line', 3): '3_2',
        ('line', 2): '2_2',
        ('line', 1): '1_2',
    }

    def __init__(self, filename, file_format=None, **kwargs):
        MeshIO.__init__(self, filename=filename, **kwargs)
        from meshio._helpers import _filetype_from_path
        import pathlib

        if file_format is None:
            file_format = _filetype_from_path(pathlib.Path(filename))

        self.file_format = file_format

    def read_bounding_box(self, ret_dim=False):
        m = meshiolib.read(self.filename, file_format=self.file_format)

        bbox = nm.vstack([nm.amin(m.points, 0), nm.amax(m.points, 0)])

        if ret_dim:
            return bbox, m.points.shape[1]
        else:
            return bbox

    def read_dimension(self):
        m = meshiolib.read(self.filename, file_format=self.file_format)
        dim = nm.sum(nm.max(m.points, axis=0)\
            - nm.min(m.points, axis=0) > 1e-15)

        return dim

    def read(self, mesh, omit_facets=False, **kwargs):
        m = meshiolib.read(self.filename, file_format=self.file_format)

        dim = nm.sum(nm.max(m.points, axis=0)\
            - nm.min(m.points, axis=0) > 1e-15)

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
        inv_cell_types = {v: k for k, v in self.cell_types.items()}

        coors, ngroups, conns, _, descs = mesh._get_io_data()

        out = {} if out is None else out

        point_data = {k: v.data for k, v in out.items() if v.mode == 'vertex'}
        cell_data_keys = [k for k, v in out.items() if v.mode == 'cell']

        if self.file_format in ['vtk', 'vtu']:
            ngkey = 'node_groups'
            cgkey = 'mat_id'
        else:
            ngkey = '%s:ref' % self.file_format
            cgkey = '%s:ref' % self.file_format

        point_data[ngkey] = ngroups
        point_sets = {str(k): nm.where(ngroups == k)[0]
            for k in nm.unique(ngroups)}

        cmesh = mesh.cmesh
        cell_groups = cmesh.cell_groups
        cgrps = nm.unique(cell_groups)

        # meshio.__version__ > 3.3.2
        cells = []
        cgroups = [ ]
        cell_data = {k: [] for k in cell_data_keys}
        cell_sets = {str(k): [] for k in cgrps}
        for ii, desc in enumerate(descs):
            cells.append(meshiolib.Cells(type=inv_cell_types[desc][0],
                                            data=conns[ii]))
            cidxs = nm.where(cmesh.cell_types == cmesh.key_to_index[desc])
            cidxs = cidxs[0].astype(nm.uint32)

            cgroups.append(cell_groups[cidxs])
            for k in cell_data_keys:
                cell_data[k].append(out[k].data[cidxs, 0, :, 0])

            for k in cgrps:
                idxs = nm.where(cell_groups[cidxs] == k)[0]
                cell_sets[str(k)].append(cidxs[idxs])

        cell_data[cgkey] = cgroups

        meshiolib.write_points_cells(filename, coors, cells,
                                     point_data=point_data,
                                     point_sets=point_sets,
                                     cell_data=cell_data,
                                     cell_sets=cell_sets,
                                     file_format=self.file_format)


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
                        aux = aux[:,(0,1,3,2)]
                        conns.append(aux)
                        descs.append('2_4')
                        is_conn = True
                    elif t_name == 'hex':
                        # Rearrange element node order to match SfePy.
                        aux = aux[:,(0,1,3,2,4,5,7,6)]
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
                nn = conn[ii] # Zero based
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
        fd.write("0 # lowest mesh point index\n\n") # Always zero in SfePy

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
    def write_mesh_to_hdf5(filename, group, mesh):
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

    def write(self, filename, mesh, out=None, ts=None, cache=None, **kwargs):
        from time import asctime

        if pt is None:
            raise ValueError('pytables not imported!')

        step = get_default_attr(ts, 'step', 0)
        if (step == 0) or not op.exists(filename):
            # A new file.
            with pt.open_file(filename, mode="w",
                              title="SfePy output file") as fd:
                mesh_group = fd.create_group('/', 'mesh', 'mesh')
                self.write_mesh_to_hdf5(fd, mesh_group, mesh)

                if ts is not None:
                    ts_group = fd.create_group('/', 'ts', 'time stepper')
                    fd.create_array(ts_group, 't0', ts.t0, 'initial time')
                    fd.create_array(ts_group, 't1', ts.t1, 'final time' )
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
                step, time, nt  = 0, 0.0, 0.0
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
                group_name = '__' + key.translate(self._tr)
                data_group = fd.create_group(step_group, group_name,
                                             '%s data' % key)
                fd.create_array(data_group, 'dname', enc(key), 'data name')
                fd.create_array(data_group, 'mode', enc(val.mode), 'mode')
                name = val.get('name', 'output_data')
                fd.create_array(data_group, 'name', enc(name), 'object name')
                if val.mode == 'custom':
                    write_to_hdf5(fd, data_group, 'data', val.data,
                                  cache=cache,
                                  unpack_markers=getattr(val, 'unpack_markers',
                                                         False))
                    continue

                shape = val.get('shape', val.data.shape)
                dofs = val.get('dofs', None)
                if dofs is None:
                    dofs = [''] * nm.squeeze(shape)[-1]
                var_name = val.get('var_name', '')

                fd.create_array(data_group, 'data', val.data, 'data')
                fd.create_array(data_group, 'dofs', [enc(ic) for ic in dofs],
                                'dofs')
                fd.create_array(data_group, 'shape', shape, 'shape')
                fd.create_array(data_group, 'var_name',
                                enc(var_name), 'object parent name')
                if val.mode == 'full':
                    fd.create_array(data_group, 'field_name',
                                    enc(val.field_name), 'field name')

                name_dict[key] = group_name

            step_group._v_attrs.name_dict = name_dict
            fd.root.last_step[0] = step

            fd.remove_node(fd.root.tstat.finished)
            fd.create_array(fd.root.tstat, 'finished', enc(asctime()),
                            'file closing time')
            fd.close()

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
            out =  (ts_group.t0.read(), ts_group.t1.read(),
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
        if fd is None: return None

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
        if fd is None: return None

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


class Mesh3DMeshIO(MeshIO):
    format = "mesh3d"

    def read(self, mesh, **kwargs):
        f = open(self.filename)
        # read the whole file:
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
        while l == "" or l[0] == "#": # comment or an empty line
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
            if not row: break
            if len(row) == 0: continue

            if (row[0] == 'NUMNP'):
                row = fd.readline().split()
                n_nod, n_el, dim = row[0], row[1], int(row[4])
                break;

        if ret_fd:
            return dim, fd
        else:
            fd.close()
            return dim

    def read(self, mesh, **kwargs):

        el = {'3_8' : [], '3_4' : [], '2_4' : [], '2_3' : []}
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
            if not row: break
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
        mesh._set_io_data(nod, None, conns, mat_ids, descs, nodal_bcs=nodal_bcs)

        return mesh

    def write(self, filename, mesh, out=None, **kwargs):
        raise NotImplementedError


class XYZMeshIO(MeshIO):
    """
    Trivial XYZ format working only with coordinates (in a .XYZ file) and the
    connectivity stored in another file with the same base name and .IEN
    suffix.
    """
    format = 'xyz'

    def _read_coors(self):
        coors = nm.loadtxt(self.filename)
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

        if ret_fd: fd = open(self.filename, 'r')
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
                          dtype=nm.int32) - 1
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
                   conns[0] + 1)

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

def any_from_filename(filename, prefix_dir=None, file_format=None):
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
            # if io_class is 'meshio':
            kwargs['file_format'] = file_format
        else:
            raise ValueError('unknown file format! (%s)' % file_format)
    else:
        ext2io = {e: (v[0], k) for k, v in supported_formats.items()
            for e in v[1] if '*' not in v[2]}
        ext = op.splitext(filename)[1].lower()
        io_class = ext2io[ext][0]
        kwargs['file_format'] = ext2io[ext][1]

    return io_table[io_class](filename, **kwargs)


insert_static_method(MeshIO, any_from_filename)
del any_from_filename
