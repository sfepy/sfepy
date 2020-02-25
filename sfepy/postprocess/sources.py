from __future__ import absolute_import
import os

import numpy as nm
import six

try:
    from enthought.tvtk.api import tvtk
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.core.dataset_manager import DatasetManager
    from enthought.pyface.timer.api import Timer

except:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi.core.dataset_manager import DatasetManager
    from pyface.timer.api import Timer

from sfepy.base.base import Struct, basestr
from sfepy.postprocess.utils import mlab
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import MeshIO, supported_formats

vtk_cell_types = {'1_1' : 1, '1_2' : 3, '2_2' : 3, '3_2' : 3,
                  '2_3' : 5, '2_4' : 9, '3_4' : 10, '3_8' : 12}

def create_file_source(filename, watch=False, offscreen=True):
    """Factory function to create a file source corresponding to the
    given file format."""

    if isinstance(filename, FileSource):
        return filename

    from distutils.version import LooseVersion

    try:
        from enthought.mayavi import version

    except:
        from mayavi import version

    # Work around a Mayavi 4.4.x issue.
    can_vtk_source = LooseVersion(version.version) < LooseVersion('4.4')

    kwargs = {'watch' : watch, 'offscreen' : offscreen}

    if isinstance(filename, basestr):
        fmt = os.path.splitext(filename)[1]
        is_sequence = False

    else: # A sequence.
        fmt = os.path.splitext(filename[0])[1]
        is_sequence = True

    fmt = fmt.lower()

    if can_vtk_source and (fmt == '.vtk'):
        # VTK is supported directly by Mayavi, no need to use MeshIO.
        if is_sequence:
            return VTKSequenceFileSource(filename, **kwargs)
        else:
            return VTKFileSource(filename, **kwargs)

    else:
        if is_sequence:
            if fmt == '.h5':
                raise ValueError('format .h5 does not support file sequences!')
            else:
                return GenericSequenceFileSource(filename, **kwargs)
        else:
            return GenericFileSource(filename, **kwargs)

class FileSource(Struct):
    """General file source."""

    def __init__(self, filename, watch=False, offscreen=True):
        """Create a file source using the given file name."""
        mlab.options.offscreen = offscreen
        self.watch = watch
        self.filename = filename
        self.reset()

    def __call__(self, step=0):
        """Get the file source."""
        if self.source is None:
            self.source = self.create_source()
            if self.watch:
                self.timer = Timer(1000, self.poll_file)

        return self.source

    def reset(self):
        """Reset."""
        self.mat_id_name = None
        self.source = None
        self.notify_obj = None
        self.steps = []
        self.times = []
        self.step = 0
        self.time = 0.0
        if self.watch:
            self.last_stat = os.stat(self.filename)

    def setup_mat_id(self, mat_id_name='mat_id', single_color=False):
        self.mat_id_name = mat_id_name
        self.single_color = single_color

    def get_step_time(self, step=None, time=None):
        """
        Set current step and time to the values closest greater or equal to
        either step or time. Return the found values.
        """
        if (step is not None) and len(self.steps):
            step = step if step >= 0 else self.steps[-1] + step + 1
            ii = nm.searchsorted(self.steps, step)
            ii = nm.clip(ii, 0, len(self.steps) - 1)

            self.step = self.steps[ii]
            if len(self.times):
                self.time = self.times[ii]

        elif (time is not None) and len(self.times):
            ii = nm.searchsorted(self.times, time)
            ii = nm.clip(ii, 0, len(self.steps) - 1)

            self.step = self.steps[ii]
            self.time = self.times[ii]

        return self.step, self.time

    def get_ts_info(self):
        return self.steps, self.times

    def get_mat_id(self, mat_id_name='mat_id'):
        """
        Get material ID numbers of the underlying mesh elements.
        """
        if self.source is not None:
            dm = DatasetManager(dataset=self.source.outputs[0])

            mat_id = dm.cell_scalars[mat_id_name]
            return mat_id

    def file_changed(self):
        pass

    def setup_notification(self, obj, attr):
        """The attribute 'attr' of the object 'obj' will be set to True
        when the source file is watched and changes."""
        self.notify_obj = obj
        self.notify_attr = attr

    def poll_file(self):
        """Check the source file's time stamp and notify the
        self.notify_obj in case it changed. Subclasses should implement
        the file_changed() method."""
        if not self.notify_obj:
            return

        s = os.stat(self.filename)
        if s[-2] == self.last_stat[-2]:
            setattr(self.notify_obj, self.notify_attr, False)
        else:
            self.file_changed()
            setattr(self.notify_obj, self.notify_attr, True)
            self.last_stat = s

class VTKFileSource(FileSource):
    """A thin wrapper around mlab.pipeline.open()."""

    def create_source(self):
        """Create a VTK file source """
        return mlab.pipeline.open(self.filename)

    def get_bounding_box(self):
        bbox = nm.array(self.source.reader.unstructured_grid_output.bounds)
        return bbox.reshape((3,2)).T

    def set_filename(self, filename, vis_source):
        self.filename = filename
        vis_source.base_file_name = filename

        # Force re-read.
        vis_source.reader.modified()
        vis_source.update()
        # Propagate changes in the pipeline.
        vis_source.data_changed = True

class VTKSequenceFileSource(VTKFileSource):
    """A thin wrapper around mlab.pipeline.open() for VTK file sequences."""

    def __init__(self, *args, **kwargs):
        FileSource.__init__(self, *args, **kwargs)

        self.steps = nm.arange(len(self.filename), dtype=nm.int32)

    def create_source(self):
        """Create a VTK file source """
        return mlab.pipeline.open(self.filename[0])

    def set_filename(self, filename, vis_source):
        self.filename = filename
        vis_source.base_file_name = filename[self.step]

class GenericFileSource(FileSource):
    """File source usable with any format supported by MeshIO classes."""

    def __init__(self, *args, **kwargs):
        FileSource.__init__(self, *args, **kwargs)

        self.io = None
        self.read_common(self.filename)

    def read_common(self, filename):
        self.io = MeshIO.any_from_filename(filename)
        self.steps, self.times, _ = self.io.read_times()

        self.mesh = Mesh.from_file(filename)
        self.n_nod, self.dim = self.mesh.coors.shape

    def create_source(self):
        """
        Create a VTK source from data in a SfePy-supported file.

        Notes
        -----
        All data need to be set here, otherwise time stepping will not
        work properly - data added by user later will be thrown away on
        time step change.
        """
        if self.io is None:
            self.read_common(self.filename)

        dataset = self.create_dataset()

        try:
            out = self.io.read_data(self.step)
        except ValueError:
            out = None

        if out is not None:
            self.add_data_to_dataset(dataset, out)

        if self.mat_id_name is not None:
            mat_id = self.mesh.cmesh.cell_groups
            if self.single_color:
                rm = mat_id.min(), mat_id.max()
                mat_id[mat_id > rm[0]] = rm[1]

            dm = DatasetManager(dataset=dataset)
            dm.add_array(mat_id, self.mat_id_name, 'cell')

        src = VTKDataSource(data=dataset)
#        src.print_traits()
#        debug()
        return src

    def get_bounding_box(self):
        bbox = self.mesh.get_bounding_box()
        if self.dim < 3:
            bbox = nm.c_[bbox, nm.zeros((2, 3 - self.dim))]
        return bbox

    def set_filename(self, filename, vis_source):
        self.filename = filename
        self.source = self.create_source()
        vis_source.data = self.source.data

    def get_mat_id(self, mat_id_name='mat_id'):
        """
        Get material ID numbers of the underlying mesh elements.
        """
        if self.source is not None:
            mat_id = self.mesh.cmesh.cell_groups
            return mat_id

    def file_changed(self):
        self.steps, self.times, _ = self.io.read_times()

    def create_dataset(self):
        """Create a tvtk.UnstructuredGrid dataset from the Mesh instance of the
        file source."""
        mesh = self.mesh
        n_nod, dim = self.n_nod, self.dim

        if dim < 3:
            nod_zz = nm.zeros((n_nod, 3 - dim), dtype=mesh.coors.dtype)
            points = nm.c_[mesh.coors, nod_zz]
        else:
            points = mesh.coors

        dataset = tvtk.UnstructuredGrid(points=points)

        cell_types = []
        cells = []
        offset = [0]
        n_cells = [1]
        for ig, desc in enumerate(mesh.descs):
            conn = mesh.get_conn(desc)
            n_cell, n_cv = conn.shape

            n_cells.append(n_cell)

            aux = nm.empty(n_cell, dtype=nm.int32)
            aux.fill(vtk_cell_types[desc])
            cell_types.append(aux)

            aux = nm.empty((n_cell, n_cv + 1), dtype=nm.int32)
            aux[:, 0] = n_cv
            aux[:, 1:] = conn
            cells.append(aux.ravel())

            offset.append(aux.shape[1])

        cells = nm.concatenate(cells)
        cell_types = nm.concatenate(cell_types)

        offset = nm.repeat(offset, n_cells)
        offset = nm.cumsum(offset)[:-1]

        cell_array = tvtk.CellArray()
        cell_array.set_cells(mesh.n_el, cells)

        dataset.set_cells(cell_types, offset, cell_array)

        return dataset

    def add_data_to_dataset(self, dataset, data):
        """Add point and cell data to the dataset."""
        dm = DatasetManager(dataset=dataset)
        for key, val in six.iteritems(data):
            if val.mode == 'vertex':
                aux = self._reshape(val.data, self.dim)
                dm.add_array(aux, key, 'point')

            elif val.mode == 'cell':
                aux = self._reshape(val.data, self.dim)
                dm.add_array(aux, key, 'cell')

    @staticmethod
    def _reshape(data, dim):
        sym = (dim + 1) * dim // 2

        if len(data.shape) == 1:
            num, nr, nc = data.shape[0], 1, 1

        elif len(data.shape) == 2:
            num, nr, nc = data.shape[0], data.shape[1], 1

        else:
            num, _, nr, nc = data.shape

        if (nr == 1) and (nc == 1):
            aux = data.reshape((num,))

        elif (nr == dim) and (nc == 1):
            if dim == 3:
                aux = data.reshape((num, dim))
            else:
                zz = nm.zeros((num, 1), dtype=data.dtype);
                aux = nm.c_[data.squeeze(), zz]

        elif ((((nr == sym) or (nr == (dim * dim))) and (nc == 1))
              or ((nr == dim) and (nc == dim))):
            data = data.squeeze()

            if dim == 3:
                if nr == sym:
                    aux = data[:,[0,3,4,3,1,5,4,5,2]]
                elif nr == (dim * dim):
                    aux = data[:,[0,3,4,6,1,5,7,8,2]]
                else:
                    aux = data.reshape((num, dim*dim))
            else:
                zz = nm.zeros((num, 1), dtype=data.dtype);
                if nr == sym:
                    aux = nm.c_[data[:,[0,2]], zz, data[:,[2,1]],
                                zz, zz, zz, zz]
                elif nr == (dim * dim):
                    aux = nm.c_[data[:,[0,2]], zz, data[:,[3,1]],
                                zz, zz, zz, zz]
                else:
                    aux = nm.c_[data[:,0,[0,1]], zz, data[:,1,[0,1]],
                                zz, zz, zz, zz]
        else:
            raise ValueError('unknown data shape! (%s)' % list(data.shape))

        return aux

class GenericSequenceFileSource(GenericFileSource):
    """File source usable with any format supported by MeshIO classes, with
    exception of HDF5 (.h5), for file sequences."""

    def read_common(self, filename):
        GenericFileSource.read_common(self, filename[self.step])

        self.steps = nm.arange(len(self.filename), dtype=nm.int32)
        self.times = nm.zeros(len(self.filename), dtype=nm.float64)

    def set_filename(self, filename, vis_source):
        self.filename = filename
        self.io = None
        self.source = self.create_source()
        vis_source.data = self.source.data

    def file_changed(self):
        pass
