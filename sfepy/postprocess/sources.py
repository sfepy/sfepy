import os

import numpy as nm

try:
    from enthought.tvtk.api import tvtk
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.pyface.timer.api import Timer

except:
    from tvtk.api import tvtk
    from mayavi.sources.vtk_data_source import VTKDataSource
    from pyface.timer.api import Timer

from dataset_manager import DatasetManager

from sfepy.base.base import Struct, basestr
from sfepy.postprocess.utils import mlab
from sfepy.fem import Mesh
from sfepy.fem.meshio import MeshIO, vtk_cell_types, supported_formats

def create_file_source(filename, watch=False, offscreen=True):
    """Factory function to create a file source corresponding to the
    given file format."""
    kwargs = {'watch' : watch, 'offscreen' : offscreen}

    if isinstance(filename, basestr):
        fmt = os.path.splitext(filename)[1]
        is_sequence = False

    else: # A sequence.
        fmt = os.path.splitext(filename[0])[1]
        is_sequence = True

    fmt = fmt.lower()

    if fmt == '.vtk':
        # VTK is supported directly by Mayavi, no need to use MeshIO.
        if is_sequence:
            return VTKSequenceFileSource(filename, **kwargs)
        else:
            return VTKFileSource(filename, **kwargs)

    elif fmt in supported_formats.keys():
        if is_sequence:
            if fmt == '.h5':
                raise ValueError('format .h5 does not support file sequences!')
            else:
                return GenericSequenceFileSource(filename, **kwargs)
        else:
            return GenericFileSource(filename, **kwargs)

    else:
        raise ValueError('unknown file format! (%s)' % fmt)

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
        self.step_range = None
        self.notify_obj = None
        if self.watch:
            self.last_stat = os.stat(self.filename)
        self.set_step()

    def setup_mat_id(self, mat_id_name='mat_id', single_color=False):
        self.mat_id_name = mat_id_name
        self.single_color = single_color

    def set_step(self, step=0):
        """Set step of a data sequence."""
        self.step = step

    def get_step_range(self):
        return self.step_range

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

    def get_step_range(self):
        return (0, 0)

class VTKSequenceFileSource(VTKFileSource):
    """A thin wrapper around mlab.pipeline.open() for VTK file sequences."""

    def create_source(self):
        """Create a VTK file source """
        return mlab.pipeline.open(self.filename[0])

    def set_filename(self, filename, vis_source):
        self.filename = filename
        vis_source.base_file_name = filename[self.step]

    def get_step_range(self):
        return (0, len(self.filename) - 1)


class GenericFileSource(FileSource):
    """File source usable with any format supported by MeshIO classes."""

    def __init__(self, *args, **kwargs):
        FileSource.__init__(self, *args, **kwargs)

        self.io = None

    def read_common(self, filename):
        self.io = MeshIO.any_from_filename(filename)
        self.step_range = (0, self.io.read_last_step())

        self.mesh = mesh = Mesh.from_file(filename)
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
            mat_id = nm.concatenate(self.mesh.mat_ids)
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
        if self.dim == 2:
            bbox = nm.c_[bbox, [0.0, 0.0]]
        return bbox

    def set_filename(self, filename, vis_source):
        self.filename = filename
        self.source = self.create_source()
        vis_source.data = self.source.data
        
    def get_step_range(self):
        if self.step_range is None:
            io = MeshIO.any_from_filename(self.filename)
            self.step_range = (0, io.read_last_step())

        return self.step_range

    def get_mat_id(self, mat_id_name='mat_id'):
        """
        Get material ID numbers of the underlying mesh elements.
        """
        if self.source is not None:
            mat_id = nm.concatenate(self.mesh.mat_ids)
            return mat_id

    def file_changed(self):
        self.step_range = (0, self.io.read_last_step())

    def create_dataset(self):
        """Create a tvtk.UnstructuredGrid dataset from the Mesh instance of the
        file source."""
        mesh = self.mesh
        n_nod, dim = self.n_nod, self.dim
        n_el, n_els, n_e_ps = mesh.n_el, mesh.n_els, mesh.n_e_ps

        if dim == 2:
            nod_zz = nm.zeros((n_nod, 1), dtype=mesh.coors.dtype)
            points = nm.c_[mesh.coors, nod_zz]
        else:
            points = mesh.coors

        dataset = tvtk.UnstructuredGrid(points=points)

        cell_types = []
        cells = []
        offset = [0]
        for ig, conn in enumerate(mesh.conns):
            cell_types += [vtk_cell_types[mesh.descs[ig]]] * n_els[ig]

            nn = nm.array([conn.shape[1]] * n_els[ig])
            aux = nm.c_[nn[:,None], conn]
            cells.extend(aux.ravel())

            offset.extend([aux.shape[1]] * n_els[ig])

        cells = nm.array(cells)
        cell_types = nm.array(cell_types)
        offset = nm.cumsum(offset)[:-1]
        
        cell_array = tvtk.CellArray()
        cell_array.set_cells(n_el, cells)

        dataset.set_cells(cell_types, offset, cell_array)

        return dataset

    def add_data_to_dataset(self, dataset, data):
        """Add point and cell data to the dataset."""
        dim = self.dim
        sym = (dim + 1) * dim / 2

        dm = DatasetManager(dataset=dataset)
        for key, val in data.iteritems():
            vd = val.data
##             print vd.shape
            if val.mode == 'vertex':
                if vd.shape[1] == 1:
                    aux = vd.reshape((vd.shape[0],))

                elif vd.shape[1] == 2:
                    zz = nm.zeros((vd.shape[0], 1), dtype=vd.dtype)
                    aux = nm.c_[vd, zz]

                elif vd.shape[1] == 3:
                    aux = vd

                else:
                    raise ValueError('unknown vertex data format! (%s)'\
                                     % vd.shape)

                dm.add_array(aux, key, 'point')

            elif val.mode == 'cell':
                ne, aux, nr, nc = vd.shape
                if (nr == 1) and (nc == 1):
                    aux = vd.reshape((ne,))

                elif (nr == dim) and (nc == 1):
                    if dim == 3:
                        aux = vd.reshape((ne, dim))
                    else:
                        zz = nm.zeros((vd.shape[0], 1), dtype=vd.dtype);
                        aux = nm.c_[vd.squeeze(), zz]

                elif (((nr == sym) or (nr == (dim * dim))) and (nc == 1)) \
                         or ((nr == dim) and (nc == dim)):
                    vd = vd.squeeze()

                    if dim == 3:
                        if nr == sym:
                            aux = vd[:,[0,3,4,3,1,5,4,5,2]]
                        elif nr == (dim * dim):
                            aux = vd[:,[0,3,4,6,1,5,7,8,2]]
                        else:
                            aux = vd.reshape((vd.shape[0], dim*dim))
                    else:
                        zz = nm.zeros((vd.shape[0], 1), dtype=vd.dtype);
                        if nr == sym:
                            aux = nm.c_[vd[:,[0,2]], zz, vd[:,[2,1]],
                                        zz, zz, zz, zz]
                        elif nr == (dim * dim):
                            aux = nm.c_[vd[:,[0,2]], zz, vd[:,[3,1]],
                                        zz, zz, zz, zz]
                        else:
                            aux = nm.c_[vd[:,0,[0,1]], zz, vd[:,1,[0,1]],
                                        zz, zz, zz, zz]

                dm.add_array(aux, key, 'cell')

class GenericSequenceFileSource(GenericFileSource):
    """File source usable with any format supported by MeshIO classes, with
    exception of HDF5 (.h5), for file sequences."""

    def create_source(self):
        """Create a VTK source from data in a SfePy-supported file."""
        if self.io is None:
            self.read_common(self.filename[self.step])

        dataset = self.create_dataset()

        src = VTKDataSource(data=dataset)
        return src

    def set_filename(self, filename, vis_source):
        self.filename = filename
        self.io = None
        self.source = self.create_source()
        vis_source.data = self.source.data
        
    def get_step_range(self):
        return (0, len(self.filename) - 1)
