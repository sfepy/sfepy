from sfepy.base.base import *
from sfepy.postprocess.utils import mlab
from sfepy.fem import Mesh
from sfepy.fem.meshio import HDF5MeshIO, vtk_cell_types
from sfepy.solvers.ts import TimeStepper

from dataset_manager import DatasetManager
from enthought.tvtk.api import tvtk
from enthought.mayavi.sources.vtk_data_source import VTKDataSource
from enthought.pyface.timer.api import Timer

def create_file_source(filename, watch=False, offscreen=True):
    """Factory function to create a file source corresponding to the
    given file format."""
    kwargs = {'watch' : watch, 'offscreen' : offscreen}

    if isinstance(filename, str):
        fmt = os.path.splitext(filename)[1]
        is_sequence = False

    else: # A sequence.
        fmt = os.path.splitext(filename[0])[1]
        is_sequence = True

    if fmt.lower() == '.vtk':
        if is_sequence:
            return VTKSequenceFileSource(filename, **kwargs)
        else:
            return VTKFileSource(filename, **kwargs)

    elif fmt.lower() == '.h5':
        return HDF5FileSource(filename, **kwargs)

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
        self.source = None
        self.step_range = None
        self.notify_obj = None
        if self.watch:
            self.last_stat = os.stat(self.filename)
        self.set_step()

    def set_step(self, step=0):
        """Set step of a data sequence."""
        self.step = step

    def get_step_range(self):
        return self.step_range

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

    def create_source(self):
        """Create a VTK file source """
        return mlab.pipeline.open(self.filename)

    def get_bounding_box(self):
        bbox = nm.array(self.source.reader.unstructured_grid_output.bounds)
        return bbox.reshape((3,2)).T

    def set_filename(self, filename, vis_source):
        self.filename = filename
        vis_source.base_file_name = filename

    def get_step_range(self):
        return (0, 0)

class VTKSequenceFileSource(VTKFileSource):

    def create_source(self):
        """Create a VTK file source """
        return mlab.pipeline.open(self.filename[0])

    def set_filename(self, filename, vis_source):
        self.filename = filename
        vis_source.base_file_name = filename[self.step]

    def get_step_range(self):
        return (0, len(self.filename) - 1)

class HDF5FileSource(FileSource):

    def __init__(self, *args, **kwargs):
        FileSource.__init__(self, *args, **kwargs)

        self.io = None

    def read_common(self, filename):
        self.io = HDF5MeshIO(filename)
        self.ts = TimeStepper(*self.io.read_time_stepper())

        self.step_range = (0, self.io.read_last_step())

        self.mesh = mesh = Mesh.from_file(self.filename)
        self.n_nod, self.dim = self.mesh.coors.shape
        
    def create_source(self):
        """Create a VTK source from data in a SfePy HDF5 file."""
        if self.io is None:
            self.read_common(self.filename)

        mesh = self.mesh
        n_nod, dim = self.n_nod, self.dim
        sym = (dim + 1) * dim / 2
        n_el, n_els, n_e_ps = mesh.n_el, mesh.n_els, mesh.n_e_ps

        out = self.io.read_data(self.step)

        if dim == 2:
            nod_zz = nm.zeros((n_nod, 1), dtype=mesh.coors.dtype)
            points = nm.c_[mesh.coors, nod_zz]
        else:
            points = mesh.coors

        data = tvtk.UnstructuredGrid(points=points)
        dm = DatasetManager(dataset=data)

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

        data.set_cells(cell_types, offset, cell_array)

        for key, val in out.iteritems():
            vd = val.data
##             print vd.shape
            if val.mode == 'vertex':
                if vd.shape[1] == 1:
                    aux = vd.reshape((vd.shape[0],))

                elif vd.shape[1] == 2:
                    aux = nm.c_[vd, nod_zz]

                elif vd.shape[1] == 3:
                    aux = vd

                else:
                    raise ValueError('unknown vertex data format! (%s)'\
                                     % vd.shape)

                dm.add_array(aux, key, 'point')

            elif val.mode == 'cell':
                cdata = data.cell_data
                ne, aux, nr, nc = val.data.shape
                if (nr == 1) and (nc == 1):
                    aux = vd.reshape((ne,))

                elif (nr == dim) and (nc == 1):
                    if dim == 3:
                        aux = vd.reshape((ne, dim))
                    else:
                        zz = nm.zeros( (vd.shape[0], 1), dtype = nm.float64 );
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
                        zz = nm.zeros( (vd.shape[0], 1), dtype = nm.float64 );
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
        
        src = VTKDataSource(data=data)
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
            io = HDF5MeshIO(self.filename)
            self.step_range = (0, io.read_last_step())

        return self.step_range

    def file_changed(self):
        self.step_range = (0, self.io.read_last_step())
