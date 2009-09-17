from sfepy.base.base import *
from sfepy.postprocess.utils import mlab
from sfepy.fem import Mesh
from sfepy.fem.meshio import HDF5MeshIO, vtk_cell_types
from sfepy.solvers.ts import TimeStepper

from dataset_manager import DatasetManager
from enthought.tvtk.api import tvtk
from enthought.mayavi.sources.vtk_data_source import VTKDataSource

def create_file_source(filename):
    """Factory function to create a file source corresponding to the
    given file format."""
    fmt = os.path.splitext(filename)[1]

    if fmt.lower() == '.vtk':
        return VTKFileSource(filename)

    elif fmt.lower() == '.h5':
        return HDF5FileSource(filename)

    else:
        raise ValueError('unknown file format! (%s)' % fmt)

class FileSource(Struct):
    
    """General file source."""
    def __init__(self, filename):
        """Create a file source using the given file name."""
        self.filename = filename
        self.reset()

    def __call__(self, step=0):
        """Get the file source."""
        if self.source is None:
            self.source = self.create_source()

        return self.source

    def reset(self):
        """Reset."""
        self.source = None
        self.set_step()

    def set_step(self, step=0):
        """Set step of a data sequence."""
        self.step = step

class VTKFileSource(FileSource):

    def create_source(self):
        """Create a VTK file source """
        return mlab.pipeline.open(self.filename)

    def get_bounding_box(self):
        bbox = nm.array(self.source.reader.unstructured_grid_output.bounds)
        return bbox.reshape((3,2)).T

    def set_source_filename(self, filename):
        self.filename = filename
        self.source.base_file_name = filename

class HDF5FileSource(FileSource):

    def create_source(self):
        """Create a VTK source from data in a SfePy HDF5 file."""
        self.io = HDF5MeshIO(self.filename)
        self.ts = TimeStepper(*self.io.read_time_stepper())

        self.mesh = mesh = Mesh.from_file(self.filename)
        n_el, n_els, n_e_ps = mesh.n_el, mesh.n_els, mesh.n_e_ps

        out = self.io.read_data(self.step)

        data = tvtk.UnstructuredGrid(points=mesh.coors)
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

                elif vd.shape[1] == 3:
                    aux = vd

                dm.add_array(vd, key, 'point')

            elif val.mode == 'cell':
                cdata = data.cell_data
                dim, sym = 3, 6
                ne, aux, nr, nc = val.data.shape
                if (nr == 1) and (nc == 1):
                    aux = vd.reshape((ne,))

                elif (nr == dim) and (nc == 1):
                    aux = vd.reshape((ne, dim))

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
        return bbox

    def set_source_filename(self, filename):
        self.filename = filename
        self.source = self.create_source()
        
