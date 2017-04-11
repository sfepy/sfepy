# coding=utf8
from __future__ import absolute_import
from sfepy import data_dir
import six

filename_meshes = ['/meshes/3d/cylinder.mesh',
                   '/meshes/3d/cylinder.vtk',
                   '/meshes/various_formats/small2d.mesh',
                   '/meshes/various_formats/small2d.vtk',
                   '/meshes/various_formats/octahedron.node',
                   '/meshes/various_formats/comsol_tri.txt',
                   '/meshes/various_formats/abaqus_hex.inp',
                   '/meshes/various_formats/abaqus_tet.inp',
                   '/meshes/various_formats/abaqus_quad.inp',
                   '/meshes/various_formats/abaqus_tri.inp',
                   '/meshes/various_formats/abaqus_quad_tri.inp',
                   '/meshes/various_formats/hex4.mesh3d',
                   '/meshes/various_formats/tetra8.mesh3d',
                   '/meshes/various_formats/cube.bdf',
                   '/meshes/various_formats/med_2d_tri_quad.med',
                   '/meshes/various_formats/med_3d_tet_hex.med',
                   '/meshes/various_formats/msh_tri.msh',
                   '/meshes/various_formats/msh_tetra.msh']
filename_meshes = [data_dir + name for name in filename_meshes]

def mesh_hook(mesh, mode):
    """
    Define a mesh programmatically.
    """
    if mode == 'read':
        nodes = [[0, 0], [1, 0], [1, 1], [0, 1]]
        nod_ids = [0, 0, 1, 1]
        conns = [[[0, 1, 2], [0, 2, 3]]]
        mat_ids = [[0, 1]]
        descs = ['2_3']

        mesh._set_io_data(nodes, nod_ids, conns, mat_ids, descs)

    elif mode == 'write':
        pass


from sfepy.discrete.fem.meshio import UserMeshIO

filename_meshes.extend([mesh_hook, UserMeshIO(mesh_hook)])

same = [(0, 1), (2, 3)]

import os.path as op
from sfepy.base.base import assert_, assert_equals
from sfepy.base.testing import TestCommon

class Test(TestCommon):
    """Write test names explicitely to impose a given order of evaluation."""
    tests = ['test_read_meshes', 'test_compare_same_meshes',
             'test_read_dimension', 'test_write_read_meshes',
             'test_hdf5_meshio'
            ]

    @staticmethod
    def from_conf(conf, options):
        return Test(conf=conf, options=options)

    def test_read_meshes(self):
        """Try to read all listed meshes."""
        from sfepy.discrete.fem import Mesh

        conf_dir = op.dirname(__file__)
        meshes = {}
        for ii, filename in enumerate(filename_meshes):
            self.report('%d. mesh: %s' % (ii + 1, filename))
            mesh = Mesh.from_file(filename, prefix_dir=conf_dir)

            assert_(mesh.dim == (mesh.coors.shape[1]))
            assert_(mesh.n_nod == (mesh.coors.shape[0]))
            assert_(mesh.n_nod == (mesh.cmesh.vertex_groups.shape[0]))
            assert_(mesh.n_el == mesh.cmesh.num[mesh.cmesh.tdim])

            self.report('read ok')
            meshes[filename] = mesh

        self.meshes = meshes

        return True

    def _compare_meshes(self, mesh0, mesh1):
        import numpy as nm

        oks = []

        ok0 = (mesh0.dim == mesh1.dim)
        if not ok0:
            self.report('dimension failed!')
        oks.append(ok0)

        ok0 = mesh0.n_nod == mesh1.n_nod
        if not ok0:
            self.report('number of nodes failed!')
        oks.append(ok0)

        ok0 = mesh0.n_el == mesh1.n_el
        if not ok0:
            self.report('number of elements failed!')
        oks.append(ok0)

        ok0 = mesh0.descs == mesh1.descs
        if not ok0:
            self.report('element types failed!')
        oks.append(ok0)

        ok0 = nm.allclose(mesh0.coors, mesh1.coors)
        if not ok0:
            self.report('nodes failed!')
        oks.append(ok0)

        ok0 = nm.all(mesh0.cmesh.vertex_groups == mesh1.cmesh.vertex_groups)
        if not ok0:
            self.report('node groups failed!')
        oks.append(ok0)

        ok0 = nm.all(mesh0.cmesh.cell_groups == mesh1.cmesh.cell_groups)
        if not ok0:
            self.report('material ids failed!')
        oks.append(ok0)

        ok0 = (nm.all(mesh0.cmesh.get_cell_conn().indices ==
                      mesh1.cmesh.get_cell_conn().indices) and
               nm.all(mesh0.cmesh.get_cell_conn().offsets ==
                      mesh1.cmesh.get_cell_conn().offsets))
        if not ok0:
            self.report('connectivities failed!')
        oks.append(ok0)

        return oks

    def test_compare_same_meshes(self):
        """
        Compare same meshes in various formats.
        """
        oks = []
        for i0, i1 in same:
            name0 = filename_meshes[i0]
            name1 = filename_meshes[i1]
            self.report('comparing meshes from "%s" and "%s"' % (name0, name1))

            mesh0 = self.meshes[name0]
            mesh1 = self.meshes[name1]
            oks = self._compare_meshes(mesh0, mesh1)

        return sum(oks) == len(oks)

    def test_read_dimension(self):
        from sfepy.discrete.fem import MeshIO

        meshes = {data_dir + '/meshes/various_formats/small2d.mesh' : 2,
                  data_dir + '/meshes/various_formats/small2d.vtk' : 2,
                  data_dir + '/meshes/various_formats/small3d.mesh' : 3}

        ok = True
        conf_dir = op.dirname(__file__)
        for filename, adim in six.iteritems(meshes):
            self.report('mesh: %s, dimension %d' % (filename, adim))
            io = MeshIO.any_from_filename(filename, prefix_dir=conf_dir)
            dim = io.read_dimension()
            if dim != adim:
                self.report('read dimension %d -> failed' % dim)
                ok = False
            else:
                self.report('read dimension %d -> ok' % dim)

        return ok

    def test_hdf5_meshio(self):
        from sfepy.discrete.fem import Mesh
        conf_dir = op.dirname(__file__)
        mesh0 = Mesh.from_file(data_dir +
                               '/meshes/various_formats/small3d.mesh',
                               prefix_dir=conf_dir)

        import numpy as np
        import scipy.sparse as ss
        import tempfile
        #import pickle
        from sfepy.discrete.fem.meshio import HDF5MeshIO
        from sfepy.base.base import Struct
        from sfepy.base.ioutils import Cached, Uncached
        from sfepy.discrete.iga.domain import IGDomain
        from sfepy.discrete.iga.domain_generators import gen_patch_block_domain

        shape = [4, 4, 4]
        dims = [5, 5, 5]
        centre = [0, 0, 0]
        degrees = [2, 2, 2]

        nurbs, bmesh, regions = gen_patch_block_domain(dims, shape, centre,
                                                       degrees,
                                                       cp_mode='greville',
                                                       name='iga')
        ig_domain = IGDomain('iga', nurbs, bmesh, regions=regions)

        data = {
            'list': range(4),
            'mesh1': mesh0,
            'mesh2': mesh0,
            'iga' : ig_domain,
            'mesh3': Uncached(mesh0),
            'cached1': Cached(1),
            'cached2': Cached(2),
            'types': ( True, False, None ),
            'tuple': ('first string', 'druhý UTF8 řetězec'),
            'struct': Struct(
                         double = np.arange(4, dtype = float),
                         int = np.array([2,3,4,7]),
                         sparse = ss.csr_matrix(np.array([1,0,0,5]).
                             reshape((2,2)))
                    )
        }

        with tempfile.NamedTemporaryFile(suffix = '.h5') as fil:
            io = HDF5MeshIO(fil.name)
            io.write(fil.name, mesh0, {
                'cdata' : Struct(
                    mode = 'custom',
                    data = data,
                    unpack_markers = True
                )
            })
            fout = io.read_data(0)
            out = fout['cdata']

            assert_(out['mesh1'] is out['mesh2'],
                'Two meshes should be in fact the same object')

            assert_(out['mesh3'] is not out['mesh2'],
                'Two meshes should be different object')

            assert_(self._compare_meshes(out['mesh1'], mesh0),
                'Failed to restore mesh')

            assert_(self._compare_meshes(out['mesh3'], mesh0),
                'Failed to restore mesh')

            assert_((out['struct'].sparse == data['struct'].sparse).todense()
                    .all(),'Sparse matrix restore failed')

        #this property is not restored
        del data['iga'].nurbs.nurbs

        #not supporting comparison
        del data['iga']._bnf
        del out['iga']._bnf

        #restoration of this property failed
        #del data['iga'].vertex_set_bcs
        #del out['iga'].vertex_set_bcs

        assert_equals( out, data )
        return True


    def test_write_read_meshes(self):
        """
        Try to write and then read all supported formats.
        """
        from sfepy.discrete.fem import Mesh
        from sfepy.discrete.fem.meshio import (supported_formats,
                                               supported_capabilities)

        conf_dir = op.dirname(__file__)
        mesh0 = Mesh.from_file(data_dir +
                               '/meshes/various_formats/small3d.mesh',
                               prefix_dir=conf_dir)

        oks = []
        for suffix, format_ in six.iteritems(supported_formats):
            if isinstance(format_, tuple):
                continue
            if 'w' not in supported_capabilities[format_]: continue

            filename = op.join(self.options.out_dir, 'test_mesh_wr' + suffix)
            self.report('%s format: %s' % (suffix, filename))

            mesh0.write(filename, io='auto')
            mesh1 = Mesh.from_file(filename)

            oks.extend(self._compare_meshes(mesh0, mesh1))

        return sum(oks) == len(oks)
