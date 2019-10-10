# coding=utf8
from __future__ import absolute_import
import os
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
                   '/meshes/various_formats/msh_tetra.msh',
                   '/meshes/various_formats/xyz_quad.xyz',
                   '/meshes/various_formats/xyz_tet.xyz']
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
from sfepy.base.base import assert_
from sfepy.base.testing import TestCommon

class Test(TestCommon):
    """Write test names explicitely to impose a given order of evaluation."""
    tests = ['test_read_meshes', 'test_compare_same_meshes',
             'test_read_dimension', 'test_write_read_meshes',
             'test_hdf5_meshio']

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

    def test_write_read_meshes(self):
        """
        Try to write and then read all supported formats.
        """
        from sfepy.discrete.fem import Mesh
        from sfepy.discrete.fem.meshio import (supported_formats,
                                               supported_capabilities)

        conf_dir = op.dirname(__file__)
        mesh0 = Mesh.from_file(data_dir
                               + '/meshes/various_formats/small3d.mesh',
                               prefix_dir=conf_dir)

        oks = []
        for suffix, format_ in six.iteritems(supported_formats):
            if isinstance(format_, tuple) or (format_ == 'xyz'):
                continue
            if 'w' not in supported_capabilities[format_]: continue

            filename = op.join(self.options.out_dir, 'test_mesh_wr' + suffix)
            self.report('%s format: %s' % (suffix, filename))

            mesh0.write(filename, io='auto')
            mesh1 = Mesh.from_file(filename)

            oks.extend(self._compare_meshes(mesh0, mesh1))

        return sum(oks) == len(oks)

    def test_hdf5_meshio(self):
        try:
            from igakit import igalib
        except ImportError:
            self.report('hdf5_meshio not-tested (missing igalib module)!')
            return True

        import tempfile
        import numpy as nm
        import scipy.sparse as sps
        from sfepy.discrete.fem.meshio import HDF5MeshIO
        from sfepy.base.base import Struct
        from sfepy.base.ioutils import Cached, Uncached, SoftLink, \
                                       DataSoftLink
        from sfepy.discrete.iga.domain import IGDomain
        from sfepy.discrete.iga.domain_generators import gen_patch_block_domain
        from sfepy.solvers.ts import TimeStepper
        from sfepy.discrete.fem import Mesh

        conf_dir = op.dirname(__file__)
        mesh0 = Mesh.from_file(data_dir +
                               '/meshes/various_formats/small3d.mesh',
                               prefix_dir=conf_dir)

        shape = [4, 4, 4]
        dims = [5, 5, 5]
        centre = [0, 0, 0]
        degrees = [2, 2, 2]

        nurbs, bmesh, regions = gen_patch_block_domain(dims, shape, centre,
                                                       degrees,
                                                       cp_mode='greville',
                                                       name='iga')
        ig_domain = IGDomain('iga', nurbs, bmesh, regions=regions)

        int_ar = nm.arange(4)

        data = {
            'list': range(4),
            'mesh1': mesh0,
            'mesh2': mesh0,
            'mesh3': Uncached(mesh0),
            'mesh4': SoftLink('/step0/__cdata/data/data/mesh2'),
            'mesh5': DataSoftLink('Mesh','/step0/__cdata/data/data/mesh1/data'),
            'mesh6': DataSoftLink('Mesh','/step0/__cdata/data/data/mesh2/data',
                mesh0),
            'mesh7': DataSoftLink('Mesh','/step0/__cdata/data/data/mesh1/data',
                True),
            'iga' : ig_domain,
            'cached1': Cached(1),
            'cached2': Cached(int_ar),
            'cached3': Cached(int_ar),
            'types': ( True, False, None ),
            'tuple': ('first string', 'druhý UTF8 řetězec'),
            'struct': Struct(
                double=nm.arange(4, dtype=float),
                int=nm.array([2,3,4,7]),
                sparse=sps.csr_matrix(nm.array([1,0,0,5]).
                                      reshape((2,2)))
             )
        }

        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as fil:
            io = HDF5MeshIO(fil.name)
            ts = TimeStepper(0,1.,0.1, 10)

            io.write(fil.name, mesh0, {
                'cdata' : Struct(
                    mode='custom',
                    data=data,
                    unpack_markers=False
                )
            }, ts=ts)
            ts.advance()

            mesh = io.read()
            data['problem_mesh'] = DataSoftLink('Mesh', '/mesh', mesh)

            io.write(fil.name, mesh0, {
                'cdata' : Struct(
                    mode='custom',
                    data=data,
                    unpack_markers=True
                )
            }, ts=ts)

            cache = {'/mesh': mesh }
            fout = io.read_data(0, cache=cache)
            fout2 = io.read_data(1, cache=cache )
            out = fout['cdata']
            out2 = fout2['cdata']

            assert_(out['mesh7'] is out2['mesh7'],
                'These two meshes should be in fact the same object')

            assert_(out['mesh6'] is out2['mesh6'],
                'These two meshes should be in fact the same object')

            assert_(out['mesh5'] is not out2['mesh5'],
                'These two meshes shouldn''t be in fact the same object')

            assert_(out['mesh1'] is out['mesh2'],
                'These two meshes should be in fact the same object')

            assert_(out['mesh1'] is out['mesh2'],
                'These two meshes should be in fact the same object')

            assert_(out['mesh4'] is out['mesh2'],
                'These two meshes should be in fact the same object')

            assert_(out['mesh5'] is not out['mesh2'],
                'These two meshes shouldn''t be in fact the same object')

            assert_(out['mesh6'] is out['mesh2'],
                'These two meshes should be in fact the same object')

            assert_(out['mesh7'] is not out['mesh2'],
                'These two meshes shouldn''t be in fact the same object')

            assert_(out['mesh3'] is not out['mesh2'],
                'These two meshes should be different objects')

            assert_(out['cached2'] is out['cached3'],
                'These two array should be the same object')

            assert_(out2['problem_mesh'] is mesh,
                'These two meshes should be the same objects')

            assert_(self._compare_meshes(out['mesh1'], mesh0),
                'Failed to restore mesh')

            assert_(self._compare_meshes(out['mesh3'], mesh0),
                'Failed to restore mesh')

            assert_((out['struct'].sparse == data['struct'].sparse).todense()
                    .all(), 'Sparse matrix restore failed')

            ts.advance()
            io.write(fil.name, mesh0, {
                    'cdata' : Struct(
                        mode='custom',
                        data=[
                            DataSoftLink('Mesh',
                                         '/step0/__cdata/data/data/mesh1/data',
                                         mesh0),
                            mesh0
                        ]
                    )
            }, ts=ts)
            out3 = io.read_data(2)['cdata']
            assert_(out3[0] is out3[1])

        os.remove(fil.name)

        #this property is not restored
        del data['iga'].nurbs.nurbs

        #not supporting comparison
        del data['iga']._bnf
        del out2['iga']._bnf

        #restoration of this property fails
        del data['iga'].vertex_set_bcs
        del out2['iga'].vertex_set_bcs

        #these soflink has no information how to unpack, so it must be
        #done manually
        data['mesh4'] = mesh0
        data['mesh5'] = mesh0
        data['mesh7'] = mesh0

        for key, val in six.iteritems(out2):
            self.report('comparing:', key)
            self.assert_equal(val, data[key])

        return True
