from sfepy import data_dir

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
                   '/meshes/various_formats/med_3d_tet_hex.med']
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

        mesh._set_data(nodes, nod_ids, conns, mat_ids, descs)

        ## mesh.write('aux.vtk', io='auto')

    elif mode == 'write':
        pass

from sfepy.fem.meshio import UserMeshIO

filename_meshes.extend([mesh_hook, UserMeshIO(mesh_hook)])

same = [(0, 1), (2, 3)]

import os.path as op
from sfepy.base.base import assert_
from sfepy.base.testing import TestCommon

##
# c: 05.02.2008
class Test( TestCommon ):
    """Write test names explicitely to impose a given order of evaluation."""
    tests = ['test_read_meshes', 'test_compare_same_meshes',
             'test_read_dimension', 'test_write_read_meshes']

    ##
    # c: 05.02.2008, r: 05.02.2008
    def from_conf( conf, options ):
        return Test( conf = conf, options = options )
    from_conf = staticmethod( from_conf )
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def test_read_meshes( self ):
        """Try to read all listed meshes."""
        from sfepy.fem import Mesh

        conf_dir = op.dirname(__file__)
        meshes = {}
        for ii, filename in enumerate( filename_meshes ):
            self.report( '%d. mesh: %s' % (ii + 1, filename) )
            mesh = Mesh.from_file(filename, prefix_dir=conf_dir)

            assert_(mesh.dim == (mesh.coors.shape[1]))
            assert_(mesh.n_nod == (mesh.coors.shape[0]))
            assert_(mesh.n_nod == (mesh.ngroups.shape[0]))
            assert_(mesh.n_el == sum(mesh.n_els))
            for ig, conn in enumerate( mesh.conns ):
                assert_(conn.shape[0] == len(mesh.mat_ids[ig]))
                assert_(conn.shape[0] == mesh.n_els[ig])
                assert_(conn.shape[1] == mesh.n_e_ps[ig])
                
            self.report( 'read ok' )
            meshes[filename] = mesh

        self.meshes = meshes

        return True

    def _compare_meshes(self, mesh0, mesh1):
        import numpy as nm

        oks = []

        ok0 = (mesh0.dim == mesh1.dim)
        if not ok0:
            self.report( 'dimension failed!' )
        oks.append( ok0 )

        ok0 = mesh0.n_nod == mesh1.n_nod
        if not ok0:
            self.report( 'number of nodes failed!' )
        oks.append( ok0 )

        ok0 = mesh0.n_el == mesh1.n_el
        if not ok0:
            self.report( 'number of elements failed!' )
        oks.append( ok0 )

        ok0 = mesh0.n_e_ps == mesh1.n_e_ps
        if not ok0:
            self.report( 'number of element points failed!' )
        oks.append( ok0 )

        ok0 = mesh0.descs == mesh1.descs
        if not ok0:
            self.report( 'element types failed!' )
        oks.append( ok0 )

        ok0 = nm.allclose( mesh0.coors, mesh1.coors )
        if not ok0:
            self.report( 'nodes failed!' )
        oks.append( ok0 )

        ok0 = nm.all( mesh0.ngroups == mesh1.ngroups )
        if not ok0:
            self.report( 'node groups failed!' )
        oks.append( ok0 )

        for ii in range( len( mesh0.mat_ids ) ):
            ok0 = nm.all( mesh0.mat_ids[ii] == mesh1.mat_ids[ii] )
            if not ok0:
                self.report( 'material ids failed!' )
            oks.append( ok0 )

        for ii in range( len( mesh0.mat_ids ) ):
            ok0 = nm.all( mesh0.conns[ii] == mesh1.conns[ii] )
            if not ok0:
                self.report( 'connectivities failed!' )
            oks.append( ok0 )

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

    ##
    # c: 03.07.2008, r: 03.07.2008
    def test_read_dimension( self ):
        from sfepy.fem import MeshIO

        meshes = {data_dir + '/meshes/various_formats/small2d.mesh' : 2,
                  data_dir + '/meshes/various_formats/small2d.vtk' : 2,
                  data_dir + '/meshes/various_formats/small3d.mesh' : 3}

        ok = True
        conf_dir = op.dirname(__file__)
        for filename, adim in meshes.iteritems():
            self.report( 'mesh: %s, dimension %d' % (filename, adim) )
            io = MeshIO.any_from_filename(filename, prefix_dir=conf_dir)
            dim = io.read_dimension()
            if dim != adim:
                self.report( 'read dimension %d -> failed' % dim )
                ok = False
            else:
                self.report( 'read dimension %d -> ok' % dim )

        return ok

    def test_write_read_meshes(self):
        """
        Try to write and then read all supported formats.
        """
        from sfepy.fem import Mesh
        from sfepy.fem.meshio import supported_formats, supported_capabilities

        conf_dir = op.dirname(__file__)
        mesh0 = Mesh.from_file(data_dir
                               + '/meshes/various_formats/small3d.mesh',
                               prefix_dir=conf_dir)

        oks = []
        for suffix, format_ in supported_formats.iteritems():
            if isinstance(format_, tuple):
                continue
            if 'w' not in supported_capabilities[format_]: continue

            filename = op.join(self.options.out_dir, 'test_mesh_wr' + suffix)
            self.report('%s format: %s' % (suffix, filename))

            mesh0.write(filename, io='auto')
            mesh1 = Mesh.from_file(filename)

            oks.extend(self._compare_meshes(mesh0, mesh1))

        return sum(oks) == len(oks)
